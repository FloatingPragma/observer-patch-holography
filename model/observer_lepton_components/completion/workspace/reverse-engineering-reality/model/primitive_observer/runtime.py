"""Finite complete CQ execution, independent of the working field.

An instrument groups Kraus maps by visible outcome. Hidden alternatives are
summed as CP maps and never become feedback symbols. A Schur representation
is available for diagonal hidden channels. Its work basis is lifted over an
untouched reference, which is always the minor tensor factor.

This runtime checks finite numerical contracts. It does not turn source
bindings into a derivation of apparatus, action, clock or physical services.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
from itertools import product
import json
import math

import numpy as np

from tools import strict_json
from ..common_source.open_protocol import _density, conjugate


TOLERANCE = 2e-9
MAX_EVENTS = 1024
MAX_APPLICATIONS = 1_000_000
MAX_BRANCHES = 256
MAX_DENSITY_ENTRIES = 2_000_000
MAX_CHECKPOINT_BYTES = 16_000_000
MAX_RESOURCE = 10**18


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def _binding(value):
    if type(value) is str:
        value = strict_json.loads(value)
    if type(value) is not dict:
        raise ValueError('a JSON object binding is required')
    try:
        return _json(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError('finite JSON source bindings required') from error


def _immutable(value):
    value = np.ascontiguousarray(value, dtype=complex)
    return np.frombuffer(value.tobytes(), dtype=complex).reshape(value.shape)


def _matrix(value, size=None):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked quantum matrices are inadmissible')
    # Validate before coercion, including mixed lists containing booleans.
    def original(item):
        if isinstance(item, (bool, np.bool_)):
            raise ValueError('Boolean quantum entries are inadmissible')
        if isinstance(item, (list, tuple)):
            for child in item:
                original(child)
    original(value)
    raw = np.asarray(value)
    if (raw.dtype.kind not in 'iufc' or raw.ndim != 2 or not raw.shape[0]
            or raw.shape[0] != raw.shape[1] or raw.size > MAX_DENSITY_ENTRIES
            or (size is not None and raw.shape != (size, size))
            or not np.isfinite(raw).all()):
        raise ValueError('bounded finite square quantum matrix required')
    with np.errstate(over='ignore', invalid='ignore'):
        result = _immutable(raw)
    if not np.isfinite(result).all():
        raise ValueError('quantum matrix must remain finite in binary64')
    return result


def _cost(value):
    if (type(value) is not tuple or any(type(row) is not tuple or len(row) != 2
                                      for row in value)):
        raise ValueError('resources must be a tuple of (name, count) pairs')
    names = tuple(row[0] for row in value)
    if (any(type(name) is not str or not name for name in names)
            or len(set(names)) != len(names)
            or any(type(count) is not int or not 0 <= count <= MAX_RESOURCE
                   for _, count in value)):
        raise ValueError('unique resource names and bounded nonnegative integer counts required')
    return tuple(sorted(value))


def add_resources(*values):
    totals = {}
    for value in values:
        for name, count in value:
            totals[name] = totals.get(name, 0)+count
    return _cost(tuple(totals.items()))


def _scaled(value, count):
    return _cost(tuple((name, amount*count) for name, amount in value))


def _array_digest(value):
    return hashlib.sha256(_json(list(value.shape)).encode()+value.tobytes()).hexdigest()


@dataclass(frozen=True, eq=False)
class Instrument:
    """Complete visible CP groups, or one hidden channel on the work algebra."""

    kraus: tuple[tuple[np.ndarray, ...], ...] = ()
    visible: bool = False
    schur: np.ndarray | None = None
    resources: tuple[tuple[str, int], ...] = ()
    binding: object = None

    def __post_init__(self):
        if type(self.visible) is not bool or type(self.kraus) is not tuple:
            raise ValueError('typed visibility and immutable Kraus groups required')
        resources = _cost(self.resources)
        binding = _binding({} if self.binding is None else self.binding)
        if self.schur is not None:
            if self.kraus or self.visible:
                raise ValueError('Schur channels are hidden and cannot also supply Kraus groups')
            schur = _matrix(self.schur)
            if (np.linalg.norm(schur-schur.conj().T, 2) > TOLERANCE
                    or np.max(np.abs(np.diag(schur)-1)) > TOLERANCE
                    or np.linalg.eigvalsh((schur+schur.conj().T)/2)[0] < -TOLERANCE):
                raise ValueError('a Schur channel requires a positive correlation matrix')
            object.__setattr__(self, 'schur', schur)
        else:
            if (not self.kraus or len(self.kraus) > MAX_BRANCHES
                    or (not self.visible and len(self.kraus) != 1)
                    or any(type(group) is not tuple or not group or len(group) > 256
                           for group in self.kraus)):
                raise ValueError('every visible outcome needs its nonempty complete CP group')
            first = _matrix(self.kraus[0][0])
            size = first.shape[0]
            groups = tuple(tuple(_matrix(k, size) for k in group) for group in self.kraus)
            total = sum((k.conj().T@k for group in groups for k in group),
                        np.zeros((size, size), complex))
            if not np.isfinite(total).all() or np.linalg.norm(total-np.eye(size), 2) > TOLERANCE:
                raise ValueError('the full instrument must be trace preserving')
            object.__setattr__(self, 'kraus', groups)
        object.__setattr__(self, 'resources', resources)
        object.__setattr__(self, 'binding', binding)

    @property
    def dimension(self):
        return self.schur.shape[0] if self.schur is not None else self.kraus[0][0].shape[0]

    @property
    def outcomes(self):
        return len(self.kraus) if self.visible else 1

    def specification(self):
        return {'visible': self.visible, 'resources': self.resources,
                'binding': strict_json.loads(self.binding),
                'kraus': [[_array_digest(k) for k in group] for group in self.kraus],
                'schur': None if self.schur is None else _array_digest(self.schur)}


@dataclass(frozen=True)
class Event:
    """A single visible/controlled instrument or a repeated hidden sequence."""

    name: str
    instruments: tuple[Instrument, ...] = ()
    control: int | None = None
    sequence: tuple[Instrument, ...] = ()
    repeat: int = 1

    def __post_init__(self):
        if type(self.name) is not str or not self.name:
            raise ValueError('a nonempty event name required')
        if (type(self.instruments) is not tuple or type(self.sequence) is not tuple
                or bool(self.instruments) == bool(self.sequence)
                or any(type(i) is not Instrument for i in (*self.instruments, *self.sequence))):
            raise ValueError('choose one immutable instrument menu or hidden sequence')
        if type(self.repeat) is not int or not 1 <= self.repeat <= MAX_APPLICATIONS:
            raise ValueError('bounded positive integer repetition count required')
        if self.sequence:
            if self.control is not None or any(i.visible for i in self.sequence):
                raise ValueError('a repeated sequence contains only hidden uncontrolled channels')
        else:
            if (self.repeat != 1 or (self.control is None and len(self.instruments) != 1)
                    or (self.control is not None and
                        (type(self.control) is not int or self.control < 0))):
                raise ValueError('controlled alternatives require an earlier visible record')
            first = self.instruments[0]
            if any((i.visible, i.outcomes) != (first.visible, first.outcomes) for i in self.instruments):
                raise ValueError('feedback alternatives must retain the same visible outcome menu')

    @property
    def visible(self):
        return bool(self.instruments and self.instruments[0].visible)

    @property
    def outcomes(self):
        return self.instruments[0].outcomes if self.instruments else 1

    def specification(self):
        return {'name': self.name, 'control': self.control, 'repeat': self.repeat,
                'instruments': [i.specification() for i in self.instruments],
                'sequence': [i.specification() for i in self.sequence]}


@dataclass(frozen=True, eq=False)
class Program:
    name: str
    work_dimension: int
    initial_density: np.ndarray
    events: tuple[Event, ...]
    preparation_resources: tuple[tuple[str, int], ...] = ()
    binding: object = None

    def __post_init__(self):
        if (type(self.name) is not str or not self.name or type(self.work_dimension) is not int
                or self.work_dimension < 1 or self.work_dimension**2 > MAX_DENSITY_ENTRIES
                or type(self.events) is not tuple or len(self.events) > MAX_EVENTS):
            raise ValueError('bounded named finite program and immutable event sequence required')
        initial = _density(self.initial_density, self.work_dimension, normalized=True)
        if initial.shape != (self.work_dimension,)*2:
            raise ValueError('default preparation must be on work alone')
        preparation_resources = _cost(self.preparation_resources)
        menus, applications, maximum_cost = [], 0, preparation_resources
        for event in self.events:
            if type(event) is not Event:
                raise ValueError('each declared event must be typed')
            if event.control is not None:
                if event.control >= len(menus) or len(event.instruments) != menus[event.control]:
                    raise ValueError('feedback can use only a complete earlier outcome menu')
            if any(i.dimension != self.work_dimension for i in (*event.instruments, *event.sequence)):
                raise ValueError('every operation acts on the same work algebra')
            if event.visible:
                menus.append(event.outcomes)
            applications += event.repeat*len(event.sequence) if event.sequence else 1
            if event.sequence:
                cost = _scaled(add_resources(*(i.resources for i in event.sequence)), event.repeat)
            else:
                keys = set(name for instrument in event.instruments for name, _ in instrument.resources)
                cost = tuple((name, max(dict(i.resources).get(name, 0) for i in event.instruments))
                             for name in sorted(keys))
            maximum_cost = add_resources(maximum_cost, cost)
        if math.prod(menus) > MAX_BRANCHES or applications > MAX_APPLICATIONS:
            raise ValueError('complete branch or application budget exceeded')
        object.__setattr__(self, 'initial_density', _immutable(initial))
        object.__setattr__(self, 'preparation_resources', preparation_resources)
        object.__setattr__(self, 'binding', _binding({} if self.binding is None else self.binding))

    def specification(self):
        # A repeated apparatus shares a small set of compiled matrices. Bind
        # each one once, then retain every event address and repeat count.
        indices, instruments, events = {}, [], []
        def address(instrument):
            key = id(instrument)
            if key not in indices:
                indices[key] = len(instruments)
                instruments.append(instrument.specification())
            return indices[key]
        for event in self.events:
            events.append({'name': event.name, 'control': event.control, 'repeat': event.repeat,
                           'instruments': [address(i) for i in event.instruments],
                           'sequence': [address(i) for i in event.sequence]})
        return {'schema': 'oph.primitive_observer.program.v1', 'name': self.name,
                'work_dimension': self.work_dimension,
                'initial_density': _array_digest(self.initial_density),
                'preparation_resources': self.preparation_resources,
                'binding': strict_json.loads(self.binding),
                'instruments': instruments, 'events': events}

    @property
    def digest(self):
        return hashlib.sha256(_json(self.specification()).encode()).hexdigest()


def validate_program(program):
    """Detach caller-owned frozen envelopes; validate each shared map once."""
    if type(program) is not Program:
        raise ValueError('a typed finite program required')
    cache = {}
    def instrument(value):
        if type(value) is not Instrument:
            raise ValueError('typed instruments required')
        key = id(value)
        if key not in cache:
            cache[key] = Instrument(value.kraus, value.visible, value.schur,
                                    value.resources, value.binding)
        return cache[key]
    if type(program.events) is not tuple or len(program.events) > MAX_EVENTS:
        raise ValueError('bounded immutable event sequence required')
    events = []
    for event in program.events:
        if type(event) is not Event:
            raise ValueError('typed events required')
        events.append(Event(event.name, tuple(instrument(i) for i in event.instruments),
                            event.control, tuple(instrument(i) for i in event.sequence), event.repeat))
    return Program(program.name, program.work_dimension, program.initial_density, tuple(events),
                   program.preparation_resources, program.binding)


@dataclass(frozen=True)
class Branch:
    record: tuple[int, ...]
    density: np.ndarray
    resources: tuple[tuple[str, int], ...]

    @property
    def probability(self):
        return float(np.trace(self.density).real)


@dataclass(frozen=True)
class Stage:
    cursor: int
    name: str
    branches: int
    probability_mass: float


@dataclass(frozen=True)
class History:
    program: Program
    cursor: int
    reference_dimension: int
    initial_density: np.ndarray
    preparation_mode: str
    branches: tuple[Branch, ...]
    expanded_tree_resources: tuple[tuple[str, int], ...]
    stages: tuple[Stage, ...]
    seal: str = ''

    @property
    def work_dimension(self):
        return self.program.work_dimension

    @property
    def probability_mass(self):
        return math.fsum(branch.probability for branch in self.branches)

    @property
    def records(self):
        return tuple(branch.record for branch in self.branches)

    @property
    def resources(self):
        return self.expanded_tree_resources


def _records(program, cursor):
    menus = [range(e.outcomes) for e in program.events[:cursor] if e.visible]
    return tuple(product(*menus))


def _header(history):
    return {'program_digest': history.program.digest, 'cursor': history.cursor,
            'reference_dimension': history.reference_dimension,
            'preparation_mode': history.preparation_mode,
            'expanded_tree_resources': history.expanded_tree_resources,
            'stages': [stage.__dict__ for stage in history.stages],
            'branches': [{'record': b.record, 'resources': b.resources} for b in history.branches]}


def _seal(history):
    return hashlib.sha256(_json(_header(history)).encode()+_immutable(history.initial_density).tobytes()
                          +b''.join(_immutable(b.density).tobytes() for b in history.branches)).hexdigest()


def validate_history(history):
    """Validate integrity/domain of an in-memory result, without replay.

    A public digest is not authentication against a caller who fabricates an
    entire object and reseals it. Use checkpoint_from_json for untrusted saved
    data; that boundary independently reruns the declared prefix.
    """
    if type(history) is not History:
        raise ValueError('a typed complete history required')
    program = validate_program(history.program)
    if (type(history.cursor) is not int or not 0 <= history.cursor <= len(program.events)
            or type(history.reference_dimension) is not int or history.reference_dimension < 1
            or history.preparation_mode not in ('declared', 'supplied')
            or type(history.branches) is not tuple or type(history.stages) is not tuple):
        raise ValueError('invalid history domain')
    records = _records(program, history.cursor)
    size = program.work_dimension*history.reference_dimension
    if len(records)*size**2 > MAX_DENSITY_ENTRIES or len(history.branches) != len(records):
        raise ValueError('complete history allocation budget exceeded')
    initial = _density(history.initial_density, program.work_dimension, normalized=True)
    if initial.shape != (size, size):
        raise ValueError('initial density has a different reference dimension')
    if history.preparation_mode == 'declared' and (
            history.reference_dimension != 1 or not np.array_equal(initial, program.initial_density)):
        raise ValueError('declared preparation differs from program')
    branches = []
    for branch, record in zip(history.branches, records, strict=True):
        if (type(branch) is not Branch or type(branch.record) is not tuple
                or any(type(x) is not int for x in branch.record) or branch.record != record
                or not isinstance(branch.density, np.ndarray) or branch.density.shape != (size, size)):
            raise ValueError('every ordered complete record and density block required')
        density = _density(branch.density, program.work_dimension, normalized=False)
        branches.append(Branch(record, _immutable(density), _cost(branch.resources)))
    if not math.isfinite(history.probability_mass) or abs(history.probability_mass-1) > TOLERANCE:
        raise ValueError('complete history lost probability mass')
    if len(history.stages) != history.cursor:
        raise ValueError('complete stage journal required')
    for cursor, stage in enumerate(history.stages, 1):
        if (type(stage) is not Stage or type(stage.cursor) is not int or stage.cursor != cursor
                or stage.name != program.events[cursor-1].name
                or type(stage.branches) is not int or stage.branches != len(_records(program, cursor))
                or type(stage.probability_mass) not in (int, float)
                or not math.isfinite(stage.probability_mass)
                or abs(stage.probability_mass-1) > TOLERANCE):
            raise ValueError('stage journal differs from the declared complete prefix')
    _cost(history.expanded_tree_resources)
    if not history.seal or _seal(history) != history.seal:
        raise ValueError('checkpoint integrity digest changed')
    return replace(history, program=program, initial_density=_immutable(initial), branches=tuple(branches),
                   stages=tuple(replace(s) for s in history.stages))


def _apply(instrument, rho, reference):
    if instrument.schur is not None:
        size = instrument.dimension
        tensor = rho.reshape(size, reference, size, reference)
        return ((tensor*instrument.schur[:, None, :, None]).reshape(rho.shape),)
    def apply(k):
        # This is the same work action. Avoid a tensor contraction planner on
        # every empirical microstep when there is no spectator reference.
        return k@rho@k.conj().T if reference == 1 else conjugate(k, rho, reference)
    return tuple(apply(group[0]) if len(group) == 1 else
                 sum((apply(k) for k in group), np.zeros_like(rho))
                 for group in instrument.kraus)


def execute(program, *, initial=None, stop_after=None, checkpoint=None, max_branches=MAX_BRANCHES):
    """Execute full unnormalized records; repetitions never reset the state.

    In-memory checkpoints are trusted results of execute or a replay-validated
    import. Their integrity is checked; importing untrusted saved data uses
    checkpoint_from_json, which additionally reconstructs the entire prefix.
    """
    program = validate_program(program)
    end = len(program.events) if stop_after is None else stop_after
    if type(end) is not int or not 0 <= end <= len(program.events):
        raise ValueError('stop_after must name an event prefix')
    if type(max_branches) is not int or not 1 <= max_branches <= MAX_BRANCHES:
        raise ValueError('bounded positive branch allocation required')
    records = _records(program, end)
    if len(records) > max_branches:
        raise ValueError('complete branch budget exceeded')
    if checkpoint is not None:
        if initial is not None:
            raise ValueError('continuation cannot replace its retained input')
        history = validate_history(checkpoint)
        if history.program.digest != program.digest or end < history.cursor:
            raise ValueError('continuation cannot change or rewind the program')
    else:
        density = (program.initial_density if initial is None else
                   _density(initial, program.work_dimension, normalized=True))
        reference = density.shape[0]//program.work_dimension
        prep = program.preparation_resources if initial is None else ()
        history = History(program, 0, reference, _immutable(density),
                          'declared' if initial is None else 'supplied',
                          (Branch((), _immutable(density), prep),), prep, ())
    if len(records)*history.initial_density.size > MAX_DENSITY_ENTRIES:
        raise ValueError('complete density/reference allocation budget exceeded')
    for index in range(history.cursor, end):
        event = program.events[index]
        branches, total = [], history.expanded_tree_resources
        for branch in history.branches:
            if event.sequence:
                rho = branch.density
                for _ in range(event.repeat):
                    for instrument in event.sequence:
                        rho = _apply(instrument, rho, history.reference_dimension)[0]
                cost = _scaled(add_resources(*(i.resources for i in event.sequence)), event.repeat)
                outputs = (rho,)
            else:
                instrument = event.instruments[0 if event.control is None else branch.record[event.control]]
                outputs, cost = _apply(instrument, branch.density, history.reference_dimension), instrument.resources
            total = add_resources(total, cost)
            for outcome, rho in enumerate(outputs):
                record = (*branch.record, outcome) if event.visible else branch.record
                if not np.isfinite(rho).all():
                    raise ArithmeticError('executed channel produced nonfinite quantum output')
                branches.append(Branch(record, _immutable(rho), add_resources(branch.resources, cost)))
        mass = math.fsum(b.probability for b in branches)
        if not math.isfinite(mass) or abs(mass-1) > TOLERANCE:
            raise ArithmeticError('complete event lost probability mass')
        stage = Stage(index+1, event.name, len(branches), mass)
        history = replace(history, cursor=index+1, branches=tuple(branches), expanded_tree_resources=total,
                          stages=(*history.stages, stage), seal='')
    # The operation boundary checks mass without an eigensolve per microstep.
    # Check the complete output's positivity once before exposing it.
    for branch in history.branches:
        _density(branch.density, program.work_dimension, normalized=False)
    return replace(history, seal=_seal(history))


def branch_experiment(history, future_events, *, name, binding=None):
    """Calculate an alternative continuation of a retained complete prefix.

    This reuses simulator arithmetic, not a physically cloned quantum state.
    The returned program contains its original preparation and every executed
    prefix event, so exporting/importing its checkpoint independently replays
    the entire experiment. An earlier visible read stays in every branch;
    this operation never replaces a measured state by its unread precursor.
    """
    history = validate_history(history)
    if type(future_events) is not tuple:
        raise ValueError('an immutable future event sequence required')
    program = Program(name, history.work_dimension, history.program.initial_density,
                      (*history.program.events[:history.cursor], *future_events),
                      history.program.preparation_resources,
                      {'prefix_declaration': strict_json.loads(history.program.binding),
                       'alternative_continuation': strict_json.loads(_binding({} if binding is None else binding)),
                       'host_prefix_reuse_is_physical_copying': False})
    rebased = replace(history, program=program, seal='')
    rebased = replace(rebased, seal=_seal(rebased))
    return execute(program, checkpoint=rebased)


def complete_distance(left, right):
    left, right = validate_history(left), validate_history(right)
    if (left.records != right.records or left.work_dimension != right.work_dimension
            or left.reference_dimension != right.reference_dimension):
        raise ValueError('matching complete records and quantum output algebras required')
    return math.fsum(float(np.abs(np.linalg.eigvalsh(
        (a.density-b.density+(a.density-b.density).conj().T)/2)).sum())/2
        for a, b in zip(left.branches, right.branches, strict=True))


def _encode(value):
    return {'shape': list(value.shape), 'real': value.real.ravel().tolist(),
            'imag': value.imag.ravel().tolist()}


def _decode(value, size):
    if (type(value) is not dict or set(value) != {'shape', 'real', 'imag'}
            or type(value['shape']) is not list or any(type(x) is not int for x in value['shape'])
            or value['shape'] != [size, size]):
        raise ValueError('serialized quantum matrix has a different shape')
    for name in ('real', 'imag'):
        if (type(value[name]) is not list or len(value[name]) != size*size
                or any(type(x) not in (int, float) or not math.isfinite(x) for x in value[name])):
            raise ValueError('serialized quantum matrix requires finite numeric entries')
    # Assign real and imaginary lanes separately: adding real+1j*imag can
    # erase signed zeros and invalidate an otherwise exact local checkpoint.
    decoded = np.empty((size, size), dtype=complex)
    decoded.real = np.asarray(value['real']).reshape(size, size)
    decoded.imag = np.asarray(value['imag']).reshape(size, size)
    return _immutable(decoded)


def checkpoint_to_json(history):
    history = validate_history(history)
    value = {'schema': 'oph.primitive_observer.checkpoint.v1', **_header(history),
             'initial_density': _encode(history.initial_density),
             'branches': [{'record': b.record, 'resources': b.resources, 'density': _encode(b.density)}
                          for b in history.branches], 'seal': history.seal}
    text = _json(value)
    if len(text.encode()) > MAX_CHECKPOINT_BYTES:
        raise ValueError('checkpoint serialization budget exceeded')
    return text


def checkpoint_from_json(program, payload):
    """Import and independently replay; a resealed fabricated state fails."""
    program = validate_program(program)
    if type(payload) is not str or len(payload.encode()) > MAX_CHECKPOINT_BYTES:
        raise ValueError('bounded checkpoint JSON text required')
    value = strict_json.loads(payload)
    keys = {'schema', 'program_digest', 'cursor', 'reference_dimension', 'preparation_mode',
            'expanded_tree_resources', 'stages', 'branches', 'initial_density', 'seal'}
    if (type(value) is not dict or set(value) != keys
            or value['schema'] != 'oph.primitive_observer.checkpoint.v1'
            or value['program_digest'] != program.digest):
        raise ValueError('checkpoint belongs to a different declared program')
    cursor, reference = value['cursor'], value['reference_dimension']
    if (type(cursor) is not int or not 0 <= cursor <= len(program.events)
            or type(reference) is not int or reference < 1):
        raise ValueError('invalid checkpoint cursor or reference dimension')
    size = program.work_dimension*reference
    records = _records(program, cursor)
    if len(records)*size**2 > MAX_DENSITY_ENTRIES:
        raise ValueError('checkpoint density allocation budget exceeded')
    if (type(value['branches']) is not list or len(value['branches']) != len(records)
            or type(value['stages']) is not list or len(value['stages']) != cursor):
        raise ValueError('complete serialized branches and stage journal required')
    def cost(rows):
        if type(rows) is not list or any(type(row) is not list for row in rows):
            raise ValueError('serialized resource rows must be arrays')
        return _cost(tuple(tuple(row) for row in rows))
    branches, stages = [], []
    for row in value['branches']:
        if type(row) is not dict or set(row) != {'record', 'resources', 'density'} or type(row['record']) is not list:
            raise ValueError('complete serialized branch required')
        branches.append(Branch(tuple(row['record']), _decode(row['density'], size), cost(row['resources'])))
    for row in value['stages']:
        if type(row) is not dict or set(row) != {'cursor', 'name', 'branches', 'probability_mass'}:
            raise ValueError('complete serialized stage required')
        stages.append(Stage(**row))
    history = History(program, cursor, reference, _decode(value['initial_density'], size),
                      value['preparation_mode'], tuple(branches), cost(value['expanded_tree_resources']),
                      tuple(stages), value['seal'])
    history = validate_history(history)
    replay = execute(program, stop_after=cursor,
                     initial=None if history.preparation_mode == 'declared' else history.initial_density)
    if (complete_distance(history, replay) > TOLERANCE
            or history.expanded_tree_resources != replay.expanded_tree_resources
            or any(a.resources != b.resources for a, b in zip(history.branches, replay.branches, strict=True))
            or any(abs(a.probability_mass-b.probability_mass) > TOLERANCE
                   for a, b in zip(history.stages, replay.stages, strict=True))):
        raise ValueError('checkpoint does not replay its declared complete prefix')
    return history
