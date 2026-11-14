"""Complete mixed-state execution of the canonical observer program.

The visible operation grammar, source maps and preparation are owned by
protocol.py.  An Environment attaches unobserved, independently refreshed
source channels at declared boundaries.  Its draws are summed as CP maps,
never as amplitudes or selected outcomes.  Work and an untouched reference
remain in each unnormalized visible-record block through continuation.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
from itertools import product
import json
import math

import numpy as np

from tools import strict_json
from . import protocol
from . import open_environment as environment


TOLERANCE = 2e-9
MAX_DENSITY_ENTRIES = 2_000_000
RESOURCE_KEYS = (*protocol.RESOURCE_KEYS, 'independent_environment_draws')
NOISE_BACKENDS = ('source', 'gaussian')


def _immutable(value):
    return protocol._immutable(value)


def _cost(value=None, **extra):
    value = {} if value is None else dict(value)
    return tuple((key, value.get(key, 0)+extra.get(key, 0)) for key in RESOURCE_KEYS)


def _add(*values):
    return tuple((key, sum(dict(value)[key] for value in values)) for key in RESOURCE_KEYS)


def _density(value, size, *, normalized):
    protocol._original_scalars(value)
    if np.ma.isMaskedArray(value):
        raise ValueError('masked density operators are inadmissible')
    raw = np.asarray(value)
    if (raw.ndim != 2 or raw.shape[0] != raw.shape[1] or not raw.shape[0]
            or raw.shape[0] % size or raw.dtype.kind not in 'iufc'
            or raw.size > MAX_DENSITY_ENTRIES or not np.isfinite(raw).all()):
        raise ValueError('finite square density on work times reference required')
    rho = np.array(raw, dtype=complex, copy=True)
    if not np.isfinite(rho).all():
        raise ValueError('density must remain finite in the execution representation')
    scale = float(np.max(np.abs(rho)))
    if np.linalg.norm(rho-rho.conj().T, 2) > TOLERANCE*scale:
        raise ValueError('density operator must be Hermitian')
    eigenvalues = np.linalg.eigvalsh((rho+rho.conj().T)/2)
    mass = float(np.trace(rho).real)
    if eigenvalues[0] < -TOLERANCE*max(scale, mass):
        raise ValueError('density operator must be positive semidefinite')
    if normalized and abs(mass-1) > TOLERANCE:
        raise ValueError('initial density operator must have unit trace')
    return _immutable(rho)


def conjugate(matrix, rho, reference_dimension):
    """Apply a work map to all reference blocks; never trace the reference."""
    size = matrix.shape[0]
    tensor = rho.reshape(size, reference_dimension, size, reference_dimension)
    result = np.einsum('ai,irjs,bj->arbs', matrix, tensor, matrix.conj(), optimize=True)
    return result.reshape(rho.shape)


@dataclass(frozen=True)
class MixedBranch:
    record: tuple[int, ...]
    density: np.ndarray
    resources: tuple[tuple[str, int], ...]

    @property
    def probability(self):
        return float(np.trace(self.density).real)


@dataclass(frozen=True)
class MixedStage:
    cursor: int
    kind: str
    branch_count: int
    probability_mass: float
    environment_applications_per_path: int


@dataclass(frozen=True)
class MixedHistory:
    program: protocol.Program
    environment: environment.Environment
    backend: str
    noise_backend: str
    cursor: int
    reference_dimension: int
    initial_density: np.ndarray
    preparation_mode: str
    preparation_resources: tuple[tuple[str, int], ...]
    branches: tuple[MixedBranch, ...]
    expanded_tree_resources: tuple[tuple[str, int], ...]
    stages: tuple[MixedStage, ...]
    seal: str = ''

    @property
    def probability_mass(self):
        return math.fsum(branch.probability for branch in self.branches)

    @property
    def records(self):
        return tuple(branch.record for branch in self.branches)


def _specification(history):
    return {
        'program': history.program.specification(),
        'environment': history.environment.specification(),
        'backend': history.backend, 'noise_backend': history.noise_backend,
        'cursor': history.cursor, 'reference_dimension': history.reference_dimension,
        'preparation_mode': history.preparation_mode,
        'preparation_resources': history.preparation_resources,
        'expanded_tree_resources': history.expanded_tree_resources,
        'stages': [stage.__dict__ for stage in history.stages],
        'branches': [{'record': branch.record, 'resources': branch.resources}
                     for branch in history.branches],
    }


def _seal(history):
    digest = hashlib.sha256(protocol._canonical(_specification(history)))
    for array in (history.initial_density, *(branch.density for branch in history.branches)):
        digest.update(protocol._canonical(list(array.shape)))
        digest.update(np.ascontiguousarray(array, dtype=complex).tobytes())
    return digest.hexdigest()


def _check_cost(value):
    if (type(value) is not tuple or len(value) != len(RESOURCE_KEYS)
            or any(type(row) is not tuple or len(row) != 2 for row in value)
            or tuple(row[0] for row in value) != RESOURCE_KEYS
            or any(type(row[1]) is not int or row[1] < 0 for row in value)):
        raise ValueError('complete nonnegative integer resource ledger required')


def _check(history, *, seal=True):
    if type(history) is not MixedHistory:
        raise ValueError('a complete mixed history required')
    program = protocol.validate_program(history.program)
    env = environment.validated(program, history.environment)
    if (type(history.backend) is not str or history.backend not in protocol.BACKENDS
            or type(history.noise_backend) is not str or history.noise_backend not in NOISE_BACKENDS):
        raise ValueError('declared operation and noise backends required')
    if (type(history.cursor) is not int or not 0 <= history.cursor <= len(program.operations)
            or type(history.reference_dimension) is not int or history.reference_dimension < 1
            or history.preparation_mode not in ('declared_vacuum', 'supplied_density')):
        raise ValueError('invalid mixed-history domain')
    size = program.work_dimension*history.reference_dimension
    reads = sum(op.kind in protocol.READS for op in program.operations[:history.cursor])
    if (type(history.branches) is not tuple or len(history.branches) != 2**reads
            or len(history.branches)*size**2 > MAX_DENSITY_ENTRIES):
        raise ValueError('mixed history exceeds its complete record or storage budget')
    if not isinstance(history.initial_density, np.ndarray) or history.initial_density.shape != (size, size):
        raise ValueError('initial reference dimension differs from its contract')
    _density(history.initial_density, program.work_dimension, normalized=True)
    _check_cost(history.preparation_resources)
    _check_cost(history.expanded_tree_resources)
    for branch, expected in zip(history.branches, product((0, 1), repeat=reads), strict=True):
        if (type(branch) is not MixedBranch or type(branch.record) is not tuple
                or any(type(bit) is not int for bit in branch.record)
                or branch.record != expected or not isinstance(branch.density, np.ndarray)
                or branch.density.shape != (size, size)):
            raise ValueError('every ordered visible record and full quantum block required')
        _density(branch.density, program.work_dimension, normalized=False)
        _check_cost(branch.resources)
    if abs(history.probability_mass-1) > TOLERANCE:
        raise ArithmeticError('complete mixed history lost probability mass')
    if type(history.stages) is not tuple or len(history.stages) != history.cursor:
        raise ValueError('mixed-history stage journal is incomplete')
    for cursor, stage in enumerate(history.stages, 1):
        count = sum(op.kind in protocol.READS for op in program.operations[:cursor])
        events = sum(event.after <= cursor for event in env.events)
        if (type(stage) is not MixedStage or type(stage.cursor) is not int or stage.cursor != cursor
                or stage.kind != program.operations[cursor-1].kind
                or type(stage.branch_count) is not int or stage.branch_count != 2**count
                or type(stage.environment_applications_per_path) is not int
                or stage.environment_applications_per_path != events
                or type(stage.probability_mass) not in (int, float)
                or not math.isfinite(stage.probability_mass)
                or abs(stage.probability_mass-1) > TOLERANCE):
            raise ValueError('stage journal differs from the declared complete history')
    if seal and (not history.seal or history.seal != _seal(history)):
        raise ValueError('mixed checkpoint integrity digest differs from its content')
    # Frozen dataclasses still have mutable backing dictionaries. Detach every
    # envelope when accepting a caller-owned checkpoint, including no-op resume.
    return replace(history, program=program, environment=env,
                   initial_density=_immutable(history.initial_density),
                   branches=tuple(replace(branch, density=_immutable(branch.density))
                                  for branch in history.branches),
                   stages=tuple(replace(stage) for stage in history.stages))


def _after(history, cursor):
    branches, total = history.branches, history.expanded_tree_resources
    for event in history.environment.events:
        if event.after != cursor:
            continue
        metadata = environment.event_metadata(history.program, event)
        resources = metadata['resources_per_realization'] if history.noise_backend == 'source' else {}
        cost = _cost(resources, independent_environment_draws=1)
        output = []
        for branch in branches:
            rho = environment.apply(branch.density, history.program, event, backend=history.noise_backend)
            output.append(MixedBranch(branch.record, _immutable(rho), _add(branch.resources, cost)))
            total = _add(total, cost)
        branches = tuple(output)
    return replace(history, branches=branches, expanded_tree_resources=total, seal='')


def execute(program, env=None, *, backend='source', noise_backend='source', initial=None,
            checkpoint=None, stop_after=None, max_branches=256):
    """Execute or resume the same program on complete density/reference blocks.

    ``initial`` is a normalized work-times-reference density matrix. It replaces
    vacuum preparation; its reference is retained. Source and Gaussian noise
    backends are distinct and included in the checkpoint identity. No noiseless
    prepared-state certificate is inherited by an arbitrary supplied density.
    """
    program = protocol.validate_program(program)
    env = environment.validated(program, environment.Environment(()) if env is None else env)
    if (type(backend) is not str or backend not in protocol.BACKENDS
            or type(noise_backend) is not str or noise_backend not in NOISE_BACKENDS):
        raise ValueError('declared operation and environment backends required')
    if type(max_branches) is not int or not 1 <= max_branches <= 2**protocol.MAX_READS:
        raise ValueError('integer visible branch allocation from 1 through 256 required')
    end = len(program.operations) if stop_after is None else stop_after
    if type(end) is not int or not 0 <= end <= len(program.operations):
        raise ValueError('stop_after must name an operation prefix')
    reads = sum(op.kind in protocol.READS for op in program.operations[:end])
    if 2**reads > max_branches:
        raise ValueError('visible record allocation budget exceeded')
    if checkpoint is not None:
        if initial is not None:
            raise ValueError('continuation cannot replace the retained input')
        history = _check(checkpoint)
        if (history.program != program or history.environment != env
                or history.backend != backend or history.noise_backend != noise_backend):
            raise ValueError('checkpoint belongs to a different program, environment or backend')
        if end < history.cursor:
            raise ValueError('continuation cannot rewind the history')
    else:
        if initial is None:
            if (2**reads)*program.work_dimension**2 > MAX_DENSITY_ENTRIES:
                raise ValueError('complete density/reference allocation budget exceeded')
            prepared = protocol.execute(program, backend=backend, stop_after=0)
            before = prepared.initial_state.ravel()
            after = prepared.branches[0].state.ravel()
            density, state = np.outer(before, before.conj()), np.outer(after, after.conj())
            prep = _cost(prepared.preparation_resources)
            resources = _cost(prepared.branches[0].resources)
            mode, reference = 'declared_vacuum', 1
        else:
            density = _density(initial, program.work_dimension, normalized=True)
            reference = density.shape[0]//program.work_dimension
            if (2**reads)*density.size > MAX_DENSITY_ENTRIES:
                raise ValueError('complete density/reference allocation budget exceeded')
            kick, kick_cost = protocol._kick(program.dimension, program.native_current, backend)
            state = conjugate(kick, density, reference)
            prep, resources, mode = _cost(), _cost(kick_cost), 'supplied_density'
        history = MixedHistory(program, env, backend, noise_backend, 0, reference,
                               _immutable(density), mode, prep,
                               (MixedBranch((), _immutable(state), resources),), resources, ())
        history = _after(history, 0)
    if (2**reads)*history.initial_density.size > MAX_DENSITY_ENTRIES:
        raise ValueError('complete density/reference allocation budget exceeded')
    for cursor in range(history.cursor, end):
        operation = program.operations[cursor]
        branches, total = [], history.expanded_tree_resources
        for branch in history.branches:
            selected = branch.record[operation.record] if operation.kind == 'feedback' else None
            maps, raw_cost = protocol._operation_maps(program, operation, selected, backend)
            cost = _cost(raw_cost)
            total = _add(total, cost)
            for bit, matrix in enumerate(maps):
                record = (*branch.record, bit) if operation.kind in protocol.READS else branch.record
                state = conjugate(matrix, branch.density, history.reference_dimension)
                branches.append(MixedBranch(record, _immutable(state), _add(branch.resources, cost)))
        history = replace(history, cursor=cursor+1, branches=tuple(branches),
                          expanded_tree_resources=total, seal='')
        history = _after(history, cursor+1)
        stage = MixedStage(cursor+1, operation.kind, len(branches), history.probability_mass,
                           sum(event.after <= cursor+1 for event in env.events))
        history = replace(history, stages=(*history.stages, stage))
        _check(history, seal=False)
    _check(history, seal=False)
    return replace(history, seal=_seal(history))


def prefix_probabilities(history, length):
    """Marginalize later outcomes without conditioning on any of them."""
    _check(history)
    count = len(history.branches[0].record)
    if type(length) is not int or not 0 <= length <= count:
        raise ValueError('an existing visible prefix length required')
    return {prefix: math.fsum(b.probability for b in history.branches if b.record[:length] == prefix)
            for prefix in product((0, 1), repeat=length)}


def complete_distance(left, right):
    """Numerical complete CQ trace distance; includes the retained reference."""
    for history in (left, right):
        _check(history)
    if (left.program != right.program or left.cursor != right.cursor
            or left.records != right.records or left.reference_dimension != right.reference_dimension):
        raise ValueError('matching programs, record menus and reference registers required')
    return math.fsum(float(np.abs(np.linalg.eigvalsh(
        (a.density-b.density + (a.density-b.density).conj().T)/2)).sum())/2
        for a, b in zip(left.branches, right.branches, strict=True))


def checkpoint_to_json(history):
    """Export simulator data. This does not physically copy quantum memory."""
    _check(history)
    value = {'schema': 'oph.common_source.mixed_checkpoint.v1', **_specification(history),
             'initial_density': protocol._encoded_array(history.initial_density),
             'branches': [{'record': b.record, 'resources': b.resources,
                           'density': protocol._encoded_array(b.density)} for b in history.branches],
             'seal': history.seal}
    result = protocol._canonical(value)
    if len(result) > protocol.MAX_CHECKPOINT_BYTES:
        raise ValueError('mixed checkpoint exceeds serialization byte budget')
    return result.decode()


def _decode_cost(value):
    if type(value) is not list or any(type(row) is not list for row in value):
        raise ValueError('serialized resource rows must be arrays')
    result = tuple(tuple(row) for row in value)
    _check_cost(result)
    return result


def checkpoint_from_json(program, env, payload, *, backend='source', noise_backend='source'):
    """Strictly import and independently replay the declared complete prefix.

    A digest alone is never accepted as evidence that its density belongs to
    the program. Replay also verifies resources and chronology after resealing.
    """
    program = protocol.validate_program(program)
    env = environment.validated(program, env)
    if type(payload) is not str or len(payload.encode()) > protocol.MAX_CHECKPOINT_BYTES:
        raise ValueError('bounded JSON checkpoint text required')
    value = strict_json.loads(payload)
    expected = {'schema', 'program', 'environment', 'backend', 'noise_backend', 'cursor',
                'reference_dimension', 'preparation_mode', 'preparation_resources',
                'expanded_tree_resources', 'stages', 'branches', 'initial_density', 'seal'}
    protocol._keys(value, expected)
    if (value['schema'] != 'oph.common_source.mixed_checkpoint.v1'
            or protocol._canonical(value['program']) != protocol._canonical(program.specification())
            or protocol._canonical(value['environment']) != protocol._canonical(env.specification())
            or value['backend'] != backend or value['noise_backend'] != noise_backend):
        raise ValueError('serialized mixed checkpoint has a different declared identity')
    reference, cursor = value['reference_dimension'], value['cursor']
    if (type(reference) is not int or reference < 1 or type(cursor) is not int
            or not 0 <= cursor <= len(program.operations)):
        raise ValueError('invalid serialized reference dimension or cursor')
    size = reference*program.work_dimension
    reads = sum(op.kind in protocol.READS for op in program.operations[:cursor])
    if (2**reads)*size**2 > MAX_DENSITY_ENTRIES:
        raise ValueError('serialized density allocation exceeds the declared budget')
    if type(value['branches']) is not list or len(value['branches']) != 2**reads:
        raise ValueError('complete serialized visible records required')
    if type(value['stages']) is not list or len(value['stages']) != cursor:
        raise ValueError('complete serialized stage journal required')
    branches, stages = [], []
    for row in value['branches']:
        protocol._keys(row, ('record', 'density', 'resources'))
        if type(row['record']) is not list:
            raise ValueError('serialized record must be an array')
        branches.append(MixedBranch(tuple(row['record']), protocol._decode_array(row['density'], size, size),
                                    _decode_cost(row['resources'])))
    for row in value['stages']:
        protocol._keys(row, ('cursor', 'kind', 'branch_count', 'probability_mass',
                             'environment_applications_per_path'))
        stages.append(MixedStage(**row))
    history = MixedHistory(program, env, backend, noise_backend, cursor, reference,
                           protocol._decode_array(value['initial_density'], size, size),
                           value['preparation_mode'], _decode_cost(value['preparation_resources']),
                           tuple(branches), _decode_cost(value['expanded_tree_resources']),
                           tuple(stages), value['seal'])
    _check(history)
    replay = execute(program, env, backend=backend, noise_backend=noise_backend, stop_after=cursor,
                     initial=None if history.preparation_mode == 'declared_vacuum' else history.initial_density)
    if (complete_distance(history, replay) > TOLERANCE
            or np.linalg.norm(history.initial_density-replay.initial_density, 2) > TOLERANCE
            or history.preparation_resources != replay.preparation_resources
            or history.expanded_tree_resources != replay.expanded_tree_resources
            or any(a.resources != b.resources for a, b in zip(history.branches, replay.branches, strict=True))
            or any(abs(a.probability_mass-b.probability_mass) > TOLERANCE
                   for a, b in zip(history.stages, replay.stages, strict=True))):
        raise ValueError('mixed checkpoint does not reproduce its declared complete prefix')
    return history
