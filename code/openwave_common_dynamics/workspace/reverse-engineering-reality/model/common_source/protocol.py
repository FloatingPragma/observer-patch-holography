"""Finite recorded programs on the existing interacting source register.

The canonical observer keeps one prepared work state through every read and
control. A branch stores its unnormalized work/reference vector and its full
prefix Kraus operator. Checkpoints are simulator data, not an implementation
of physical quantum-state storage, transport, or cloning.

The source backend contracts the existing emitted source words. The product
and continuous backends use independent target matrices for comparison.
Exact rational specifications describe the ideal maps. Numerical backends
evaluate their angles as binary64; diagnostics are not roundoff certificates.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import json
import math

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar import compiler, field
from . import experiment as e, scalar_coupling as scalar, smooth_protocol as smooth


KINDS = frozenset(('evolve', 'gauge_read', 'smooth_read', 'sharp_read', 'feedback'))
READS = frozenset(('gauge_read', 'smooth_read', 'sharp_read'))
BACKENDS = frozenset(('source', 'product', 'continuous'))
TOLERANCE = 2e-9
MAX_OPERATIONS = 64
MAX_READS = 8
MAX_CHECKPOINT_BYTES = 16_000_000
RESOURCE_KEYS = (*e.COUNT_KEYS, 'fresh_pointer_blanks', 'committed_classical_bits',
                 'work_preparations')


def _rational(value, name):
    if type(value) not in (int, F):
        raise ValueError(f'{name} must be an exact integer or Fraction')
    return F(int(value.numerator), int(value.denominator)) if type(value) is F else F(value)


def _finite_float(value):
    try:
        return math.isfinite(float(value))
    except (OverflowError, ValueError):
        return False


@dataclass(frozen=True)
class Operation:
    """One declared operation; feedback reads an earlier zero-based record."""

    kind: str
    time: F = F(0)
    record: int | None = None
    enabled: bool = True

    def __post_init__(self):
        if type(self.kind) is not str or self.kind not in KINDS:
            raise ValueError('an admitted observer operation is required')
        time = _rational(self.time, 'model time')
        if time < 0 or (self.kind != 'evolve' and time != 0):
            raise ValueError('only evolution admits a nonnegative model duration')
        if type(self.enabled) is not bool or (self.kind != 'feedback' and not self.enabled):
            raise ValueError('only feedback admits an explicit enable switch')
        if self.kind == 'feedback':
            if type(self.record) is not int or self.record < 0:
                raise ValueError('feedback requires an earlier binary record index')
        elif self.record is not None:
            raise ValueError('only feedback selects a prior record')
        object.__setattr__(self, 'time', time)


@dataclass(frozen=True)
class Program:
    """One preparation, fixed action and a finite causal sequence of operations.

    Durations, readout meanings, and feedback addresses are chosen before the
    run. The operation count is finite; execution separately enforces its
    explicit branch-allocation budget. No operation resets the working state.
    """

    dimension: int = 2
    steps: int = 64
    coupling: F = F(1, 4)
    native_current: int = 0
    operations: tuple[Operation, ...] = ()

    def __post_init__(self):
        if type(self.dimension) is not int or self.dimension not in (2, 4):
            raise ValueError('executed onsite dimension 2 or 4 required')
        if type(self.steps) is not int or self.steps < 1:
            raise ValueError('a positive integer split count required')
        coupling = _rational(self.coupling, 'coupling')
        if not _finite_float(coupling):
            raise ValueError('coupling must have a finite execution representation')
        if type(self.native_current) is not int or self.native_current not in (-1, 0, 1):
            raise ValueError('native current -1, 0 or 1 required')
        if type(self.operations) is not tuple:
            raise ValueError('the finite operation sequence must be an immutable tuple')
        if len(self.operations) > MAX_OPERATIONS:
            raise ValueError('program exceeds the 64-operation execution budget')
        reads = 0
        for operation in self.operations:
            if type(operation) is not Operation:
                raise ValueError('every program entry must be an Operation')
            if not _finite_float(operation.time):
                raise ValueError('time must have a finite execution representation')
            if operation.kind == 'feedback' and operation.record >= reads:
                raise ValueError('feedback cannot read a future or absent record')
            reads += operation.kind in READS
        if reads > MAX_READS:
            raise ValueError('program exceeds the eight-read execution budget')
        object.__setattr__(self, 'coupling', coupling)

    @property
    def work_dimension(self):
        return 4*self.dimension**2

    @property
    def digest(self):
        return hashlib.sha256(_canonical(self.specification())).hexdigest()

    def specification(self):
        return {'schema': 'oph.common_source.program.v1',
                'dimension': self.dimension, 'steps': self.steps,
                'coupling': str(self.coupling), 'native_current': self.native_current,
                'preparation': 'source-prepared projected free vacuum and blank gauge bits',
                'operations': [{'kind': op.kind, 'time': str(op.time),
                                'record': op.record, 'enabled': op.enabled}
                               for op in self.operations],
                'bit_meanings': {'gauge_read': 'species-0 right occupation',
                                 'smooth_read': '0: (I+sin Q0)/2; 1: (I-sin Q0)/2',
                                 'sharp_read': '1: Q0>0; 0: Q0<0'},
                'feedback': 'selected committed bit 1 applies the existing gauge X pulse',
                'clock': 'model durations; fixed Gibbs reference gives t=-beta*s on unitary segments'}


def validate_program(program):
    """Revalidate and detach the complete canonical program before execution.

    Frozen dataclasses still expose mutable backing dictionaries. Rebuilding
    every operation checks its own contract before Program checks scheduling,
    and reconstructs Fraction fields instead of retaining caller-owned values.
    The tuple and operation budgets are checked before traversing its entries.
    """
    if type(program) is not Program:
        raise ValueError('a typed canonical Program required')
    operations = program.operations
    if type(operations) is not tuple or len(operations) > MAX_OPERATIONS:
        raise ValueError('a bounded immutable tuple of program operations required')
    if any(type(operation) is not Operation for operation in operations):
        raise ValueError('every program entry must be a canonical Operation')
    return replace(program, operations=tuple(replace(operation) for operation in operations))


def smooth_program(dimension=2, steps=64, *, native_current=0, feedback=True, coupling=F(1, 4)):
    """The existing smooth history as an explicit composable program."""
    return Program(dimension, steps, coupling, native_current, (
        Operation('evolve', F(1)), Operation('gauge_read'),
        Operation('feedback', record=0, enabled=feedback),
        Operation('evolve', F(1, 2)), Operation('smooth_read')))


def sharp_program(dimension=2, steps=64, *, native_current=0, feedback=True, coupling=F(1, 4)):
    """Compatibility replay only; sharp readout has different refinement scope."""
    base = smooth_program(dimension, steps, native_current=native_current,
                          feedback=feedback, coupling=coupling)
    return replace(base, operations=(*base.operations[:-1], Operation('sharp_read')))


def continuation_program(dimension=2, steps=64, *, native_current=0, coupling=F(1, 4), clock_time=F(1, 2)):
    """Two smooth reads and a later gauge clock probe on one retained state."""
    base = smooth_program(dimension, steps, native_current=native_current, coupling=coupling)
    return replace(base, operations=(*base.operations,
        Operation('evolve', F(1, 2)), Operation('smooth_read'),
        Operation('feedback', record=2), Operation('evolve', clock_time),
        Operation('gauge_read')))


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def _immutable(array):
    array = np.ascontiguousarray(array, dtype=np.complex128)
    # An immutable bytes owner prevents callers from re-enabling NumPy writes.
    return np.frombuffer(array.tobytes(), dtype=np.complex128).reshape(array.shape)


def _cost(metadata=None, **extra):
    source = (metadata or {}).get('resources', {})
    return tuple((key, int(source.get(key, 0))+extra.get(key, 0)) for key in RESOURCE_KEYS)


def _sum_cost(*costs):
    return tuple((key, sum(dict(cost).get(key, 0) for cost in costs)) for key in RESOURCE_KEYS)


@dataclass(frozen=True)
class Branch:
    record: tuple[int, ...]
    state: np.ndarray
    kraus: np.ndarray
    resources: tuple[tuple[str, int], ...]

    @property
    def probability(self):
        return float(np.vdot(self.state, self.state).real)


@dataclass(frozen=True)
class Stage:
    cursor: int
    kind: str
    record_count: int
    branch_count: int
    completeness_error: float
    probability_mass: float


@dataclass(frozen=True)
class History:
    program: Program
    backend: str
    cursor: int
    branches: tuple[Branch, ...]
    initial_state: np.ndarray
    preparation_resources: tuple[tuple[str, int], ...]
    expanded_tree_resources: tuple[tuple[str, int], ...]
    stages: tuple[Stage, ...]
    preparation_mode: str
    seal: str = ''

    @property
    def probability_mass(self):
        return math.fsum(branch.probability for branch in self.branches)

    @property
    def completeness_error(self):
        total = sum((b.kraus.conj().T@b.kraus for b in self.branches),
                    np.zeros((self.program.work_dimension,)*2, complex))
        return float(np.linalg.norm(total-np.eye(self.program.work_dimension), 2))

    @property
    def records(self):
        return tuple(branch.record for branch in self.branches)


def _seal(history):
    digest = hashlib.sha256(_canonical({
        'program': history.program.specification(), 'backend': history.backend,
        'cursor': history.cursor, 'preparation_mode': history.preparation_mode,
        'preparation_resources': history.preparation_resources,
        'expanded_tree_resources': history.expanded_tree_resources,
        'stages': [stage.__dict__ for stage in history.stages]}))
    for array in (history.initial_state,):
        digest.update(_canonical(list(array.shape)))
        digest.update(np.ascontiguousarray(array, dtype=np.complex128).tobytes())
    for branch in history.branches:
        digest.update(_canonical([branch.record, branch.resources]))
        for array in (branch.state, branch.kraus):
            digest.update(_canonical(list(array.shape)))
            digest.update(np.ascontiguousarray(array, dtype=np.complex128).tobytes())
    return digest.hexdigest()


def _checked(history):
    if (not math.isfinite(history.completeness_error)
            or history.completeness_error > TOLERANCE):
        raise ArithmeticError('complete prefix Kraus maps are not trace preserving')
    if not math.isfinite(history.probability_mass) or abs(history.probability_mass-1) > TOLERANCE:
        raise ArithmeticError('complete unnormalized history lost probability mass')
    return replace(history, seal=_seal(history))


def validate_checkpoint(program, checkpoint, backend='source'):
    """Check simulator checkpoint identity, integrity, and complete prefix.

    This detects accidental or malicious edits to the saved envelope. It is
    not authentication against an actor who also recomputes its public digest.
    """
    if type(checkpoint) is not History:
        raise ValueError('a typed program and complete History checkpoint required')
    program = validate_program(program)
    saved_program = validate_program(checkpoint.program)
    if backend not in BACKENDS or checkpoint.backend != backend or saved_program != program:
        raise ValueError('checkpoint belongs to a different program or backend')
    if type(checkpoint.cursor) is not int or not 0 <= checkpoint.cursor <= len(program.operations):
        raise ValueError('checkpoint cursor is not an admitted program prefix')
    reads = sum(op.kind in READS for op in program.operations[:checkpoint.cursor])
    if checkpoint.records != tuple(product((0, 1), repeat=reads)):
        raise ValueError('checkpoint does not retain every ordered binary prefix')
    shape = checkpoint.initial_state.shape
    if len(shape) != 2 or shape[0] != program.work_dimension or shape[1] < 1:
        raise ValueError('checkpoint initial work/reference shape is invalid')
    if not np.isfinite(checkpoint.initial_state).all():
        raise ValueError('checkpoint initial state contains nonfinite quantum data')
    if len(checkpoint.stages) != checkpoint.cursor:
        raise ValueError('checkpoint stage journal does not match its prefix')
    for branch in checkpoint.branches:
        if branch.state.shape != shape or branch.kraus.shape != (program.work_dimension,)*2:
            raise ValueError('checkpoint work/reference or Kraus shape is invalid')
        if not np.isfinite(branch.state).all() or not np.isfinite(branch.kraus).all():
            raise ValueError('checkpoint contains nonfinite quantum data')
        if np.linalg.norm(branch.state-branch.kraus@checkpoint.initial_state) > TOLERANCE:
            raise ValueError('checkpoint quantum state is detached from its prefix map')
    if not checkpoint.seal or _seal(checkpoint) != checkpoint.seal:
        raise ValueError('checkpoint integrity digest does not match its contents')
    if (not math.isfinite(checkpoint.completeness_error)
            or not math.isfinite(checkpoint.probability_mass)
            or checkpoint.completeness_error > TOLERANCE
            or abs(checkpoint.probability_mass-1) > TOLERANCE):
        raise ValueError('checkpoint is not a complete normalized recorded instrument')
    # The resumed/exported envelope must not share its program with the
    # caller's checkpoint, even when no further operation will execute.
    return replace(checkpoint, program=program)


def _original_scalars(value):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or value is None:
        raise ValueError('masked, missing, or Boolean quantum amplitudes are inadmissible')
    if isinstance(value, (tuple, list)):
        for element in value:
            _original_scalars(element)
    elif isinstance(value, np.ndarray) and value.dtype.kind not in 'iufc':
        raise ValueError('quantum input arrays require an original non-Boolean numeric dtype')


def _input_state(initial, size):
    _original_scalars(initial)
    if np.ma.isMaskedArray(initial):
        raise ValueError('masked input is not a work/reference state')
    value = np.asarray(initial)
    if (value.dtype.kind not in 'iufc' or value.ndim not in (1, 2)
            or value.shape[0] != size or not np.isfinite(value).all()):
        raise ValueError('finite numeric work vector or work/reference columns required')
    value = np.asarray(value, dtype=complex).reshape(size, -1)
    if value.shape[1] < 1 or abs(float(np.vdot(value, value).real)-1) > TOLERANCE:
        raise ValueError('the complete input state must have unit norm without renormalization')
    return _immutable(value)


def _gauge_x(size):
    result = np.zeros((size, size))
    for index in range(size):
        result[index ^ 1, index] = 1
    return result


@lru_cache(maxsize=48)
def _kick(dimension, current, backend):
    size = 4*dimension**2
    if backend == 'source':
        word = [compiler.pulse('X', -current*e.FEEDBACK_ANGLE, [0])]
        matrix, _ = compiler.stabilize_executed_unitary(compiler.execute(word, scalar.qubits(dimension)))
        cost = _cost({'resources': e._word_resources(word)})
    else:
        matrix = expm(-1j*current*e.FEEDBACK_ANGLE*_gauge_x(size))
        cost = _cost()
    return _immutable(matrix), cost


@lru_cache(maxsize=256)
def _operation_maps(program, operation, selected_bit, backend):
    """Complete local maps; independent backends never call source execution."""
    d, size = program.dimension, program.work_dimension
    kind = operation.kind
    metadata = None
    if kind == 'evolve':
        time = float(operation.time)
        if backend == 'source':
            matrix, metadata = scalar.source_unitary(d, program.steps, time, float(program.coupling))
        elif backend == 'product':
            matrix = scalar.independent_split_unitary(d, program.steps, time, float(program.coupling))
        else:
            matrix = expm(-1j*time*scalar.ideal_hamiltonian(d, float(program.coupling)))
        maps = (matrix,)
    elif kind == 'feedback':
        if selected_bit not in (0, 1):
            raise ValueError('a committed binary feedback symbol required')
        if backend == 'source':
            matrix, metadata = e.feedback(d, selected_bit, operation.enabled)
        else:
            matrix = expm(-1j*(e.FEEDBACK_ANGLE if operation.enabled and selected_bit else 0.)*_gauge_x(size))
        maps = (matrix,)
    elif kind == 'gauge_read':
        if backend == 'source':
            maps, metadata = e.pointer(d, 0)
        else:
            maps = tuple(e.projector(size, 0, bit) for bit in (0, 1))
    elif kind == 'smooth_read':
        if backend == 'source':
            maps, metadata = smooth.pointer(d)
        else:
            maps = smooth.ideal_kraus(d)
    elif kind == 'sharp_read':
        if backend == 'source':
            maps, metadata = e.scalar_pointer(d)
        else:
            maps = e.scalar_projectors(d)
    else:
        raise ValueError('unrecognized operation')
    if len(maps) != (2 if kind in READS else 1):
        raise ArithmeticError('operation omitted or added an outcome in its declared menu')
    for matrix in maps:
        raw = np.asarray(matrix)
        if (raw.shape != (size, size) or raw.dtype.kind not in 'iufc'
                or not np.isfinite(raw).all()):
            raise ArithmeticError('local Kraus map must be finite on the declared working algebra')
    maps = tuple(_immutable(matrix) for matrix in maps)
    complete = sum((matrix.conj().T@matrix for matrix in maps), np.zeros((size, size), complex))
    if float(np.linalg.norm(complete-np.eye(size), 2)) > TOLERANCE:
        raise ArithmeticError('local operation has an incomplete Kraus family')
    extra = {'fresh_pointer_blanks': int(kind in READS),
             'committed_classical_bits': int(kind in READS)}
    return maps, _cost(metadata, **extra)


def execute(program, *, backend='source', initial=None, stop_after=None,
            checkpoint=None, max_branches=256):
    """Execute or resume a complete CQ history without resetting its work state.

    ``initial`` is a normalized work vector or work-by-reference amplitude
    matrix. It explicitly replaces vacuum preparation and thus does not inherit
    its occupation certificate. Prefix Kraus maps act on this initial input;
    default preparation is counted once, before that instrument. All outcomes,
    including exact zero vectors, remain in lexical record order.
    """
    program = validate_program(program)
    if type(backend) is not str or backend not in BACKENDS:
        raise ValueError('a typed program and declared execution backend required')
    if type(max_branches) is not int or not 1 <= max_branches <= 2**MAX_READS:
        raise ValueError('an integer branch allocation budget from 1 through 256 required')
    end = len(program.operations) if stop_after is None else stop_after
    if type(end) is not int or not 0 <= end <= len(program.operations):
        raise ValueError('stop_after must be an operation-prefix index')
    reads = sum(op.kind in READS for op in program.operations[:end])
    if 2**reads > max_branches:
        raise ValueError('program prefix exceeds the explicit branch allocation budget')
    if checkpoint is not None:
        if initial is not None:
            raise ValueError('resumption cannot replace its retained quantum input')
        history = validate_checkpoint(program, checkpoint, backend)
        if end < history.cursor:
            raise ValueError('resumption cannot rewind the retained history')
    else:
        if initial is None:
            if backend == 'source':
                prepared, metadata = scalar.preparation(program.dimension)
                preparation = _cost(metadata, work_preparations=1)
            else:
                prepared = np.kron(field.vacuum_vector(program.dimension)[0], [1., 0., 0., 0.])
                preparation = _cost(work_preparations=1)
            initial_state, mode = _input_state(prepared, program.work_dimension), 'declared_vacuum'
        else:
            initial_state, preparation, mode = _input_state(initial, program.work_dimension), _cost(), 'supplied_input'
        kick, kick_cost = _kick(program.dimension, program.native_current, backend)
        resources = _sum_cost(preparation, kick_cost)
        branch = Branch((), _immutable(kick@initial_state), kick, resources)
        history = _checked(History(program, backend, 0, (branch,), initial_state,
                                   preparation, resources, (), mode))
    for cursor in range(history.cursor, end):
        operation = program.operations[cursor]
        output, tree_cost = [], history.expanded_tree_resources
        for branch in history.branches:
            selected = branch.record[operation.record] if operation.kind == 'feedback' else None
            maps, cost = _operation_maps(program, operation, selected, backend)
            # The retained tree accounts for each node once. It is distinct
            # from costs along one realized path, and from host matrix caching.
            tree_cost = _sum_cost(tree_cost, cost)
            for bit, matrix in enumerate(maps):
                record = (*branch.record, bit) if operation.kind in READS else branch.record
                output.append(Branch(record, _immutable(matrix@branch.state),
                                     _immutable(matrix@branch.kraus),
                                     _sum_cost(branch.resources, cost)))
        temporary = replace(history, cursor=cursor+1, branches=tuple(output),
                            expanded_tree_resources=tree_cost, seal='')
        stage = Stage(cursor+1, operation.kind, len(output[0].record), len(output),
                      temporary.completeness_error, temporary.probability_mass)
        history = _checked(replace(temporary, stages=(*history.stages, stage)))
    return history


def prefix_effects(history, length):
    """All-input effects of the retained old record, summing later outcomes."""
    count = len(history.branches[0].record)
    if type(length) is not int or not 0 <= length <= count:
        raise ValueError('an existing prefix length required')
    effects = {prefix: np.zeros((history.program.work_dimension,)*2, complex)
               for prefix in product((0, 1), repeat=length)}
    for branch in history.branches:
        effects[branch.record[:length]] += branch.kraus.conj().T@branch.kraus
    return effects


def execution_prefixes(program, *, backend='source', initial=None):
    """Execute once while retaining each operation checkpoint for stage audits."""
    program = validate_program(program)
    history = execute(program, backend=backend, initial=initial, stop_after=0)
    prefixes = [history]
    for cursor in range(1, len(program.operations)+1):
        history = execute(program, backend=backend, stop_after=cursor, checkpoint=history)
        prefixes.append(history)
    return tuple(prefixes)


def compare_histories(left, right):
    """Compare every full branch map and complete work/reference output."""
    if left.program != right.program or left.cursor != right.cursor or left.records != right.records:
        raise ValueError('comparison requires the same program and complete record prefix')
    if left.initial_state.shape != right.initial_state.shape:
        raise ValueError('comparison requires matching reference registers')
    return {
        'maximum_kraus_operator_error': max(float(np.linalg.norm(a.kraus-b.kraus, 2))
                                          for a, b in zip(left.branches, right.branches)),
        'complete_flagged_vector_error': math.sqrt(sum(float(np.linalg.norm(a.state-b.state)**2)
                                          for a, b in zip(left.branches, right.branches))),
        'classical_total_variation': math.fsum(abs(a.probability-b.probability)
                                          for a, b in zip(left.branches, right.branches))/2,
        'left_completeness_error': left.completeness_error,
        'right_completeness_error': right.completeness_error,
    }


def _encoded_array(array):
    flat = np.asarray(array).ravel()
    return {'shape': list(array.shape),
            'values': [[float(z.real), float(z.imag)] for z in flat]}


def checkpoint_to_json(history):
    """Canonical simulator checkpoint; no physical quantum transfer is implied."""
    history = validate_checkpoint(history.program, history, history.backend)
    payload = {'schema': 'oph.common_source.checkpoint.v1',
               'program': history.program.specification(), 'backend': history.backend,
               'cursor': history.cursor, 'preparation_mode': history.preparation_mode,
               'initial_state': _encoded_array(history.initial_state),
               'preparation_resources': history.preparation_resources,
               'expanded_tree_resources': history.expanded_tree_resources,
               'stages': [stage.__dict__ for stage in history.stages],
               'branches': [{'record': branch.record, 'state': _encoded_array(branch.state),
                             'kraus': _encoded_array(branch.kraus), 'resources': branch.resources}
                            for branch in history.branches], 'seal': history.seal}
    value = _canonical(payload)
    if len(value) > MAX_CHECKPOINT_BYTES:
        raise ValueError('checkpoint exceeds the explicit serialization byte budget')
    return value.decode()


def _keys(value, expected):
    if type(value) is not dict or set(value) != set(expected):
        raise ValueError('checkpoint object has missing or unexpected fields')


def _decode_cost(value):
    if (type(value) is not list or len(value) != len(RESOURCE_KEYS)
            or any(type(pair) is not list or len(pair) != 2 for pair in value)
            or tuple(pair[0] for pair in value) != RESOURCE_KEYS
            or any(type(pair[1]) is not int or pair[1] < 0 for pair in value)):
        raise ValueError('checkpoint resources must be complete nonnegative integer counts')
    return tuple((key, count) for key, count in value)


def _decode_array(value, rows, columns=None):
    _keys(value, ('shape', 'values'))
    shape, values = value['shape'], value['values']
    if (type(shape) is not list or len(shape) != 2
            or any(type(item) is not int or item < 1 for item in shape)
            or shape[0] != rows or (columns is not None and shape[1] != columns)
            or type(values) is not list or len(values) != shape[0]*shape[1]):
        raise ValueError('checkpoint array shape differs from its declared work/reference register')
    if any(type(pair) is not list or len(pair) != 2
           or any(type(x) not in (int, float) or not _finite_float(x) for x in pair)
           for pair in values):
        raise ValueError('checkpoint amplitudes must be finite real/imaginary numeric pairs')
    result = np.asarray(values, dtype=float)
    # Componentwise assignment preserves IEEE signed zeros. Arithmetic such
    # as real+1j*imag can change their bits and invalidate a valid byte seal.
    decoded = np.empty(len(result), dtype=np.complex128)
    decoded.real = result[:, 0]
    decoded.imag = result[:, 1]
    return _immutable(decoded.reshape(shape))


def checkpoint_from_json(program, payload, *, backend='source'):
    """Restore a strict, bounded JSON envelope for exactly this program/backend.

    No path or executable code is accepted. A public digest binds serialized
    simulator data, not independent experimental provenance or a proof that a
    source device executed it. Mathematical composition is checked separately.
    """
    program = validate_program(program)
    if type(payload) is not str:
        raise ValueError('a typed program and JSON text are required')
    if len(payload) > MAX_CHECKPOINT_BYTES or len(payload.encode()) > MAX_CHECKPOINT_BYTES:
        raise ValueError('checkpoint exceeds the explicit serialization byte budget')

    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate checkpoint object key')
            result[key] = value
        return result

    def reject_constant(value):
        raise ValueError(f'nonfinite JSON constant {value}')

    value = json.loads(payload, object_pairs_hook=pairs, parse_constant=reject_constant)
    _keys(value, ('schema', 'program', 'backend', 'cursor', 'preparation_mode', 'initial_state',
                  'preparation_resources', 'expanded_tree_resources', 'stages', 'branches', 'seal'))
    if (value['schema'] != 'oph.common_source.checkpoint.v1'
            or _canonical(value['program']) != _canonical(program.specification())
            or type(value['backend']) is not str or value['backend'] != backend
            or type(value['cursor']) is not int or not 0 <= value['cursor'] <= len(program.operations)
            or value['preparation_mode'] not in ('declared_vacuum', 'supplied_input')
            or type(value['seal']) is not str):
        raise ValueError('checkpoint specification, cursor, backend, or preparation is incompatible')
    initial = _decode_array(value['initial_state'], program.work_dimension)
    stages = value['stages']
    if type(stages) is not list or len(stages) != value['cursor']:
        raise ValueError('checkpoint must retain each executed stage')
    stage_rows = []
    read_count = 0
    for index, item in enumerate(stages):
        _keys(item, ('cursor', 'kind', 'record_count', 'branch_count', 'completeness_error', 'probability_mass'))
        operation = program.operations[index]
        read_count += operation.kind in READS
        if (type(item['cursor']) is not int or item['cursor'] != index+1
                or item['kind'] != operation.kind
                or type(item['record_count']) is not int or item['record_count'] != read_count
                or type(item['branch_count']) is not int or item['branch_count'] != 2**read_count
                or any(type(item[key]) not in (int, float) or not _finite_float(item[key])
                       for key in ('completeness_error', 'probability_mass'))):
            raise ValueError('checkpoint stage journal does not match the program prefix')
        stage_rows.append(Stage(**item))
    rows = value['branches']
    expected_records = tuple(product((0, 1), repeat=read_count))
    if type(rows) is not list or len(rows) != len(expected_records):
        raise ValueError('checkpoint omitted a complete recorded branch')
    branches = []
    for row, expected in zip(rows, expected_records):
        _keys(row, ('record', 'state', 'kraus', 'resources'))
        record = row['record']
        if (type(record) is not list or any(type(bit) is not int for bit in record)
                or tuple(record) != expected):
            raise ValueError('checkpoint records are not the complete ordered binary menu')
        branches.append(Branch(expected,
            _decode_array(row['state'], program.work_dimension, initial.shape[1]),
            _decode_array(row['kraus'], program.work_dimension, program.work_dimension),
            _decode_cost(row['resources'])))
    history = History(program, value['backend'], value['cursor'], tuple(branches), initial,
                      _decode_cost(value['preparation_resources']),
                      _decode_cost(value['expanded_tree_resources']), tuple(stage_rows),
                      value['preparation_mode'], value['seal'])
    history = validate_checkpoint(program, history, backend)
    # Imported JSON is a trust boundary: replay its declared prefix so a
    # consistently forged state/Kraus pair with a recomputed public digest is
    # not accepted as a checkpoint of that program. In-memory resumption keeps
    # the cheaper integrity check; import pays this once before continuing.
    replay = execute(program, backend=backend, stop_after=history.cursor,
                     initial=None if history.preparation_mode == 'declared_vacuum' else history.initial_state)
    comparison = compare_histories(history, replay)
    if (comparison['maximum_kraus_operator_error'] > TOLERANCE
            or comparison['complete_flagged_vector_error'] > TOLERANCE
            or np.linalg.norm(history.initial_state-replay.initial_state) > TOLERANCE
            or history.preparation_resources != replay.preparation_resources
            or history.expanded_tree_resources != replay.expanded_tree_resources
            or any(a.resources != b.resources for a, b in zip(history.branches, replay.branches))):
        raise ValueError('checkpoint does not reproduce the declared prefix and resource ledger')
    for saved, actual in zip(history.stages, replay.stages):
        if (abs(saved.completeness_error-actual.completeness_error) > TOLERANCE
                or abs(saved.probability_mass-actual.probability_mass) > TOLERANCE):
            raise ValueError('checkpoint stage diagnostics do not reproduce its prefix')
    return history
