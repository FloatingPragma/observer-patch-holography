"""Compile finite Petz recovery and observation services into native source words.

The reference, channel design and initialized work registers are declared.
Every first-failure Kraus branch is retained by its literal prefix circuit;
streamed Gram diagnostics reduce host storage without changing the channel.
"""
from __future__ import annotations
from dataclasses import dataclass
from types import MappingProxyType
import math
import numpy as np
from ..observer.source_scalar import compiler
from ..reference_observation.observation import ReferenceObservation, frozen
from ..field_thermal.thermal import ThermalField, car_restriction

TOLERANCE = 2e-8
MAX_PULSES = 1_000_000


def _array(value, shape=None):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked arrays are not admitted')
    a = np.asarray(value)
    if (a.dtype.kind not in 'iufc' or not np.isfinite(a).all()
            or (shape is not None and a.shape != shape)):
        raise ValueError('finite array on the complete declared domain required')
    return np.asarray(a, complex)


def compile_isometry(value):
    """V:D_in->D_out, at most four input columns and six output work bits.

    Gray-row QR eliminates only the admitted input columns. The resulting
    diagonal phases, initial Gray permutation and inverse elimination gates
    are all emitted. Other input columns see the same unitary completion.
    There is no target matrix primitive or success normalization.
    """
    v = _array(value)
    if v.ndim != 2:
        raise ValueError('rectangular isometry required')
    d, n = v.shape
    bits = compiler._bits(d, 6)
    if n not in (2, 4) or n > d:
        raise ValueError('two or four input columns required')
    if np.linalg.norm(v.conj().T@v-np.eye(n), 2) > 2e-10:
        raise ValueError('complete isometry required')
    gray = [j ^ (j >> 1) for j in range(d)]
    residual = v[gray].copy()
    real = not np.any(v.imag)
    moves = []
    for column in range(n):
        for row in range(d-1, column, -1):
            x, y = residual[row-1, column], residual[row, column]
            if real:
                delta = 0.
                theta = math.atan2(float(y.real), float(x.real)) if x or y else 0.
            else:
                delta = float(np.angle(y)-np.angle(x)) if x or y else 0.
                theta = math.atan2(abs(y), abs(x)) if x or y else 0.
                residual[row] *= np.exp(-1j*delta)
            c, s = math.cos(theta), math.sin(theta)
            a, b = residual[row-1].copy(), residual[row].copy()
            residual[row-1], residual[row] = c*a+s*b, -s*a+c*b
            moves.append((gray[row-1], gray[row], theta, delta))
    expected = np.zeros_like(residual)
    phases = np.zeros(d)
    for j in range(n):
        phases[gray[j]] = np.angle(residual[j, j])
        expected[j, j] = np.exp(1j*phases[gray[j]])
    if np.linalg.norm(residual-expected, 2) > 3e-10:
        raise ArithmeticError('isometry elimination exceeds the numerical design budget')
    word = [e for target in range(bits-1) for e in compiler.cnot(target+1, target)]
    word += compiler.compile_diagonal(phases)
    for a, b, theta, delta in reversed(moves):
        word += compiler._gray_rotation(a, b, -theta, bits)
        if not real:
            phases = np.zeros(d); phases[b] = delta
            word += compiler.compile_diagonal(phases)
    if len(word) > MAX_PULSES:
        raise ValueError('isometry source schedule exceeds the bounded pulse domain')
    return word


def _failure_vectors(state, full, targets, bits):
    d = 1 << len(targets)
    axes = [bits-1-j for j in targets]
    axes += [j for j in range(bits) if j not in axes]+[bits]
    arranged = state.reshape([2]*bits+[state.shape[1]]).transpose(axes)
    return (full[d:, :d]@arranged.reshape(d, -1)).reshape(6-d, -1, state.shape[1])


@dataclass(frozen=True)
class CompleteRun:
    """Literal complete output family with lazy first-failure branch access.

    The word, input and actual full primitive matrices are immutable snapshots.
    failure(event) reconstructs every complement Kraus vector for that first
    failure. Those Kraus vectors are summed incoherently at the common abort.
    The environment retained by the final success isometry is also explicit.
    """
    word: tuple
    bits: int
    initial: np.ndarray
    primitives: object
    success: np.ndarray
    failure_weights: np.ndarray
    failure_gram: np.ndarray

    def failure(self, event):
        if type(event) is not int or not 0 <= event < len(self.word):
            raise ValueError('existing first-failure event required')
        state = self.initial.copy()
        for i, (kind, angle, targets) in enumerate(self.word):
            full = self.primitives[kind, angle]
            if i == event:
                return tuple(frozen(x) for x in _failure_vectors(state, full, targets, self.bits))
            d = 1 << len(targets)
            state = compiler._apply(state, full[:d, :d], targets, self.bits)
        raise ArithmeticError('unreachable failure event')

    @property
    def failure_probability(self):
        return float(np.trace(self.failure_gram).real)


def execute_complete(word, bits, initial):
    """Execute all native pulses; stream norms but preserve the complete tape.

    This is host storage reduction only. The full branch family is defined by
    the frozen prefix matrices, including zero-weight and complement branches.
    """
    if type(bits) is not int or not 1 <= bits <= 9:
        raise ValueError('one through nine bounded work bits required')
    a = _array(initial)
    if a.ndim == 1:
        a = a[:, None]
    if a.ndim != 2 or a.shape[0] != 1 << bits or not 1 <= a.shape[1] <= 16:
        raise ValueError('complete work register and one through sixteen reference columns required')
    tape = compiler._validate(word, bits)
    if len(tape) > MAX_PULSES:
        raise ValueError('bounded complete source word required')
    word = tuple((e['kind'], e['angle'], tuple(e['targets'])) for e in tape)
    primitives = {}
    for kind, angle, _ in word:
        if (kind, angle) not in primitives:
            full = _array(compiler.native_pulse(kind, angle), (6, 6))
            if np.linalg.norm(full.conj().T@full-np.eye(6), 2) > 2e-11:
                raise ArithmeticError('complete source primitive is not unitary')
            primitives[kind, angle] = frozen(full)
    state = a.copy()
    weights = np.zeros((len(word), 4))  # J uses two complement slots; last two are padding.
    gram = np.zeros((a.shape[1], a.shape[1]), complex)
    for i, (kind, angle, targets) in enumerate(word):
        full = primitives[kind, angle]
        failed = _failure_vectors(state, full, targets, bits)
        for j, branch in enumerate(failed):
            gram += branch.conj().T@branch
            weights[i, j] = np.vdot(branch, branch).real
        d = 1 << len(targets)
        state = compiler._apply(state, full[:d, :d], targets, bits)
    diagnostics = (np.linalg.norm(state.conj().T@state+gram-a.conj().T@a, 2),
                   abs(float(np.sum(weights))-float(np.trace(gram).real)))
    if not np.isfinite(state).all() or not all(math.isfinite(x) for x in diagnostics):
        raise ArithmeticError('nonfinite complete output diagnostic')
    if max(diagnostics) > TOLERANCE*max(1., float(np.linalg.norm(a, 2))**2):
        raise ArithmeticError('complete output/reference Gram budget exceeded')
    return CompleteRun(word, bits, frozen(a), MappingProxyType(primitives), frozen(state),
                       frozen(weights), frozen(gram))


class RecoveryService:
    """Source realization of the adjoint of a bounded reference observation.

    The actual CAR state restriction determines the Kraus ordering. Input
    low bits carry the regional system; all other work bits start at zero.
    Output low bits carry the global field, followed by retained environment.
    """
    def __init__(self, reference, restriction):
        self.observation = ReferenceObservation(reference, restriction)
        self.kraus = tuple(frozen(k.conj().T) for k in self.observation.kraus)
        self.isometry = frozen(np.vstack(self.kraus))
        self.output_dimension = self.observation.source_dimension
        self.input_dimension = self.observation.target_dimension
        self.bits = compiler._bits(len(self.isometry), 6)
        self.word = compile_isometry(self.isometry)

    def execute(self, columns=None, *, pointer_mode=None):
        n = self.input_dimension
        columns = np.eye(n, dtype=complex) if columns is None else _array(columns)
        if columns.ndim != 2 or columns.shape[0] != n or not 1 <= columns.shape[1] <= 16:
            raise ValueError('complete regional input/reference columns required')
        bits, word = self.bits, self.word
        if pointer_mode is not None:
            modes = self.output_dimension.bit_length()-1
            if type(pointer_mode) is not int or not 0 <= pointer_mode < modes:
                raise ValueError('existing global field mode required')
            word = word+compiler.cnot(pointer_mode, bits); bits += 1
        initial = np.zeros((1 << bits, columns.shape[1]), complex)
        initial[:n] = columns
        run = execute_complete(word, bits, initial)
        target = self.isometry@columns
        if pointer_mode is not None:
            split = np.zeros_like(run.success)
            for address, row in enumerate(target):
                split[address+(((address >> pointer_mode) & 1) << self.bits)] = row
            target = split
        complete_error = math.sqrt(float(np.linalg.norm(run.success-target, 2))**2+
                                   max(0., float(np.linalg.norm(run.failure_gram, 2))))
        if not math.isfinite(complete_error) or complete_error > TOLERANCE*max(1., float(np.linalg.norm(columns, 2))):
            raise ArithmeticError('source recovery exceeds complete channel design budget')
        return run

    def resources(self):
        return {**compiler.resources(self.word), 'input_bits': self.input_dimension.bit_length()-1,
                'initialized_work_bits': self.bits-(self.input_dimension.bit_length()-1),
                'retained_environment_bits': self.bits-(self.output_dimension.bit_length()-1),
                'lazy_failure_family_retained': True, 'physical_duration_or_energy_calibrated': False}


def field_service(case=ThermalField(), kept=(0, 1)):
    target = case.target()
    return RecoveryService(target['reference'], car_restriction(case.geometry.modes, kept))


def nested_word(outer, inner):
    """Compose reference recoveries and retain both Stinespring environments."""
    if not isinstance(outer, RecoveryService) or not isinstance(inner, RecoveryService):
        raise ValueError('two compiled reference services required')
    if outer.input_dimension != inner.output_dimension:
        raise ValueError('nested recovery dimensions do not compose')
    if np.linalg.norm(outer.observation.target_reference-inner.observation.source_reference, 2) > 1e-12:
        raise ValueError('one identical intermediate reference required')
    middle_bits = outer.input_dimension.bit_length()-1
    extra = inner.bits-middle_bits
    bits = outer.bits+extra
    if bits > 9:
        raise ValueError('retained nested environments exceed the bounded work domain')
    targets = list(range(middle_bits))+list(range(outer.bits, bits))
    return compiler.remap(inner.word, targets)+outer.word, bits


def reference_preparation(service, *, pointer_mode=0):
    """Source-prepare the regional reference, recover, and read the global mode.

    The regional purifier and recovery environment both remain in the output.
    This supplies the same global Gibbs marginal by an independent source
    route, without assuming access to a pre-existing correlated global input.
    """
    if not isinstance(service, RecoveryService):
        raise ValueError('compiled reference service required')
    local_bits = service.input_dimension.bit_length()-1
    bits = service.bits+local_bits+1
    if bits > 9 or type(pointer_mode) is not int or not 0 <= pointer_mode < service.output_dimension.bit_length()-1:
        raise ValueError('bounded retained purifier and valid field pointer required')
    from ..reference_observation.observation import powers
    root, _ = powers(service.observation.target_reference)
    purification = root.T.reshape(-1)
    prep = compiler.compile_real_state(np.abs(purification))+compiler.compile_diagonal(np.angle(purification))
    targets = list(range(local_bits))+list(range(service.bits, service.bits+local_bits))
    word = compiler.remap(prep, targets)+service.word+compiler.cnot(pointer_mode, bits-1)
    initial = np.zeros((1 << bits, 1), complex); initial[0, 0] = 1
    run = execute_complete(word, bits, initial)
    # Compare with the complete mathematical isometry of that same purification.
    expected = np.zeros_like(initial)
    for pur in range(service.input_dimension):
        image = service.isometry@root[:, pur]
        for address, amplitude in enumerate(image):
            pointer = (address >> pointer_mode) & 1
            expected[address+(pur << service.bits)+(pointer << (bits-1)), 0] = amplitude
    error = math.sqrt(float(np.linalg.norm(run.success-expected))**2+run.failure_probability)
    if not math.isfinite(error) or error > TOLERANCE:
        raise ArithmeticError('complete prepared-reference observation exceeds admission budget')
    return run, {'complete_preparation_observation_error': error,
                 **compiler.resources(word), 'initialized_work_bits': bits,
                 'retained_regional_purifier_bits': local_bits,
                 'retained_recovery_environment_bits': service.bits-(service.output_dimension.bit_length()-1),
                 'all_pointer_and_export_outputs_retained': True,
                 'physical_initialization_duration_or_energy_calibrated': False}


def report():
    service = field_service()
    run = service.execute()
    prepared, resources = reference_preparation(service)
    p = prepared.success.reshape(2, -1)
    return {'scope': 'BOUNDED_SAME_REFERENCE_SOURCE_RECOVERY_AND_OBSERVATION',
            'field_sites': 2, 'field_modes': 4,
            'reference_beta': 1., 'chemical_potential': 0.,
            'recovery_isometry_operator_error': float(np.linalg.norm(run.success-service.isometry, 2)),
            'recovery_first_export_failure_weight_on_complete_input_basis': run.failure_probability,
            'recovery_resources': service.resources(),
            'prepared_reference_observation': resources,
            'complete_pointer_weights': [float(np.vdot(x, x).real) for x in p],
            'preparation_or_read_abort_probability': prepared.failure_probability,
            'generic_state_recovery_or_regional_retraction_claimed': False,
            'physical_clock_or_global_band_reference_supplied': False,
            'canonical_physical_premises_discharged': 0}


if __name__ == '__main__':
    import json
    print(json.dumps(report(), indent=2, allow_nan=False))
