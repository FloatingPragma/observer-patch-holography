"""Even reservoir spectroscopy of the existing charged dimer candidates.

The field generator and its supplied half-filled ground are unchanged.  A
supplied empty/occupied CAR reservoir and a neutral Ramsey pointer replace
controlled odd field probes.  Number and parity are conserved by completed
logical detector operations, not asserted for each encoded native pulse.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

import numpy as np
from scipy.linalg import expm

from ..field_thermal.thermal import execute_with_exports, failure_probability
from ..full_fock_field.field import FullFock
from ..full_fock_field.local_simulation import split_bound
from ..observer.source_scalar import compiler
from ..particle_response import dimer, source as spectroscopy
from ..particle_response.spectrum import FieldSpectrum
from ..reference_observation.observation import frozen
from . import source


MODES = 4
RESERVOIR = 4
POINTER = 5
BITS = 6
PRIMITIVE_GATE = 2e-9
COUNT_KEYS = ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
              'zero_angle_requests', 'identity_phase_pulses',
              'proper_code_load_events', 'proper_code_export_events',
              'proper_code_service_events')


def _index(index):
    if type(index) is not int or not 0 <= index < MODES:
        raise ValueError('an integer field mode from zero through three required')
    return index


def _bit(value, name):
    if type(value) is not int or value not in (0, 1):
        raise ValueError(f'{name} must be the integer zero or one')
    return value


def _admit(spectrum):
    if (not isinstance(spectrum, FieldSpectrum) or spectrum.dimension != 16
            or spectrum.geometry.volume != 2 or spectrum.geometry.speed != .125
            or spectrum.geometry.interaction != 1.
            or spectrum.reference_kind != 'ground' or spectrum.chemical_potential != .5):
        raise ValueError('the existing canonical pure half-filled dimer required')
    axis = spectrum.geometry.sides.index(2)
    if (np.linalg.norm(spectrum.reference@spectrum.reference-spectrum.reference, 2) > 2e-11
            or np.linalg.norm((spectrum.number-2*np.eye(16))@spectrum.reference, 2) > 2e-11):
        raise ArithmeticError('detector reference is not the declared pure number-two state')
    return axis


def _duration(duration, steps):
    if (type(duration) not in (int, F) or abs(duration) > 1
            or type(steps) is not int or not 1 <= steps <= 16):
        raise ValueError('exact duration of absolute value at most one and one through sixteen steps required')
    return F(duration), steps


def reflection(index=0, phase=0):
    """Literal full32 F_theta, with theta=phase*pi, including double occupation.

    F_theta = 1-n_j-n_r + exp(i theta)c_j* r + exp(-i theta)r* c_j.
    In particular the doubly occupied block is -1, not the bosonic SWAP.
    """
    index, phase = _index(index), _bit(phase, 'reflection phase')
    space = FullFock(5)
    a, r = (space.annihilator(j).toarray() for j in (index, RESERVOIR))
    value = (np.eye(32)-a.conj().T@a-r.conj().T@r
             +(-1)**phase*(a.conj().T@r+r.conj().T@a))
    return frozen(value)


def _controlled_rotation(label, angle):
    """Control an exp(i angle P) with its scalar phase retained."""
    base = list(label)+['I']*(BITS-len(label))
    coupled = base.copy()
    coupled[POINTER] = 'Z'
    return (compiler.pauli_rotation(base, angle/2)
            +compiler.pauli_rotation(coupled, -angle/2))


def _controlled_number_phase(mode, angle):
    identity = ['I']*BITS
    z = identity.copy(); z[mode] = 'Z'
    return (_controlled_rotation(identity, angle/2)
            +_controlled_rotation(z, -angle/2))


def controlled_reflection_word(index=0, phase=0):
    """Actual native pulses for C(F_theta), including every Jordan--Wigner sign.

    exp(-i*pi*(c*j r+r*j c)/2) exp(i*pi*(n_j+n_r)/2) is F_0.
    Reversing only the hopping angle gives F_pi.  The full pair operation
    commutes with total CAR number even though individual Pauli factors need
    not.  No target matrix is used in execution.
    """
    index, phase = _index(index), _bit(phase, 'reflection phase')
    word = []
    for letter in ('X', 'Y'):
        label = ['I']*BITS
        label[index] = label[RESERVOIR] = letter
        for middle in range(index+1, RESERVOIR):
            label[middle] = 'Z'
        word += _controlled_rotation(label, -(-1)**phase*math.pi/4)
    word += _controlled_number_phase(index, math.pi/2)
    word += _controlled_number_phase(RESERVOIR, math.pi/2)
    return word


def _validate_failures(word, bits, failures, references=1, vector=True):
    """Zero-probability complements remain part of the complete instrument."""
    if type(failures) not in (tuple, list) or len(failures) != len(word):
        raise ArithmeticError('complete first-export failure event grid required')
    for event, branches in zip(word, failures, strict=True):
        dimension = 1 << len(event['targets'])
        if type(branches) not in (tuple, list) or len(branches) != 6-dimension:
            raise ArithmeticError('complete labeled export complement required')
        shape = ((1 << bits)//dimension,) if vector else ((1 << bits)//dimension, references)
        for branch in branches:
            if (np.ma.isMaskedArray(branch) or np.asarray(branch).shape != shape
                    or np.asarray(branch).dtype.kind not in 'iufc' or not np.isfinite(branch).all()):
                raise ArithmeticError('finite unnormalized full failure/reference columns required')


def _counts(word):
    values = compiler.resources(word)
    return {key: values[key] for key in COUNT_KEYS}


def operator_certificate():
    """Full operator identity, hence valid with arbitrary untouched references.

    Reservoir one gives C_theta=<n(t)n>+exp(-i theta)<a(t)a*>.
    Reservoir zero gives C_theta=<(1-n(t))(1-n)>+exp(i theta)<a*(t)a>.
    The last charged term is the conjugate of the current reverse-Green
    convention.  No stationary-state or factorization assumption on the
    field/reference input is used; the supplied reservoir is a product state.
    """
    space = FullFock(5)
    number = sum(space.number(j).toarray() for j in range(5))
    parity = np.diag((-1.)**np.diag(number))
    columns, identities = 0, 0
    for index in range(4):
        for phase in (0, 1):
            f = reflection(index, phase)
            if not all(np.array_equal(x, y) for x, y in (
                    (f.conj().T, f), (f@f, np.eye(32)),
                    (f@number, number@f), (f@parity, parity@f))):
                raise ArithmeticError('reservoir reflection lost its exact complete CAR identities')
            word = controlled_reflection_word(index, phase)
            actual, failures = execute_with_exports(word, BITS, np.eye(64))
            _validate_failures(word, BITS, failures, 64, False)
            expected = np.zeros((64, 64), complex)
            expected[:32, :32] = np.eye(32)
            expected[32:, 32:] = f
            error = math.sqrt(float(np.linalg.norm(actual-expected))**2+failure_probability(failures))
            if not math.isfinite(error) or error > PRIMITIVE_GATE:
                raise ArithmeticError('executed controlled reflection differs on complete input columns')
            columns += 64
    # Three actual spatial orientations, with an independent exponential.
    for axis in range(3):
        sides = [1, 1, 1]; sides[axis] = 2
        from ..full_fock_field.field import FieldGeometry
        spectrum = dimer.setup(FieldGeometry(tuple(sides)))
        unitary = expm(-.25j*spectrum.generator)
        lifted = np.kron(np.eye(2), unitary)
        for index in range(4):
            a = spectrum.space.annihilator(index).toarray()
            evolved = unitary.conj().T@a@unitary
            n = a.conj().T@a; nt = evolved.conj().T@evolved
            for reservoir in (0, 1):
                for phase in (0, 1):
                    actual = lifted.conj().T@reflection(index, phase)@lifted@reflection(index, 0)
                    block = actual[16*reservoir:16*(reservoir+1), 16*reservoir:16*(reservoir+1)]
                    expected = (nt@n+(-1)**phase*evolved@a.conj().T if reservoir else
                                (np.eye(16)-nt)@(np.eye(16)-n)+(-1)**phase*evolved.conj().T@a)
                    if np.linalg.norm(block-expected, 2) > 2e-11:
                        raise ArithmeticError('arbitrary-input reservoir Green identity failed')
                    identities += 1
    return dict(full_controlled_input_columns_checked=columns,
                complete_field_operator_identities_checked=identities,
                complete_reflection_native_word_gate='2e-9',
                exact_reflection_number_and_parity_conservation=True,
                arbitrary_untouched_reference_extension=True,
                supplied_reservoir_product_state_required=True,
                no_stationarity_assumption_for_operator_identity=True,
                conservation_applies_to_completed_logical_operations=True,
                primitive_pulse_CAR_number_conservation_claimed=False)


def correlator(spectrum, index=0, reservoir=1, phase=0, quadrature=0,
               duration=F(1, 4), steps=4):
    """A fresh complete Ramsey setting; public callers cannot replace preparation."""
    axis = _admit(spectrum)
    index = _index(index)
    reservoir, phase, quadrature = (_bit(value, name) for value, name in
        ((reservoir, 'reservoir occupation'), (phase, 'reflection phase'), (quadrature, 'quadrature')))
    duration, steps = _duration(duration, steps)
    prepared = source.prepare('ground', axis=axis)
    _validate_failures(prepared.word, MODES, prepared.failures)
    if (np.asarray(prepared.state).shape != (16,) or np.asarray(prepared.target).shape != (16,)
            or not np.isfinite(prepared.state).all() or not np.isfinite(prepared.target).all()
            or np.linalg.norm(np.outer(prepared.target, prepared.target.conj())-spectrum.reference, 2) > 2e-11
            or any(prepared.resources[key] != value for key, value in _counts(prepared.word).items())):
        raise ArithmeticError('same ground preparation or its actual resource account changed')
    initial = np.zeros(64, complex)
    initial[16*reservoir:16*(reservoir+1)] = prepared.state
    hadamard = compiler.remap(compiler.compile_real_orthogonal(spectroscopy.HADAMARD), [POINTER])
    forward = spectroscopy.forward_word(spectrum, duration, steps)
    pointer_phase = -quadrature*math.pi/2
    word = (hadamard+controlled_reflection_word(index, 0)+forward
            +controlled_reflection_word(index, phase)+compiler.inverse(forward)
            +spectroscopy._pointer_phase(POINTER, pointer_phase)+hadamard)
    actual, failures = execute_with_exports(word, BITS, initial)
    _validate_failures(word, BITS, failures)
    outputs = tuple(frozen(row) for row in actual.reshape(2, 32))
    probabilities = [float(np.vdot(row, row).real) for row in outputs]
    prefix_abort = failure_probability(prepared.failures)
    continuation_abort = failure_probability(failures)
    complete = sum(probabilities)+prefix_abort+continuation_abort
    target = np.zeros(32, complex)
    target[16*reservoir:16*(reservoir+1)] = prepared.target
    unitary = np.kron(np.eye(2), expm(-1j*float(duration)*spectrum.generator))
    operation = unitary.conj().T@reflection(index, phase)@unitary@reflection(index, 0)
    changed = operation@target
    expected = ((target+np.exp(1j*pointer_phase)*changed)/2,
                (target-np.exp(1j*pointer_phase)*changed)/2)
    error = math.sqrt(sum(float(np.linalg.norm(a-b))**2 for a, b in zip(outputs, expected, strict=True))
                      +prefix_abort+continuation_abort)
    eta = 2*split_bound(spectrum.geometry.sides, duration, steps,
                       F(spectrum.geometry.speed), F(spectrum.geometry.interaction), refined=True)
    correlation = complex(np.vdot(target, changed))
    expected_mean = float((np.exp(1j*pointer_phase)*correlation).real)
    forbidden = [mask for mask in range(32) if mask.bit_count() != 2+reservoir]
    leakage = sum(float(np.vdot(row[forbidden], row[forbidden]).real) for row in outputs)
    if (not all(math.isfinite(x) for x in (complete, error, leakage, expected_mean))
            or abs(complete-1) > PRIMITIVE_GATE or error > float(eta)+4e-9 or leakage > 1e-16):
        raise ArithmeticError('complete number-preserving source detector failed admission')
    return dict(pointer_outputs=outputs, first_export_failures=failures,
                preparation_first_export_failures=prepared.failures,
                pointer_probabilities=probabilities, preparation_abort_probability=prefix_abort,
                continuation_abort_probability=continuation_abort, complete_probability=complete,
                measured_mean=probabilities[0]-probabilities[1], ideal_mean=expected_mean,
                ideal_correlation=correlation, complete_vector_error=error,
                ideal_operator_error_bound=str(eta), number_sector_leakage=leakage,
                resources=_counts(word), preparation_resources=_counts(prepared.word),
                first_export_event_count=len(word), preparation_export_event_count=len(prepared.word),
                retained_continuation_complement_count=sum(len(row) for row in failures),
                retained_preparation_complement_count=sum(len(row) for row in prepared.failures),
                success_or_particle_postselection=False)


def readout(spectrum=None, index=0, duration=F(1, 4), steps=4):
    """Eight complete settings reconstruct both current charged Green functions."""
    spectrum = dimer.setup() if spectrum is None else spectrum
    axis = _admit(spectrum); index = _index(index)
    duration, steps = _duration(duration, steps)
    correlations, records = {}, []
    totals = dict.fromkeys(COUNT_KEYS, 0)
    worst_error = 0.
    for reservoir in (0, 1):
        for phase in (0, 1):
            values = []
            for quadrature in (0, 1):
                row = correlator(spectrum, index, reservoir, phase, quadrature, duration, steps)
                values.append(row['measured_mean'])
                worst_error = max(worst_error, row['complete_vector_error'])
                for key in COUNT_KEYS:
                    totals[key] += row['resources'][key]+row['preparation_resources'][key]
                records.append(dict(reservoir_occupation=reservoir, reflection_phase_in_pi=phase,
                    pointer_quadrature=quadrature,
                    **{key: row[key] for key in ('pointer_probabilities', 'preparation_abort_probability',
                        'continuation_abort_probability', 'complete_probability', 'resources', 'preparation_resources',
                        'first_export_event_count', 'preparation_export_event_count',
                        'retained_continuation_complement_count', 'retained_preparation_complement_count')}))
            correlations[reservoir, phase] = complex(*values)
    addition = (correlations[1, 0]-correlations[1, 1])/2
    removal = ((correlations[0, 0]-correlations[0, 1])/2).conjugate()
    creation, charge = spectrum.probe('single', index)
    measure = spectrum.measure(creation, charge)
    projectors = dimer.candidate_projectors(spectrum)
    expected_residue = (1+1/math.sqrt(17))/4
    candidate_residues = {}
    for kind, operation in (('charge_plus_one_particle', creation),
                            ('charge_minus_one_hole', creation.conj().T)):
        residue = float(np.trace(projectors[kind]@operation@spectrum.reference@operation.conj().T).real)
        if not math.isfinite(residue) or abs(residue-expected_residue) > 2e-11:
            raise ArithmeticError('even detector probe lost the existing charged-candidate residue')
        candidate_residues[kind] = residue
    expected = tuple(measure.correlation(float(duration), direction) for direction in ('forward', 'reverse'))
    residual = max(abs(addition-expected[0]), abs(removal-expected[1]))
    if not math.isfinite(residual) or residual > 4*worst_error+4e-9:
        raise ArithmeticError('even detector disagrees with the complete same-candidate Lehmann response')
    encode = lambda z: [float(z.real), float(z.imag)]
    return dict(scope='SAME_DIMER_EVEN_RESERVOIR_SINGLE_PARTICLE_RESPONSE', axis=axis,
                index=index, duration=str(duration), split_steps=steps,
                reference='same supplied pure N=2 ground of K=H-N/2',
                generator_and_candidate_projectors_changed=False,
                observed_forward=encode(addition), observed_reverse=encode(removal),
                complete_Lehmann_forward=encode(expected[0]), complete_Lehmann_reverse=encode(expected[1]),
                same_lower_charged_candidate_probe_residues=candidate_residues,
                lower_candidate_gap_formula='sqrt(17)/8-1/16',
                lower_candidate_residue_formula='(1+1/sqrt(17))/4',
                maximum_correlation_residual=float(residual), maximum_complete_vector_error=worst_error,
                ideal_per_setting_operator_error_bound=row['ideal_operator_error_bound'],
                complete_records=records, independent_preparations_required=8,
                total_resources_including_preparations=totals,
                initialized_field_bits=4, supplied_initial_field_number=2,
                initialized_field_registers_across_settings=8,
                initialized_reservoir_modes_per_setting=1, reservoir_occupations=[0, 1],
                initialized_reservoir_modes_across_settings=8,
                supplied_reservoir_particles_across_settings=4,
                initialized_neutral_pointer_bits_per_setting=1, maximum_work_bits=6,
                initialized_neutral_pointer_bits_across_settings=8,
                retained_purifier_bits=0, pointer_outcomes_per_setting=2,
                reservoir_retained_in_every_pointer_output=True,
                subtraction_is_linear_record_processing=True,
                conditioned_on_success_or_particle_outcome=False,
                controlled_odd_field_operators_used=False,
                completed_logical_detector_preserves_total_CAR_number=True,
                native_primitive_number_conservation_claimed=False,
                physical_electromagnetic_reservoir_or_detector_derived=False,
                physical_clock_or_particle_identity_derived=False)


def report():
    """Small deterministic receipt: execute first, then round display observations."""
    algebra = operator_certificate()
    observed = readout()
    observed.pop('maximum_correlation_residual')
    observed.pop('maximum_complete_vector_error')
    observed['numerical_admission'] = {
        'complete_probability_gate': '2e-9',
        'complete_vector_gate': 'twice exact split bound plus 4e-9',
        'output_wrong_number_probability_gate': '1e-16',
        'floating_checks_are_interval_or_hardware_certificates': False}
    def stable(value):
        if type(value) is dict:
            return {key: stable(item) for key, item in value.items()}
        if type(value) in (list, tuple):
            return [stable(item) for item in value]
        if type(value) is float:
            result = round(value, 10)
            return 0. if result == 0 else result
        return value
    return {'operator_contract': algebra, 'complete_source_readout': stable(observed)}


certificate = report
