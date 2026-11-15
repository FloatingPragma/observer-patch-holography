"""Fixed-reference modular time on the actual common working algebra.

The Gibbs reference is a declared mathematical reference built from the
existing finite Hamiltonian. It is never substituted for the observer's
measured state. For this reference sigma_s(A)=tau_H,-beta*s(A); source words
approximate the same evolution on every retained terminal branch.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import math

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar.experiment import encoded
from . import experiment as e, operational as o, scalar_coupling as scalar


RECORDS = ((0, 0), (0, 1), (1, 0), (1, 1))
RECORD_WEIGHTS = (F(1, 10), F(2, 10), F(3, 10), F(4, 10))
TOLERANCE = 2e-9


@dataclass(frozen=True)
class ThermalReference:
    dimension: int
    beta: float
    coupling: float
    hamiltonian: np.ndarray
    density: np.ndarray
    energies: np.ndarray
    basis: np.ndarray
    probabilities: np.ndarray
    log_probabilities: np.ndarray
    log_partition: float


def reference(d, beta=1., coupling=.25):
    """Stable Gibbs spectrum of the actual finite H, with no kernel filling."""
    e._register(d)
    beta = scalar._real(beta, 'positive dimensionless inverse temperature')
    coupling = scalar._real(coupling, 'coupling')
    if beta <= 0:
        raise ValueError('a positive dimensionless inverse temperature is required')
    h = scalar.ideal_hamiltonian(d, coupling)
    energies, basis = np.linalg.eigh(h)
    with np.errstate(over='ignore', under='ignore'):
        shifted = -beta*(energies-energies[0])
        weights = np.exp(shifted)
    log_shifted_partition = math.log(float(sum(weights)))
    log_p = shifted-log_shifted_partition
    with np.errstate(under='ignore'):
        p = np.exp(log_p)
    if not np.isfinite(log_p).all() or not np.all(p > 0):
        raise ValueError('Gibbs spectrum is faithful analytically but not representable at this precision')
    rho = (basis*p)@basis.conj().T
    if np.linalg.eigvalsh(rho)[0] <= 0:
        raise ValueError('numerical Gibbs matrix lost faithfulness; no artificial floor is inserted')
    log_z = log_shifted_partition-beta*float(energies[0])
    for item in (h, energies, basis, p, log_p, rho):
        item.flags.writeable = False
    return ThermalReference(d, beta, coupling, h, rho, energies, basis, p, log_p, log_z)


def modular_unitary(ref, parameter):
    parameter = scalar._real(parameter, 'modular parameter')
    return (ref.basis*np.exp(1j*parameter*ref.log_probabilities))@ref.basis.conj().T


def modular_flow(ref, observable, parameter):
    observable = np.asarray(observable, dtype=complex)
    if observable.shape != ref.density.shape or not np.isfinite(observable).all():
        raise ValueError('a finite observable on the full working algebra required')
    u = modular_unitary(ref, parameter)
    return u@observable@u.conj().T


def _record_weights(weights):
    if not isinstance(weights, (tuple, list)) or len(weights) != 4:
        raise ValueError('four strictly positive normalized reference record weights required')
    numeric = tuple(scalar._real(x, 'reference record weight') for x in weights)
    if min(numeric) <= 0 or abs(math.fsum(numeric)-1) > 1e-12:
        raise ValueError('four strictly positive normalized reference record weights required')
    return numeric


def block_modular_flow(ref, observables, parameter, weights=RECORD_WEIGHTS):
    """Flow on direct-sum record blocks, without allocating a fourfold matrix."""
    weights = _record_weights(weights)
    parameter = scalar._real(parameter, 'modular parameter')
    if set(observables) != set(RECORDS):
        raise ValueError('all four committed record blocks are required')
    result = {}
    for record, weight in zip(RECORDS, weights, strict=True):
        observable = np.asarray(observables[record], dtype=complex)
        if observable.shape != ref.density.shape or not np.isfinite(observable).all():
            raise ValueError('a finite observable on each full working block required')
        # log(weight) is central. Retain it to check its exact cancellation.
        u = (ref.basis*np.exp(1j*parameter*(ref.log_probabilities+math.log(weight))))@ref.basis.conj().T
        result[record] = u@observable@u.conj().T
    return result


def _matrix(matrix):
    return [encoded(row) for row in matrix]


def _expectation(vectors, observable):
    return math.fsum(float(np.vdot(v, observable@v).real) for v in vectors.values())


def _selection_audit(ref):
    size = len(ref.density)
    identity = np.eye(size)
    tracial = identity/size
    centered_h = ref.hamiltonian-np.trace(ref.hamiltonian).real*identity/size
    k = (ref.basis*(-ref.log_probabilities))@ref.basis.conj().T
    reconstructed_h = (k-ref.log_partition*identity)/ref.beta
    probe = e.projector(size, 0, 1)
    trace_flow_error = float(np.linalg.norm(o.modular_flow(tracial, probe, .5)-probe, 2))
    pinched = sum(p@ref.density@p for p in (probe, identity-probe))
    control = e._independent_feedback(size, 1, True)
    return {'constructed_from_supplied_hamiltonian': True,
            'generated_by_existing_repair': False,
            'reference_selection_proved': False,
            'gibbs_identity_is_a_new_dynamical_selection_theorem': False,
            'reference_first_equivalent_definition': 'H=(-log(sigma)-log(Z)*I)/beta',
            'reference_first_hamiltonian_reconstruction_error': float(np.linalg.norm(reconstructed_h-ref.hamiltonian, 2)),
            'tracial_reference_modular_flow_error_from_identity': trace_flow_error,
            'hamiltonian_spectral_spread': float(ref.energies[-1]-ref.energies[0]),
            'centered_hamiltonian_hilbert_schmidt_norm': float(np.linalg.norm(centered_h)),
            'gibbs_minus_tracial_density_operator_norm': float(np.linalg.norm(ref.density-tracial, 2)),
            'gauge_record_pinching_reference_displacement_operator_norm': float(np.linalg.norm(pinched-ref.density, 2)),
            'conditional_gauge_control_reference_displacement_operator_norm': float(np.linalg.norm(control@ref.density@control.conj().T-ref.density, 2)),
            'tracial_flow_can_match_this_nontrivial_H_at_nonzero_time_rescaling': False,
            'equivalent_reference_first_choice_removes_physical_input': False}


def receipt(result, beta=1., record_weights=RECORD_WEIGHTS):
    """Attach a fixed Gibbs reference and source time to the same full history."""
    weights = _record_weights(record_weights)
    ref = reference(result['dimension'], beta, result['coupling'])
    vectors = o.terminal_vectors(result)
    size = len(ref.density)
    identity = np.eye(size)
    h = ref.hamiltonian
    omega = sum(np.outer(v, v.conj()) for v in vectors.values())
    probe = e.projector(size, 0, 1)
    pointer, pointer_meta = e.pointer(result['dimension'], 0)
    settings = []
    for time in (0., .5, 1.):
        s = -time/ref.beta
        exact_schrodinger = expm(-1j*time*h)
        exact_heisenberg = exact_schrodinger.conj().T
        modular = modular_unitary(ref, s)
        # The central Gibbs normalizing phase has no observable effect.
        identity_error = float(np.linalg.norm(
            np.exp(1j*s*ref.log_partition)*modular-exact_heisenberg, 2))
        moved_probe = modular_flow(ref, probe, s)
        direct_modular_probe = o.modular_flow(ref.density, probe, s)
        op_error = float(np.linalg.norm(direct_modular_probe-moved_probe, 2))
        ideal_expectation = _expectation(vectors, moved_probe)
        exact_state = exact_schrodinger@omega@exact_schrodinger.conj().T
        stationary_error = float(np.linalg.norm(
            exact_schrodinger@ref.density@exact_schrodinger.conj().T-ref.density, 2))
        if time:
            executed, metadata = scalar.source_unitary(result['dimension'], result['steps_per_segment'],
                                                       time=time, coupling=result['coupling'])
            independent_split = scalar.independent_split_unitary(result['dimension'], result['steps_per_segment'],
                                                                  time=time, coupling=result['coupling'])
            certificate = scalar.time_error_bound(result['dimension'], result['steps_per_segment'],
                                                 F(time), F(result['coupling']))
            ideal_bound = float(F(certificate['ideal_split_operator_error_upper']))
            resources = {key: metadata['resources'][key]+pointer_meta['resources'][key] for key in e.COUNT_KEYS}
            execution_error = float(np.linalg.norm(executed-independent_split, 2))
        else:
            executed = identity
            independent_split = identity
            metadata = {'no_appended_evolution': True}
            certificate = {'ideal_split_operator_error_upper': '0',
                           'ideal_channel_half_diamond_error_upper': '0'}
            ideal_bound, execution_error = 0., 0.
            resources = {key: pointer_meta['resources'][key] for key in e.COUNT_KEYS}
        rows, marginal = [], [0., 0.]
        old_error = 0.
        for record in RECORDS:
            old = vectors[record]
            before_pointer = executed@old
            subtotal = 0.
            for bit, k in enumerate(pointer):
                state = k@before_pointer
                probability = float(np.vdot(state, state).real)
                subtotal += probability
                marginal[bit] += probability
                rows.append({'old_record': list(record), 'clock_bit': bit,
                             'probability': probability, 'unnormalized_terminal_state': encoded(state)})
            old_error = max(old_error, abs(subtotal-float(np.vdot(old, old).real)))
        record_blocks = {record: probe for record in RECORDS}
        block_flow = block_modular_flow(ref, record_blocks, s, record_weights)
        weighted_block_error = max(float(np.linalg.norm(block_flow[record]-moved_probe, 2)) for record in RECORDS)
        record_projector_error = 0.
        for selected in RECORDS:
            central = {record: identity if record == selected else np.zeros_like(identity) for record in RECORDS}
            after = block_modular_flow(ref, central, s, record_weights)
            record_projector_error = max(record_projector_error,
                                        max(float(np.linalg.norm(after[r]-central[r], 2)) for r in RECORDS))
        target_error = float(np.linalg.norm(executed-exact_schrodinger, 2))
        settings.append({'model_coordinate_time': time, 'modular_parameter': s,
                         'full_implementer_identity_error_after_central_phase': identity_error,
                         'all_observable_identity_error_bound_per_unit_operator_norm': 2*identity_error,
                         'existing_operational_modular_flow_probe_error': op_error,
                         'reference_stationarity_operator_error': stationary_error,
                         'actual_state_displacement_hilbert_schmidt': float(np.linalg.norm(exact_state-omega)),
                         'reference_gauge_one_probability': float(np.trace(ref.density@moved_probe).real),
                         'actual_ideal_gauge_one_probability': ideal_expectation,
                         'source_clock_bit_marginal': marginal,
                         'source_probability_error_from_modular_target': abs(marginal[1]-ideal_expectation),
                         'old_record_marginal_error': old_error,
                         'record_reference_weights_cancel_operator_error': weighted_block_error,
                         'central_record_projectors_fixed_operator_error': record_projector_error,
                         'branches': rows,
                         'source_evolution': metadata, 'source_pointer': pointer_meta,
                         'appended_source_resources_per_realized_history': {
                             **resources, 'fresh_pointer_blanks': 1, 'new_committed_classical_bits': 1},
                         'error_budget': {
                             'ideal_temporal_certificate': certificate,
                             'source_versus_exact_modular_schrodinger_operator_error': target_error,
                             'independent_split_versus_exact_operator_error': float(np.linalg.norm(independent_split-exact_schrodinger, 2)),
                             'source_minus_independent_split_operator_error': execution_error,
                             'ideal_bound_plus_measured_execution_defect': ideal_bound+execution_error,
                             'measured_execution_defect_is_a_rigorous_roundoff_certificate': False,
                             'physical_service_or_SI_clock_error_included': False}})
    actual_probabilities = [row['actual_ideal_gauge_one_probability'] for row in settings]
    source_probabilities = [row['source_clock_bit_marginal'][1] for row in settings]
    actual_rate = float(np.linalg.norm(h@omega-omega@h))
    selection = _selection_audit(ref)
    checks = {
        'fixed_gibbs_reference_is_faithful': bool(np.all(ref.probabilities > 0)),
        'full_modular_model_time_identity': max(row['full_implementer_identity_error_after_central_phase'] for row in settings) < TOLERANCE,
        'existing_modular_convention_agrees': max(row['existing_operational_modular_flow_probe_error'] for row in settings) < TOLERANCE,
        'reference_state_stationary': max(row['reference_stationarity_operator_error'] for row in settings) < TOLERANCE,
        'actual_observation_state_is_not_stationary': actual_rate > 1e-5,
        'actual_detector_probabilities_change': max(actual_probabilities)-min(actual_probabilities) > 1e-5,
        'source_detector_probabilities_change': max(source_probabilities)-min(source_probabilities) > 1e-5,
        'both_committed_records_preserved': max(max(row['old_record_marginal_error'], row['central_record_projectors_fixed_operator_error']) for row in settings) < TOLERANCE,
        'reference_weights_are_central': max(row['record_reference_weights_cancel_operator_error'] for row in settings) < TOLERANCE,
        'source_error_accounted_for_separately': all(row['error_budget']['source_versus_exact_modular_schrodinger_operator_error'] <= row['error_budget']['ideal_bound_plus_measured_execution_defect']+TOLERANCE for row in settings),
        'reference_first_equivalence': selection['reference_first_hamiltonian_reconstruction_error'] < TOLERANCE,
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    return {'schema': 'oph.common_source.fixed_thermal_reference_clock.v1',
            'dimension': result['dimension'], 'steps_per_segment': result['steps_per_segment'],
            'native_current': result['native_current'], 'coupling': result['coupling'],
            'beta': ref.beta, 'reference_density': _matrix(ref.density),
            'reference_energy_spectrum': ref.energies.tolist(),
            'reference_density_spectrum': ref.probabilities.tolist(),
            'reference_log_density_spectrum': ref.log_probabilities.tolist(),
            'reference_log_partition': ref.log_partition,
            'actual_complete_work_density': _matrix(omega),
            'actual_state_commutator_with_H_hilbert_schmidt': actual_rate,
            'actual_history_branches': result['branches'],
            'record_reference_weights': [str(weight) for weight in record_weights],
            'accessible_algebra': f'direct sum of four M_{size} full working blocks',
            'reference_definition': 'sigma_beta=exp(-beta H)/Z; each reference record block is p_ab sigma_beta',
            'convention': 'sigma_s(A)=sigma_beta^(is) A sigma_beta^(-is)=exp(i t H) A exp(-i t H), t=-beta s',
            'actual_state_definition': 'Omega_ab=|v_ab><v_ab| from the retained source history; not the Gibbs reference',
            'settings': settings, 'reference_selection_audit': selection,
            'checks': checks, 'all_pass': True,
            'scope': {'beta_and_reference_are_declared_construction_inputs': True,
                      'reference_prepared_or_substituted_for_actual_state': False,
                      'full_working_algebra_access_supplied': True,
                      'fixed_reference_is_updated_from_changing_actual_state': False,
                      'earlier_reduced_state_clock_or_negative_result_changed': False,
                      'identity_covers_unitary_segments_not_measurement_feedback_maps': True,
                      'fresh_complete_history_per_setting': True,
                      'monotone_unique_clock_readout_or_single_shot_time_claimed': False,
                      'physical_time_or_temperature_calibrated': False,
                      'repair_selection_of_thermal_reference_proved': False,
                      'all_existing_repair_channels_preserve_gibbs_reference': False,
                      'common_reference_Lean_consensus_tower_inhabited': False,
                      'reference_spectral_values_have_outward_interval_certificate': False}}


def report(steps=512):
    rows = []
    for d in (2, 4):
        history = e.history(d, steps)
        rows.append({'dimension': d, 'actual_history': history, 'thermal_clock': receipt(history)})
    return {'schema': 'oph.common_source.thermal_clock.v1', 'levels': rows, 'all_pass': True}
