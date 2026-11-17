"""Complete native detector execution and preparation-independent predictions."""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..anharmonic_observations import source as primitives
from ..anharmonic_transfer import calibration
from ..empirical_diagnostics import experiment as preparation
from ..observer.source_scalar import compiler
from ..primitive_observer import anharmonic, runtime
from ..scalar_observations.experimental_readout import executable_source as native
from . import bounds


FOLDER = Path(__file__).resolve().parent
ROOT = FOLDER.parents[1]


def channel_matrix(flip_zero, flip_one):
    values = [bounds._finite(x, 'flip probability') for x in (flip_zero, flip_one)]
    if min(values) < 0 or max(values) > 1:
        raise ValueError('flip probabilities in [0,1] required')
    return np.array([[1-values[0], values[1]], [values[0], 1-values[1]]])


def compile_detectors(dimension, settings, h):
    operations, metadata = {}, {}
    for setting in settings:
        name = setting['name']
        if name in operations:
            raise ValueError('unique setting names required')
        operations[name], metadata[name] = primitives.detector(
            dimension, 0., setting['multiplier']*h, setting['phase_half_pi']*math.pi/2)
    return operations, metadata


def complete_raw_operators(dimension, operations):
    spin, _ = anharmonic._spin_record(dimension)
    return {name: tuple(tuple(np.kron(np.eye(2), detector[bit]) @ spin[choice]
                             for choice in (0, 1)) for bit in (0, 1))
            for name, detector in operations.items()}


def classical_channel_law(rho, spin, detector, detector_channel):
    """Independently contract every spin/true-bit outcome, then confuse bits."""
    ideal = np.empty((2, 2))
    for choice in (0, 1):
        for bit in (0, 1):
            operator = np.kron(np.eye(2), detector[bit]) @ spin[choice]
            ideal[choice, bit] = float(np.trace(operator@rho@operator.conj().T).real)
    return ideal @ detector_channel.T


def complete_law_residual(actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if actual.shape != (2, 2) or expected.shape != (2, 2):
        raise ValueError('complete 2x2 spin/bit laws required')
    return float(np.max(abs(actual-expected)))


def protocol_sampling_designs(protocol):
    design = protocol['sampling_design']
    errors = protocol['declared_detector_channel']
    return [bounds.sampling_design(n, settings=design['number_of_settings'],
                failure_probability=design['failure_probability'],
                gain=1-errors['flip_zero_to_one']-errors['flip_one_to_zero'])
            for n in design['fresh_trials_per_setting']]


def execute(dimension, alpha, thermal, operations, detector_channel):
    """Hidden readout errors are retained as CP branches, not postselection."""
    states, weights, metadata = native.preparation(dimension, alpha, thermal)
    rho = np.zeros((2*dimension, 2*dimension), complex)
    rho[:dimension, :dimension] = (states*weights) @ states.conj().T
    spin, spin_metadata = anharmonic._spin_record(dimension)
    records = {}
    for name, detector in operations.items():
        base = tuple(np.kron(np.eye(2), detector[bit]) @ spin[choice]
                     for choice in (0, 1) for bit in (0, 1))
        group = tuple(tuple(math.sqrt(detector_channel[observed, true])*base[2*choice+true]
                            for true in (0, 1))
                      for choice in (0, 1) for observed in (0, 1))
        instrument = runtime.Instrument(kraus=group, visible=True,
            binding={'detector_setting': name, 'classical_bit_channel': detector_channel.tolist(),
                     'all_true_outcomes_summed_as_CP_branches': True})
        program = runtime.Program('complete_readout_bounds', 2*dimension, rho,
            (runtime.Event('spin_and_pointer', instruments=(instrument,)),),
            binding={'preparation': metadata, 'evolution_duration': 0})
        history = runtime.execute(program)
        law = np.zeros(4)
        for branch in history.branches:
            law[branch.record[-1]] += branch.probability
        records[name] = {'joint_probabilities': law.reshape(2, 2).tolist(),
                         'mean': float(law[0]-law[1]+law[2]-law[3]),
                         'mass': float(law.sum()), 'quantum_output_branches': len(history.branches)}
    return {'preparation': metadata, 'spin_pointer': spin_metadata, 'records': records}


def _mixture(first, second, weight):
    result = {}
    for key in first:
        law = ((1-weight)*np.asarray(first[key]['joint_probabilities'])
               +weight*np.asarray(second[key]['joint_probabilities']))
        result[key] = {'joint_probabilities': law.tolist(),
                       'mean': float(law[0, 0]-law[0, 1]+law[1, 0]-law[1, 1]),
                       'mass': float(law.sum()), 'quantum_output_branches': 4}
    return result


def analyze(records, defect):
    channel = bounds.identify_bit_channel(records['null_plus']['mean'], records['null_minus']['mean'])
    corrected = {name: bounds.corrected_mean(row['mean'], channel) for name, row in records.items()}
    z = complex(corrected['h_real'], corrected['h_imag'])
    w = complex(corrected['2h_real'], corrected['2h_imag'])
    center, radius = bounds.moment_disk(z)
    controls = {'null_imag_minus_offset': records['null_imag']['mean']-channel['offset'],
                'h_real_reversal_sum_minus_two_offsets': records['h_real']['mean']+records['h_real_reverse']['mean']-2*channel['offset'],
                'h_imag_reversal_sum_minus_two_offsets': records['h_imag']['mean']+records['h_imag_reverse']['mean']-2*channel['offset']}
    return {'inferred_stationary_bit_channel': channel, 'corrected_means': corrected,
            'first_moment': [z.real, z.imag], 'second_moment': [w.real, w.imag],
            'ideal_second_moment_disk': {'center': [center.real, center.imag], 'radius': radius,
                'real_projection': [max(-1., center.real-radius), min(1., center.real+radius)],
                'imaginary_projection': [max(-1., center.imag-radius), min(1., center.imag+radius)]},
            'moment_witness': bounds.moment_witness(z, w, operator_defect=defect),
            'affine_channel_control_residuals': controls}


def source_bindings():
    paths = set(preparation._source_bindings()) | {
        'model/measurement_bounds/bounds.py', 'model/measurement_bounds/experiment.py',
        'model/measurement_bounds/__main__.py', 'model/measurement_bounds/protocol.json',
        'Lean/EventAlgebra/StateExpectation.lean', 'Lean/EventAlgebra/Basic.lean'}
    return {(ROOT/path).relative_to(ROOT).as_posix(): hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
            for path in sorted(paths)}


def produce():
    protocol = strict_json.load(FOLDER/'protocol.json')
    transfer_protocol = strict_json.load(ROOT/'model/anharmonic_transfer/protocol.json')
    parameters = calibration.transferred_parameters(transfer_protocol)
    alpha, thermal = parameters['alpha_real'], parameters['thermal_occupation']
    h, dimension = protocol['probe_h'], protocol['dimension']
    weight = preparation.coherent_mixture_weight(alpha, thermal, h)
    operations, metadata = compile_detectors(dimension, protocol['settings'], h)
    certificates = bounds.operator_certificates(complete_raw_operators(dimension, operations))
    errors = protocol['declared_detector_channel']
    channel = channel_matrix(errors['flip_zero_to_one'], errors['flip_one_to_zero'])
    components = {name: execute(dimension, a, n, operations, channel)
                  for name, a, n in (('gaussian', alpha, thermal), ('vacuum', 0., 0.),
                                      ('displaced_vacuum', -1.5, 0.))}
    mixture = _mixture(components['vacuum']['records'], components['displaced_vacuum']['records'], weight)
    complete = {'gaussian': components['gaussian']['records'], 'non_gaussian_mixture': mixture}
    analyses = {name: analyze(records, certificates['moment_inequality_defect_upper'])
                for name, records in complete.items()}
    # This deliberately incompatible moment tuple is a software negative
    # control, not a measured trajectory or a positive quantum preparation.
    invalid = bounds.moment_witness(.9+0j, -.8+0j)
    # A setting-dependent offset is a valid classical noise model that violates
    # the stationary-channel assumption. The reversal test reveals it without
    # changing any oscillator state or evolution law.
    control_delta = .04
    changed_channel = channel_matrix(errors['flip_zero_to_one']-control_delta/2,
                                     errors['flip_one_to_zero']+control_delta/2)
    changed = execute(dimension, alpha, thermal,
        {key: operations[key] for key in ('h_imag', 'h_imag_reverse')}, changed_channel)
    changed_residual = (changed['records']['h_imag']['mean']+changed['records']['h_imag_reverse']['mean']
                        -2*analyses['gaussian']['inferred_stationary_bit_channel']['offset'])
    setting_dependent = {
        'injected_offset_shift_at_both_h_imaginary_phases': control_delta,
        'expected_reversal_sum_residual': 2*control_delta,
        'executed_reversal_sum_residual': changed_residual,
        'executed_complete_laws': changed['records'],
        'interpretation': 'a declared detector-channel negative control, not empirical data'}
    source_residuals = []
    spin, _ = anharmonic._spin_record(dimension)
    for component in components.values():
        states, weights, _ = native.preparation(dimension, component['preparation']['alpha_real'],
            component['preparation']['thermal']['mean_thermal_occupation_parameter'])
        state = np.zeros((2*dimension, 2*dimension), complex)
        state[:dimension, :dimension] = (states*weights) @ states.conj().T
        for name, operator in operations.items():
            expected = classical_channel_law(state, spin, operator, channel)
            actual = component['records'][name]['joint_probabilities']
            source_residuals.append(complete_law_residual(actual, expected))
    checks = {
        'every_complete_native_law_normalized': all(abs(row['mass']-1) < 2e-9
            for component in components.values() for row in component['records'].values()),
        'native_classical_channel_matches_full_CP_execution': max(source_residuals) < 2e-12,
        'stationary_bit_channel_recovered_for_both_preparations': all(
            abs(value['inferred_stationary_bit_channel'][key]-errors[key]) < 2e-12
            for value in analyses.values() for key in errors),
        'both_preparations_satisfy_universal_moment_bound': all(
            not value['moment_witness']['violates_with_allowance'] for value in analyses.values()),
        'state_independent_affine_controls_hold': all(abs(v) < 2e-12
            for value in analyses.values() for v in value['affine_channel_control_residuals'].values()),
        'invalid_moment_control_fails': invalid['violates_with_allowance'],
        'executed_setting_dependent_offset_control_fails_stationarity': abs(changed_residual-2*control_delta) < 2e-12,
        'finite_operator_correction_small': certificates['moment_inequality_defect_upper'] < 1e-8}
    if not all(checks.values()):
        raise ArithmeticError(f'measurement bounds failed: {checks}')
    return {'schema': 'oph.measurement_bounds.results.v1', 'all_pass': True, 'checks': checks,
            'protocol': protocol, 'source_bindings': source_bindings(), 'source_controls': compiler.source_summary(),
            'compiled_detector_metadata': metadata, 'finite_operator_certificates': certificates,
            'source_executions': components, 'complete_non_gaussian_mixture_laws': mixture,
            'non_gaussian_weight': weight, 'measurement_analyses': analyses,
            'maximum_full_CP_vs_classical_channel_residual': max(source_residuals),
            'incompatible_moment_negative_control': invalid,
            'setting_dependent_channel_control': setting_dependent,
            'sampling_designs': protocol_sampling_designs(protocol),
            'theorem_attachment': {
                'existing': ['EventAlgebra.expectation_nonneg', 'EventAlgebra.expectation_one'],
                'new_argument': 'positive expectation of X*X gives Gram Cauchy-Schwarz; use X=A*-conj(<A>)I, Y=A-<A>I and bound actual operator defects',
                'new_Lean_theorem_claimed': False},
            'scope': {'new_empirical_observations_loaded': False, 'old_score_changed': False,
                      'Gaussian_preparation_required_for_bound': False,
                      'occupation_limit_required_for_ideal_unitary_bound': False,
                      'stationary_classical_bit_channel_and_phase_controls_supplied': True,
                      'absolute_probe_scale_or_clock_derived': False,
                      'floating_state_propagation_and_physical_service_errors_certified': False}}
