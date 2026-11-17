"""Fresh native-interface certification and exact count-analysis controls."""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np

from tools import strict_json
from ..measurement_bounds import bounds as original_bounds
from ..measurement_bounds import experiment as source
from .bounds import FIRST, SETTINGS, InterfaceBounds, analyze_counts, json_safe, sampling_radius


FOLDER = Path(__file__).resolve().parent
ROOT = FOLDER.parents[1]


def certify_operators(operators, guarantees):
    """Verify fixed rational outward guarantees on this fresh source execution.

    The older interval kernel is reused only after concrete ndarray and menu
    validation. Returned Kraus arrays are not modified by this packet. Their
    upstream detector producer uses its existing guarded polar representation;
    the raw pulse-product-to-representation error is not certified here.
    """
    required = set(SETTINGS) | {'null_imag', 'h_real_reverse', 'h_imag_reverse'}
    if type(operators) is not dict or set(operators) != required:
        raise ValueError('the complete source operator menu is required')
    shape = None
    for groups in operators.values():
        if type(groups) is not tuple or len(groups) != 2:
            raise ValueError('both detector outcomes required')
        for group in groups:
            if type(group) is not tuple or len(group) != 2:
                raise ValueError('both spin outcomes required')
            for matrix in group:
                if type(matrix) is not np.ndarray:
                    raise ValueError('original base ndarray source matrices required')
                if (matrix.ndim != 2 or matrix.shape[0] < 1 or matrix.shape[0] != matrix.shape[1]
                        or matrix.dtype.kind not in 'fc' or matrix.dtype.itemsize > (16 if matrix.dtype.kind == 'c' else 8)
                        or not np.isfinite(matrix).all()):
                    raise ValueError('finite square binary source matrices required')
                if shape is None:
                    shape = matrix.shape
                if matrix.shape != shape:
                    raise ValueError('all settings must act on the same carrier')
    if type(guarantees) is not dict or set(guarantees) != {
            'effect_sum_norm_upper', 'null_observable_norm_upper', 'moment_defect_upper'}:
        raise ValueError('complete exact interface guarantees required')
    mass, null, moment = (guarantees[key] for key in
        ('effect_sum_norm_upper', 'null_observable_norm_upper', 'moment_defect_upper'))
    interface = InterfaceBounds((mass,)*6, (null,)*2, moment)
    measured = original_bounds.operator_certificates(operators)
    if (any(F(value) > mass for value in measured['complete_instrument_effect_sum_defects_upper'].values())
            or F(measured['moment_inequality_defect_upper']) > moment
            or any(F(measured['state_independent_control_defects_upper'][key]) > null
                   for key in ('null_plus_minus_identity', 'null_minus_plus_identity'))):
        raise ArithmeticError('fresh source matrices exceeded a declared interface guarantee')
    return interface


def fresh_interface():
    protocol = strict_json.load(FOLDER/'protocol.json')
    original = strict_json.load(source.FOLDER/'protocol.json')
    if protocol['measured_settings'] != list(SETTINGS):
        raise ValueError('the fixed six-setting confidence family changed')
    operations, metadata = source.compile_detectors(protocol['dimension'], original['settings'], protocol['probe_h'])
    complete = source.complete_raw_operators(protocol['dimension'], operations)
    guarantees = {key: F(value) for key, value in protocol['finite_interface_guarantees'].items()}
    return protocol, certify_operators(complete, guarantees), operations, metadata


def count_fixture(means, trials):
    """Deterministic rounded mean fixtures, never observations or random draws.

    Spin-zero labels are assigned explicitly for this software table; both
    spin-one counts remain present as zeros. Inference itself sums all rows.
    """
    if type(trials) is not int or not 1 <= trials <= 10**12 or type(means) is not dict or set(means) != set(SETTINGS):
        raise ValueError('fixed integer trials and all six fixture means required')
    result = {}
    for name in SETTINGS:
        if isinstance(means[name], (bool, np.bool_)):
            raise ValueError('Boolean fixture means are inadmissible')
        mean = F(means[name])
        if not -1 <= mean <= 1:
            raise ValueError('binary fixture mean outside [-1,1]')
        zeros = round(trials*(1+mean)/2)
        result[name] = [zeros, trials-zeros, 0, 0]
    return result


def source_bindings():
    paths = set(source.source_bindings()) | {
        'model/calibrated_measurement/__init__.py', 'model/calibrated_measurement/__main__.py',
        'model/calibrated_measurement/bounds.py', 'model/calibrated_measurement/experiment.py',
        'model/calibrated_measurement/protocol.json', 'Lean/Dynamics/ChoiCPTP.lean',
        'model/common_source/output_certificate.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def produce():
    protocol, interface, operations, metadata = fresh_interface()
    transfer = strict_json.load(ROOT/'model/anharmonic_transfer/protocol.json')
    parameters = source.calibration.transferred_parameters(transfer)
    errors = protocol['detector_channel_example']
    execution = source.execute(protocol['dimension'], parameters['alpha_real'], parameters['thermal_occupation'],
        {key: operations[key] for key in SETTINGS},
        source.channel_matrix(errors['flip_zero_to_one'], errors['flip_one_to_zero']))
    means = {name: row['mean'] for name, row in execution['records'].items()}
    zero_service = dict.fromkeys(SETTINGS, F(0))
    alpha = F(protocol['failure_probability'])
    examples = []
    for trials in protocol['illustrative_trial_budgets']:
        counts = count_fixture(means, trials)
        examples.append({'kind': 'deterministic software fixture from rounded source expectations; not sampled data',
            'counts': counts, 'analysis': analyze_counts(counts, interface,
                service_mean_errors=zero_service, failure_probability=alpha)})
    invalid = count_fixture({'null_plus': F(9, 10), 'null_minus': F(-4, 5),
        'h_real': F(163, 200), 'h_imag': F(1, 20),
        '2h_real': F(-63, 100), '2h_imag': F(1, 20)}, 20000)
    negative = analyze_counts(invalid, interface, service_mean_errors=zero_service, failure_probability=alpha)
    uncertainty = analyze_counts(invalid, interface,
        service_mean_errors=dict.fromkeys(SETTINGS, F(1, 4)), failure_probability=alpha)
    missing = analyze_counts(invalid, interface, service_mean_errors=dict.fromkeys(SETTINGS), failure_probability=alpha)
    prefix = analyze_counts({key: examples[-1]['counts'][key] for key in FIRST}, interface,
        service_mean_errors=zero_service, failure_probability=alpha)
    checks = {
        'fresh_native_interface_guarantees_pass': True,
        'all_source_fixtures_not_rejected': all(not row['analysis']['rejects_joint_premises'] for row in examples),
        'incompatible_moment_fixture_rejected': negative['rejects_joint_premises'],
        'large_service_uncertainty_prevents_rejection': not uncertainty['rejects_joint_premises'],
        'missing_service_bounds_leave_obligation_open': missing['status'] == 'unbounded_services',
        'prospective_prediction_uses_only_first_four_settings': prefix['status'] == 'prediction_only'
            and prefix['second_observed_mean_predictions'] == examples[-1]['analysis']['second_observed_mean_predictions'],
    }
    if not all(checks.values()):
        raise ArithmeticError(f'calibrated measurement checks failed: {checks}')
    return json_safe({'schema': 'oph.calibrated_measurement.results.v1', 'all_pass': True,
        'checks': checks, 'protocol': protocol, 'source_bindings': source_bindings(),
        'source_detector_metadata': metadata, 'finite_interface': {
            'mass_bounds': interface.mass, 'null_bounds': interface.null, 'moment_bound': interface.moment,
            'complete_instrument_observable_corrections': {name: interface.completion(name) for name in SETTINGS},
            'returned_kraus_arrays_modified_by_this_packet': False,
            'inherited_guarded_polar_representation': True,
            'raw_pulse_product_to_returned_representation_error_certified': False},
        'source_expectation_execution_diagnostic': execution,
        'illustrative_count_analyses': examples, 'prospective_four_setting_prediction': prefix,
        'incompatible_count_fixture': {'counts': invalid, 'analysis': negative},
        'service_uncertainty_control': uncertainty, 'missing_service_control': missing,
        'sampling_designs': [{'trials_per_setting': n, 'simultaneous_raw_mean_radius_upper': sampling_radius(n, alpha)}
                             for n in protocol['illustrative_trial_budgets']],
        'scope': {'fresh_returned_binary_Kraus_operator_errors_certified': True,
            'count_to_decision_scalar_arithmetic_exact_or_outward': True,
            'Gaussian_preparation_required': False, 'unknown_fixed_affine_detector_calibrated': True,
            'physical_service_mean_errors': None, 'source_expectation_propagation_certified': False,
            'new_empirical_or_random_sample_outcomes': False, 'old_anharmonic_score_changed': False,
            'stationary_detector_and_stable_probe_preparation_assumed': True,
            'physical_probe_scale_clock_or_action_derived': False,
            'distinctive_OPH_empirical_prediction': False},
        'theorem_attachment': ['finite positive state expectation and Gram Cauchy-Schwarz',
            'complete Kraus instruments and spectral calculus for positive T',
            'Hoeffding 1963 Theorem 2, applied to independent bounded trials, plus a six-setting union bound']})
