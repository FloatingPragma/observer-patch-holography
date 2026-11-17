"""Execute calibration, record-only recovery and withheld finite-field tests."""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np

from . import native_geometry, recovery, source
from .bounds import coherent_history_cutoff


ROOT = Path(__file__).resolve().parents[2]
CALIBRATION_DURATIONS = (.24, .48)
# These settings are fixed before either occupation-level calibration runs.
WITHHELD = tuple((time, x, y, gamma, eta)
                 for time in (.31, .61)
                 for x, y in ((0, 0), (0, 3), (2, 4), (5, 7))
                 for gamma, eta in ((.3, .4), (.2, .3)))
SOURCE_PATHS = (
    'model/observer_reconstruction/source.py',
    'model/observer_reconstruction/recovery.py',
    'model/observer_reconstruction/native_geometry.py',
    'model/observer_reconstruction/bounds.py',
    'model/observer_reconstruction/experiment.py',
    'model/observer/source_scalar/compiler.py',
    'model/relativistic_scalar/lattice.py',
    'model/scalar_observations/bounds.py',
    'model/scalar_observations/experimental_readout/canonical_source.py',
    'model/observer/native.py',
    'model/common_source/slow_geometry.py',
)


def serializable(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, np.ndarray):
        return serializable(value.tolist())
    if isinstance(value, np.generic):
        return serializable(value.item())
    if isinstance(value, complex):
        return [value.real, value.imag]
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(item) for item in value]
    return value


def _probabilities(chi):
    zero = np.array([(1+chi.real)/2, (1-chi.imag)/2])
    return np.stack((zero, 1-zero), axis=-1)


def level(dimension):
    clock_executions = [source.clock_readout(t) for t in CALIBRATION_DURATIONS]
    clock_rows = [row[0] for row in clock_executions]
    calibration = [source.response_records(dimension, t) for t in CALIBRATION_DURATIONS]
    # The inverse receives only complete-law marginals, interventions and
    # declared alias calibration. It has no backend matrix/coordinate access.
    inferred = recovery.recover(
        calibration[0]['probabilities'][..., 0], calibration[1]['probabilities'][..., 0],
        clock_rows[0][:, 0], clock_rows[1][:, 0], gamma=.3, eta=.4,
        frequency_upper_bound=3.)
    stiffness = inferred['stiffness']
    heldout = []
    for duration, x, y, gamma, eta in WITHHELD:
        clock, clock_execution = source.clock_readout(duration)
        reading = recovery.clock_time(clock[:, 0])
        # Prediction is computed before the withheld field circuit executes.
        predicted = _probabilities(recovery.predict_chi(
            stiffness, reading, x, y, gamma=gamma, eta=eta))
        actual, resources = source.field_readout(dimension, duration, x, y, gamma, eta)
        heldout.append({'kick_site': x, 'read_site': y, 'kick': gamma, 'probe': eta,
            'clock_probabilities': clock, 'clock_reading': reading,
            'clock_execution': clock_execution,
            'predicted_probabilities': predicted, 'executed_probabilities': actual,
            'maximum_probability_error': float(np.max(abs(predicted-actual))),
            'execution': resources})
    # Only AFTER fitting and all withheld predictions do audits open the
    # backend target. These errors are not available to an internal observer.
    target = source.supplied_target()
    fitted_error = float(np.linalg.norm(stiffness-target, 2))
    cutoff = coherent_history_cutoff(dimension, amplitude_squared_upper=F(49, 2000))
    max_error = max(row['maximum_probability_error'] for row in heldout)
    checks = {
        'record_only_inverse_accepted': True,
        'rank_three_cube_recovered': inferred['spatial_dimension'] == 3,
        'independent_generator_audit': fitted_error < 2e-8,
        'withheld_probability_agreement': max_error < 2e-9,
        'all_field_outcomes_normalized': all(
            np.max(abs(row['probabilities'].sum(axis=-1)-1)) < 2e-10 for row in calibration),
    }
    # Coordinates in a degenerate eigenspace have arbitrary orientation.
    # Bind their invariant Gram, not a platform-specific eigenbasis.
    recovered = {key: value for key, value in inferred.items() if key != 'coordinates'}
    return {'occupation_per_normal_mode': dimension,
            'calibration_clock_probabilities': clock_rows,
            'calibration_clock_execution': [row[1] for row in clock_executions],
            'calibration_records': calibration, 'recovered': recovered,
            'withheld_predictions': heldout,
            'maximum_withheld_probability_error': max_error,
            'independent_audit_stiffness_operator_error': fitted_error,
            'ideal_output_occupation_certificate': cutoff,
            'conditional_vacuum_wightman_at_first_withheld_clock':
                recovery.free_wightman(stiffness, heldout[0]['clock_reading']),
            'checks': checks, 'all_pass': all(checks.values())}


def controls(stiffness):
    forward = source.weyl_loop(16)
    reverse = source.weyl_loop(16, reverse=True)
    chi_forward = complex(recovery.characteristic(forward[:, 0]))
    chi_reverse = complex(recovery.characteristic(reverse[:, 0]))
    phase = .3*.4
    loop_residual = max(abs(chi_forward-np.exp(1j*phase)), abs(chi_reverse-np.exp(-1j*phase)))
    # The commuting Q/P control predicts a unit loop, not a statement about
    # all possible classical hidden-state or simulator models.
    commuting = _probabilities(1.+0j)
    zero_q = source.weyl_loop(16, a=0.)
    commuting_margin = float(np.max(abs(forward-commuting)))
    nonlinear = []
    for duration, x, y in ((.31, 0, 0), (.61, 0, 3), (.31, 2, 4)):
        clock, _ = source.clock_readout(duration)
        predicted = _probabilities(recovery.predict_chi(
            stiffness, recovery.clock_time(clock[:, 0]), x, y, gamma=.3, eta=.4))
        actual, _ = source.field_readout(16, duration, x, y, .3, .4, anharmonic=.25)
        nonlinear.append({'clock_probabilities': clock, 'kick_site': x, 'read_site': y,
            'free_prediction': predicted, 'executed_anharmonic_probabilities': actual,
            'maximum_probability_disagreement': float(np.max(abs(predicted-actual)))})
    # A periodic clock alone cannot certify its own unwrapped time label.
    first_clock = source.clock_readout(.24)[0]
    wrapped_clock = source.clock_readout(.24+2*np.pi)[0]
    wrap_error = float(np.max(abs(first_clock-wrapped_clock)))
    return {'weyl_loop': {'forward_probabilities': forward, 'reverse_probabilities': reverse,
                'inferred_forward_characteristic': chi_forward,
                'inferred_reverse_characteristic': chi_reverse,
                'canonical_phase_target': phase, 'target_diagnostic_error': loop_residual,
                'commuting_control_probabilities': commuting,
                'commuting_control_is_analytic_hypothesis_prediction': True,
                'executed_zero_Q_loop_probabilities': zero_q,
                'commuting_control_probability_margin': commuting_margin,
                'ideal_output_occupation_certificate': coherent_history_cutoff(
                    16, amplitude_squared_upper=F(128, 125), modes=1),
                'untrusted_device_quantumness_certified': False},
            'changed_action': {'normal_mode_quartic_coefficient': .25,
                'is_local_site_phi4': False,
                'free_occupation_certificate_applies': False,
                'rows': nonlinear},
            'periodic_clock': {'first_probabilities': first_clock,
                'one_cycle_later_probabilities': wrapped_clock,
                'probability_alias_error': wrap_error,
                'unwrapped_window_is_a_protocol_premise': True},
            'checks': {'noncommuting_phase_reproduced': loop_residual < 2e-9,
                       'executed_zero_Q_loop_is_identity': float(np.max(abs(zero_q-commuting))) < 2e-12,
                       'specified_commuting_control_rejected': commuting_margin > .05,
                       'changed_action_disagrees_with_free_prediction': max(
                           r['maximum_probability_disagreement'] for r in nonlinear) > 1e-3,
                       'clock_cycle_alias_exposed': wrap_error < 2e-12}}


def run():
    levels = [level(d) for d in (8, 16)]
    control = controls(levels[-1]['recovered']['stiffness'])
    native = native_geometry.report()
    result = {
        'schema': 'oph.observer_reconstruction.v1',
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                          for path in SOURCE_PATHS},
        'source_control_bindings': source.c.source_summary()['source_sha256'],
        'protocol': {'gamma': '3/10', 'eta': '2/5', 'frequency_upper_bound': 3.,
            'positive_clock_window': '0 < t < pi; supplied unwrapped phase window',
            'calibration_time_ratio': 2, 'edge_length_convention': 1,
            'probabilities_are_floating_evaluations_of_ensemble_laws': True,
            'finite_shots_or_empirical_data': False,
            'preparation': 'normal-mode vacuum; fresh complete experiment per setting',
            'basis_cutoff': 'normal-mode occupation, distinct from local-site occupation',
            'backend_target_opened_only_after_fit_and_withheld_predictions': True},
        'levels': levels, 'controls': control, 'native_geometry': native,
        'scope': {'field_clocks_and_propagation_share_one_executed_protocol': True,
            'inverse_arithmetic_compiled_to_native_source': False,
            'autonomous_clock_embedded_in_continuing_field_history': False,
            'inferred_observer_time_is_in_supplied_clock_gap_units': True,
            'geometry_recovered_within_tested_cube_hypothesis': True,
            'vacuum_quantum_covariance_requires_declared_CCR_and_state': True,
            'native_carrier_identified_with_field_grid': False,
            'native_repair_selects_field_action': False,
            'inverse_fit_error_or_finite_shot_uncertainty_certified': False,
            'numerical_or_physical_service_error_certified': False,
            'interacting_relativistic_QFT_recovered': False,
            'responsive_gravity_recovered': False,
            'physical_completion_or_new_empirical_prediction': False},
        'all_pass': all(row['all_pass'] for row in levels) and native['all_pass']
                    and all(control['checks'].values()),
    }
    if not result['all_pass']:
        raise ArithmeticError('observer reconstruction failed its declared acceptance checks')
    return serializable(result)
