"""Infer a quartic action from executed records, then predict coupled histories."""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
from itertools import combinations
from pathlib import Path

import numpy as np

from ..observer_reconstruction.source import clock_readout
from ..observer_reconstruction.recovery import clock_time
from ..observer_reconstruction.experiment import serializable
from ..relativistic_scalar import bounds as regulator_bounds, exact_history
from ..relativistic_scalar.experiment import native_records, tv
from ..relativistic_scalar.lattice import ScalarParameters
from . import bounds, execution, recovery, source
from .prediction import RecoveredField


ROOT = Path(__file__).resolve().parents[2]
CALIBRATION_DURATIONS = (F(1, 2048), F(1, 1024))
INTERVAL = F(1, 128)
STEPS = 128
LEVELS = ((8,)+(2,)*7, (8,)+(4,)*7)
OMEGA = 2.
SOURCE_PATHS = (
    'code/source_selection_model/response.json',
    'code/m1_source_realization/receipt.json',
    'model/interacting_reconstruction/source.py',
    'model/interacting_reconstruction/execution.py',
    'model/interacting_reconstruction/recovery.py',
    'model/interacting_reconstruction/prediction.py',
    'model/interacting_reconstruction/bounds.py',
    'model/interacting_reconstruction/experiment.py',
    'model/observer_reconstruction/source.py',
    'model/observer_reconstruction/recovery.py',
    'model/observer_reconstruction/experiment.py',
    'model/observer/source_scalar/compiler.py',
    'model/observer/native.py',
    'model/relativistic_scalar/source.py',
    'model/relativistic_scalar/lattice.py',
    'model/relativistic_scalar/experiment.py',
    'model/relativistic_scalar/bounds.py',
    'model/relativistic_scalar/exact.py',
    'model/relativistic_scalar/exact_history.py',
)


_FLOAT_DEPENDENT_EXACT_BOUND_KEYS = frozenset({
    'finite_time_radius_exact', 'probability_radius_exact', 'radius_exact',
    'exact_hamiltonian_norm_bound', 'exact_complete_history_tv_bound',
})


def portable_bounds(value):
    """Keep numerical bounds without freezing exact expansions of float inputs.

    The inverse decodes floating clock laws, and its propagated radii inherit
    those last bits.  Exact rational arithmetic still certifies each bound
    for those represented inputs, but its long rational string is not a
    portable replay invariant.  The public bound helpers retain that exact
    output and regenerate it on request.  Only these five derived expansions
    are omitted from receipts; fixed rational splitting and occupation proofs,
    source bindings, integer resource counts and all numerical bounds remain.
    """
    if isinstance(value, dict):
        return {key: portable_bounds(item) for key, item in value.items()
                if key not in _FLOAT_DEPENDENT_EXACT_BOUND_KEYS}
    if isinstance(value, list):
        return [portable_bounds(item) for item in value]
    if isinstance(value, tuple):
        return tuple(portable_bounds(item) for item in value)
    return value


def field(inferred, dimensions, *, remove_quartic=False):
    return RecoveredField(dimensions, inferred['onsite_quadratic'],
        np.zeros(8) if remove_quartic else inferred['onsite_quartic'],
        inferred['pair_weights'], omega=OMEGA)


def level(inferred, dimensions):
    clock, clock_execution = clock_readout(float(INTERVAL))
    reading = clock_time(clock[:, 0])
    predicted_field = field(inferred, dimensions)
    quadratic_control = field(inferred, dimensions, remove_quartic=True)
    predictions, controls, actual, execution_metadata = [], [], [], []
    # Complete predictions are made before any corresponding withheld field
    # circuit. Neither source coordinates nor source matrices enter them.
    for sign in (-1, 1):
        predictions.append(predicted_field.history(reading, sign=sign))
        controls.append(quadratic_control.history(reading, sign=sign))
    for sign in (-1, 1):
        rows, metadata = execution.coupled_history(dimensions, INTERVAL, sign=sign, steps=STEPS)
        actual.append(rows)
        execution_metadata.append(metadata)
    parameter = bounds.parameter_uncertainty_bound(dimensions,
        inferred['onsite_quadratic_radius'], inferred['onsite_quartic_radius'],
        {(i, j): inferred['pair_weight_radii'][i, j] for i, j in combinations(range(8), 2)},
        2*INTERVAL, omega=OMEGA)
    splitting = min(F(1), 2*regulator_bounds.product_error(dimensions, INTERVAL, STEPS))
    total = min(1., parameter['complete_history_tv_upper']+float(splitting))
    discrepancies = [tv(a, b) for a, b in zip(predictions, actual)]
    control_errors = [tv(a, b) for a, b in zip(controls, actual)]
    # The supplied homogeneous cube has an automorphism taking the actual
    # opaque read site to the certificate's canonical origin. It preserves
    # the unique eight-level read code, seven equal neighbor cutoffs, blank,
    # onsite/bond action, pointer phase, and first-bit feedback.
    occupation = [exact_history.comparison_certificate(dimensions, time=INTERVAL,
                    degree=2, instrument_sign=sign) for sign in (-1, 1)]
    cutoff = max(row['complete_history_tv_upper'] for row in occupation)
    native_prediction = native_records(predictions)
    native_actual = native_records(actual)
    native_tv = sum(abs(a['probability']-b['probability'])
                    for a, b in zip(native_prediction, native_actual))/2
    checks = {
        'all_history_outcomes_normalized': all(abs(sum(r['probability'] for r in h)-1) < 2e-9
                                               for h in predictions+actual+controls),
        'coupled_prediction_within_ideal_bound_plus_numerical_guard': max(discrepancies) < total+2e-8,
        'coupled_prediction_numerically_accurate': max(discrepancies) < 2e-8,
        'quartic_deletion_disagrees_beyond_interacting_prediction_budget': min(control_errors) > total+2e-8,
        'native_joint_law_conserves_probability': abs(sum(r['probability'] for r in native_actual)-1) < 2e-9,
        'native_attachment_does_not_increase_tv': native_tv <= max(discrepancies)+2e-12,
    }
    return {
        'local_dimensions': dimensions, 'tensor_dimension': predicted_field.dimension,
        'read_site_opaque_label': 0, 'model_interval_exact': str(INTERVAL),
        'clock_probabilities': clock, 'clock_reading': reading, 'clock_execution': clock_execution,
        'predicted_histories_by_native_phase': predictions,
        'executed_histories_by_native_phase': actual,
        'quadratic_deletion_predictions_by_native_phase': controls,
        'prediction_history_tv_by_phase': discrepancies,
        'quadratic_deletion_history_tv_by_phase': control_errors,
        'execution_by_phase': execution_metadata,
        'parameter_error_certificate': parameter,
        'ideal_splitting_history_tv_upper_exact': str(splitting),
        'ideal_inferred_to_source_history_tv_upper': total,
        'occupation_comparisons_by_phase': occupation,
        'ideal_source_to_untruncated_tv_upper': float(min(F(1), cutoff+splitting)),
        'ideal_inferred_continuous_to_untruncated_tv_upper': min(1., float(cutoff)+parameter['complete_history_tv_upper']),
        'native_predicted_history': native_prediction, 'native_executed_history': native_actual,
        'native_history_count': len(native_actual), 'native_prediction_tv': native_tv,
        'checks': checks, 'all_pass': all(checks.values()),
    }


def run():
    clocks = [clock_readout(float(t)) for t in CALIBRATION_DURATIONS]
    calibration = [source.calibration_records(t) for t in CALIBRATION_DURATIONS]
    kwargs = dict(omega=OMEGA, norm_bound_onsite=16., norm_bound_bond=2., edge_threshold=.1)
    inferred = recovery.recover(*calibration, *(row[0][:, 0] for row in clocks), **kwargs)
    levels = [level(inferred, dimensions) for dimensions in LEVELS]
    # Extra measurement disturbance test: same recovered action, no feedback.
    no_feedback_prediction = field(inferred, LEVELS[0]).history(
        levels[0]['clock_reading'], sign=-1, feedback_phase=0.)
    no_feedback_actual, no_feedback_execution = execution.coupled_history(
        LEVELS[0], INTERVAL, sign=-1, steps=STEPS, feedback_phase=0.)
    feedback_change = tv(no_feedback_actual, levels[0]['executed_histories_by_native_phase'][0])
    # This is a new source execution with lambda zero, not a relabelled law.
    free_elements = []
    for transition in (0, 1):
        records = [source.onsite_probe(0, transition, t, ScalarParameters(coupling=0))[0]
                   for t in CALIBRATION_DURATIONS]
        free_elements.append(recovery.coefficient_from_records(np.stack(records),
            *inferred['clock_readings'], norm_bound=16.))
    free_local = recovery.onsite_parameters(
        *(r['matrix_element'] for r in free_elements), *(r['radius'] for r in free_elements), omega=OMEGA)
    sensitivity = recovery.recover(*calibration, *(row[0][:, 0] for row in clocks),
                                  probability_error=1e-7, **kwargs)
    # Audit opens the backend only after all predictions/executions. These
    # numbers never select or tune the inferred action or withheld settings.
    audit = execution.target_audit()
    audit_errors = {'mass_squared': abs(inferred['bare_mass_squared']-audit['mass_squared']),
                    'quartic_g': abs(inferred['quartic']-audit['quartic_g']),
                    'bond_kappa': abs(inferred['bond_kappa']-audit['bond_kappa'])}
    checks = {
        'all_heldout_levels_pass': all(row['all_pass'] for row in levels),
        'cube_graph_matches_independent_audit': np.array_equal(inferred['adjacency'], audit['adjacency']),
        'recovered_bare_coefficients_match_audit': max(audit_errors.values()) < 2e-8,
        'supplied_norm_caps_cover_actual_calibration_generators': audit['onsite_norm'] <= 16 and audit['pair_norm'] <= 2,
        'positive_quartic_interval_excludes_quadratic_action': inferred['homogeneous_intervals']['quartic']['lower'] > 0,
        'executed_quartic_deletion_local_interval_contains_zero': abs(free_local['quartic']) <= free_local['quartic_radius'],
        'feedback_deletion_changes_history': feedback_change > 1e-3,
        'feedback_deletion_is_also_predicted': tv(no_feedback_prediction, no_feedback_actual) < 2e-8,
    }
    result = {
        'schema': 'oph.interacting_reconstruction.v1',
        'scope': 'Record-only identification of a supplied finite local quartic scalar action, with coupled disturbed-history predictions. No interacting continuum, physical clock, or native action selection.',
        'calibration_records': calibration, 'calibration_clock_probabilities': [row[0] for row in clocks],
        'calibration_clock_execution': [row[1] for row in clocks],
        'protocol': {'basis_frequency': OMEGA, 'onsite_occupation': 4, 'bond_occupation_per_site': 2,
                     'norm_bound_onsite': 16., 'norm_bound_bond': 2., 'edge_separation_threshold': .1,
                     'record_laws_are_exact_probability_executions_not_samples': True,
                     'isolation_controls_and_clock_calibration_supplied': True,
                     'inference_arithmetic_compiled_into_source': False,
                     'time_parameter': 'separately calibrated phase-clock model-coordinate reading'},
        'recovered': {key: value for key, value in inferred.items() if key != 'coordinates'},
        'levels': levels,
        'controls': {'executed_free_onsite_parameters': free_local,
                     'no_feedback_prediction': no_feedback_prediction,
                     'no_feedback_execution': no_feedback_actual,
                     'no_feedback_resources': no_feedback_execution,
                     'feedback_deletion_history_tv': feedback_change},
        'finite_sampling_sensitivity': {
            'coefficient_intervals_if_each_calibration_law_has_error_at_most_1e_minus_7': sensitivity['homogeneous_intervals'],
            'clock_laws_treated_as_exact': True,
            'independent_Bernoulli_budget': bounds.simultaneous_shot_count(F(1, 10**7), 176),
            'finite_sampling_was_executed': False},
        'independent_backend_audit': {**audit, 'absolute_coefficient_errors': audit_errors},
        'source_binding_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_PATHS},
        'floating_or_physical_service_error_certified': False,
        'checks': checks, 'all_pass': all(checks.values()),
    }
    if not result['all_pass']:
        raise ArithmeticError(f'interacting reconstruction checks failed: {checks}')
    return serializable(portable_bounds(result))
