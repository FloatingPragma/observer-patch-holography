"""Whole-field calibration followed by complete disturbed-history predictions."""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np

from ..observer_reconstruction.source import clock_readout
from ..observer_reconstruction.experiment import serializable
from ..interacting_reconstruction.experiment import level, portable_bounds, SOURCE_PATHS as BASE_PATHS
from ..interacting_reconstruction.execution import target_audit
from . import bounds, recovery, source


ROOT = Path(__file__).resolve().parents[2]
CALIBRATION_DURATIONS = (F(1, 16384), F(1, 8192))
PRODUCT_NORM_CAP = 96.
LEVELS = ((8,)+(2,)*7, (8,)+(4,)*7)
SOURCE_PATHS = tuple(sorted(set(BASE_PATHS) | {
    'model/coupled_reconstruction/source.py',
    'model/coupled_reconstruction/recovery.py',
    'model/coupled_reconstruction/bounds.py',
    'model/coupled_reconstruction/experiment.py',
    'model/interacting_reconstruction/experiment.py',
}))


def run():
    clocks = [clock_readout(float(t)) for t in CALIBRATION_DURATIONS]
    calibration = [source.calibration_records(t) for t in CALIBRATION_DURATIONS]
    inferred = recovery.recover(*calibration, *(row[0][:, 0] for row in clocks),
        omega=2., product_generator_norm_sum_cap=PRODUCT_NORM_CAP, edge_threshold=.1)
    # Same predeclared coupled measurement protocols as the earlier packet;
    # withheld from this new calibration, not a prospective empirical trial.
    levels = [level(inferred, dimensions) for dimensions in LEVELS]
    audit = target_audit()
    cap = bounds.product_generator_cap_certificate(declared_cap=PRODUCT_NORM_CAP)
    errors = {'mass_squared': abs(inferred['bare_mass_squared']-audit['mass_squared']),
              'quartic_g': abs(inferred['quartic']-audit['quartic_g']),
              'bond_kappa': abs(inferred['bond_kappa']-audit['bond_kappa'])}
    checks = {
        'all_disturbed_history_predictions_pass': all(row['all_pass'] for row in levels),
        'graph_matches_post_prediction_audit': np.array_equal(inferred['adjacency'], audit['adjacency']),
        'bare_coefficients_match_post_prediction_audit': max(errors.values()) < 2e-7,
        'positive_quartic_interaction_identified': inferred['homogeneous_intervals']['quartic']['lower'] > 0,
        'declared_product_norm_cap_covers_actual_schedule': cap['declared_cap_covers_source_schedule'],
    }
    result = {
        'schema': 'oph.coupled_reconstruction.v1',
        'scope': 'Interacting action identification from whole-field source schedules and marginal detector laws. Every onsite and bond term participates; no isolated-generator calibration. Supplied action family, controls and calibrated clock; no interacting continuum or native action selection.',
        'calibration_records': calibration,
        'calibration_clock_probabilities': [row[0] for row in clocks],
        'calibration_clock_execution': [row[1] for row in clocks],
        'protocol': {'calibration_local_dimensions': [4]*8,
            'onsite_transition_preparations': [0, 1],
            'pair_preparation': [0, 0], 'product_steps_per_trial': 1,
            'product_generator_norm_sum_cap': PRODUCT_NORM_CAP,
            'every_onsite_and_bond_term_participates': True,
            'generator_isolation_required': False,
            'same_ordered_schedule_scaled_at_both_durations': True,
            'product_family_is_exact_continuous_semigroup': False,
            'separate_calibration_splitting_error_needed': False,
            'clock_calibration_and_unwrapped_window_supplied': True,
            'probabilities_are_source_executions_not_physical_samples': True,
            'inference_arithmetic_compiled_into_source': False},
        'recovered': {key: value for key, value in inferred.items() if key != 'coordinates'},
        'levels': levels,
        'post_prediction_backend_audit': {**audit, 'absolute_coefficient_errors': errors,
                                          'product_generator_cap_certificate': cap},
        'source_binding_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                  for name in SOURCE_PATHS},
        'floating_clock_and_physical_service_error_certified': False,
        'checks': checks, 'all_pass': all(checks.values()),
    }
    if not result['all_pass']:
        raise ArithmeticError(f'coupled calibration checks failed: {checks}')
    return portable_bounds(serializable(result))
