"""Retained native fluctuation-action and complete current-record predictions."""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
from pathlib import Path

from . import source
from ..primitive_observer.transfer import native_counts


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    'model/native_fluctuations/source.py', 'model/native_fluctuations/experiment.py',
    'model/native_fluctuations/__main__.py', 'model/primitive_observer/transfer.py',
    'model/primitive_observer/runtime.py', 'model/common_source/native_fields.py',
    'model/common_source/open_protocol.py', 'model/common_source/open_environment.py',
    'model/common_source/protocol.py',
    'model/common_source/slow_geometry.py', 'model/observer/native.py',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py',
    'model/requirements.txt', 'evidence/native_geometric_source_20260925/spec.json',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
    'Lean/SourceLaw/ReflectionPositivity.lean',
)


def run():
    kernel = native_counts()
    a, c, d = source.source_moments(kernel)
    law = source.word_law(kernel)
    score = source.score_distribution(law)
    variance, third, fourth_cumulant = source.score_moments(score)
    moments = [sum((word[0]*word[lag]*sum(row) for word, row in law.items()), F(0))
               for lag in range(3)]
    r = tuple(tuple(F(i == j)-a[i][j] for j in range(12)) for i in range(12))
    cr = source.matmul(c, r)
    crr = source.matmul(cr, r)
    crra = source.matmul(crr, a)
    expected = [2*cr[0][0], -crr[0][0], -crra[0][0]]
    if moments != expected or moments[0] != F(5, 198) or moments[1] != -F(1, 792):
        raise ArithmeticError('recorded signed currents fail the native history covariance')
    if variance != sum(source.WEIGHTS[i]*source.WEIGHTS[j]*moments[abs(i-j)]
                       for i in range(3) for j in range(3)):
        raise ArithmeticError('the executed contrast variance fails the derived Gaussian history action')
    if fourth_cumulant == 0:
        raise ArithmeticError('the finite native non-Gaussian control disappeared')
    numerical = source.execute_current_check(kernel)
    return {
        'schema': 'oph.native_fluctuations.v1', 'all_pass': True,
        'sources': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'source_attachment': {
            'native_configurations': 66, 'primitive_branches': 3960,
            'conditional_mean_checked_on_every_configuration': True,
            'lag_covariances_checked_by_finite_state_propagation': [0, 1, 2, 3, 4],
            'mean_transition_A': [list(map(str, row)) for row in a],
            'stationary_covariance_C': [list(map(str, row)) for row in c],
            'innovation_covariance_D': [list(map(str, row)) for row in d],
            'C_identity': '(5/33)*(I-J/12)',
            'D_identity': 'C-A*C*A^T, positive on the zero-total-load sector',
        },
        'Gaussian_history_action': source.gaussian_action_check(a, c, d),
        'executed_signed_current_history': numerical,
        'accessible_current_predictions': {
            'observable': 'signed native load change at port zero, computed from selected seam and current',
            'one_attempt_variance': str(moments[0]),
            'successive_attempt_covariance': str(moments[1]),
            'successive_attempt_correlation': str(moments[1]/moments[0]),
            'two_interval_covariance': str(moments[2]),
            'no_initial_configuration_measurement_required': True,
        },
        'finite_replica_prediction': {
            'score_weights_on_three_signed_records': list(source.WEIGHTS),
            'complete_one_carrier_score_law': [{'score': z, 'probability': str(p)} for z, p in score.items()],
            'one_carrier_variance': str(variance),
            'one_carrier_absolute_third_moment': str(third),
            'one_carrier_fourth_cumulant': str(fourth_cumulant),
            'scaling': 'M^(-1/2) times the sum of M independent complete native histories',
            'comparisons': [source.characteristic_comparison(score, m) for m in (1, 4, 16, 64, 256)],
            'product_law_is_factorized_not_a_dense_tensor_simulation': True,
            'all_product_outcomes_included_without_tail_cut': True,
        },
        'scope': {
            'independent_stationary_product_histories_supplied': True,
            'uniform_seam_scheduler_and_fair_coin_supplied': True,
            'Gaussian_limit_and_action_derived_from_native_histories': True,
            'finite_dimensional_characteristic_bound_is_analytic': True,
            'new_Lean_CLT_claimed': False,
            'single_native_carrier_is_Gaussian': False,
            'physical_clock_or_spatial_continuum_derived': False,
            'interacting_quantum_field_or_real_time_evolution_derived': False,
            'source_selects_complete_quantum_instrument': False,
            'new_empirical_dataset_tested': False,
            'native_pulse_words_emitted': False,
            'numerical_diagnostics_are_outward_roundoff_certificates': False,
        },
    }
