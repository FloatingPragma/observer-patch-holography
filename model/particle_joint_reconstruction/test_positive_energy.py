from copy import deepcopy
from fractions import Fraction as F
import json

import pytest
import sympy as s

from model.particle_common_contract import contract
from . import positive_energy as positive


def test_actual_selected_tail_uses_reused_finite_theorem_and_same_ids():
    for level, expected in ((0, F(1, 4)), (2, F(1, 16)), (6, F(1, 256))):
        member = contract.member(level)
        result = positive.selected_tail_bound(member)
        assert F(result['original_negative_spectral_weight_upper']) == expected
        assert result['source_id'] == member['ids']['source_id']
        assert result['clock_id'] == member['ids']['clock_id']
        assert result['source_density_stationary']
        assert not result['idle_reference_purification_assumed_stationary']
        assert not result['bare_charged_probe_uses_this_floor']
        assert not result['original_weights_normalized_by_probe_success']


def test_zero_and_large_probe_caps_do_not_normalize_the_original_measure():
    member = contract.member(0)
    assert positive.selected_tail_bound(member, operator_norm_squared='0')['original_negative_spectral_weight_upper'] == '0'
    result = positive.selected_tail_bound(member, threshold='1/128', operator_norm_squared='9')
    assert result['original_negative_spectral_weight_upper'] == '9'
    assert F(2)**(3*6)*F(member['state']['target_ground_excess']) > 1


def test_independent_rotated_matrix_control_retains_all_transitions_and_signs():
    result = positive.rotated_three_level_control('1/16')
    h, rho, a = (result[key] for key in ('H', 'rho', 'A'))
    assert h[0, 1] != 0
    assert rho*h == h*rho
    assert s.trace(rho*h) == s.Rational(1, 16)
    assert len(result['all_nine_transition_rows']) == 9
    nonzero = {omega: weight for omega, weight in result['spectral_measure'].items() if weight}
    assert nonzero == {3: F(31, 32), -1: F(1, 64), -2: F(1, 64)}
    assert sum(nonzero.values()) == 1
    # Independent moment from the literal commutator, not the emitted rows.
    assert s.trace(rho*a.H*(h*a-a*h)) == sum(omega*weight for omega, weight in nonzero.items())
    assert result['idle_purification_energy_variance'] == F(39, 256)
    assert not result['is_selected_spatial_source']


def test_forged_source_or_positive_energy_result_is_rejected():
    member = contract.member(0)
    receipt = positive.bind_member(member)
    positive.verify_bound(member, receipt)
    changed = deepcopy(member)
    changed['state']['source_stationary'] = False
    with pytest.raises(ValueError):
        positive.selected_tail_bound(changed)
    changed = deepcopy(receipt)
    changed['result']['original_negative_spectral_weight_upper'] = '0'
    changed['result_sha256'] = contract.digest(changed['result'])
    with pytest.raises(ValueError):
        positive.verify_bound(member, changed)


def test_energy_escape_has_unit_original_weight_and_unbounded_energy():
    for n in (1, 3, 10):
        result = positive.escaping_energy_control(n)
        assert len(result['excited_energies'])*result['each_original_weight'] == 1
        assert sum(result['excited_energies'])*result['each_original_weight'] == result['mean_excitation_energy']
        assert result['source_preparation_ground_excess'] == 0
        assert not result['continuous_limiting_characteristic_function']


def test_compact_temporal_filter_controls_exact_first_absolute_moment():
    member = contract.member(0)
    result = positive.triangular_filter_bounds(member, half_width='3/2', operator_norm_squared='4')
    assert result['spectral_first_absolute_moment_upper'] == '16/3'
    assert result['correlation_uniform_Lipschitz_constant'] == '16/3'
    assert not result['nonzero_or_nonzero_energy_response_proved']
    assert not result['same_point_time_observation_domain']
    # The scalar domination proving the moment bound holds on both sides
    # of its crossover, including zero. No floating Fourier approximation.
    T = F(3, 2)
    for x in (F(0), F(1, 10), F(4, 3), F(10), F(1000)):
        weighted = F(0) if x == 0 else x*min(F(1), (2/(T*x))**2)
        assert weighted <= 2/T


def test_export_is_json_safe_and_does_not_claim_a_continuum_or_gap():
    result = positive.certificate()
    assert len(json.dumps(result, allow_nan=False)) < 20000
    assert result['continuous_limit_positive_support_is_conditional']
    assert not result['joint_refinement_or_noncollapse_proved']
    assert not result['interacting_relativistic_field_or_mass_gap_proved']


@pytest.mark.parametrize('delta', ['0', '-1', True, 0.25])
def test_invalid_exact_threshold_rejected(delta):
    with pytest.raises(ValueError):
        positive.selected_tail_bound(contract.member(0), threshold=delta)
