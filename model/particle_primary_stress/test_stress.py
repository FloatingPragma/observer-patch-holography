"""Independent source attachment, local allocation and current controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import pytest
import sympy as s

from model.particle_common_contract import contract
from . import stress


@pytest.fixture(scope='module')
def first(): return contract.member(0)


def test_every_open_source_term_is_allocated_once_with_actual_owner_supports():
    for side in (2, 3, 4):
        result = stress.incidence_control(side)
        assert result['allocated_terms'] == 4*side**3+3*side**2*(side-1)+3*side*(side-1)**2
        assert result['each_source_term_allocated_once']
        assert result['maximum_overlap_graph_distance'] <= 4
    bulk = stress.incidence_control(3)
    assert bulk['maximum_plaquettes_per_owner_cell'] == 9
    assert bulk['maximum_cell_energy_norm_upper'] == F(15315, 14)
    assert bulk['maximum_overlap_graph_distance'] == 4


def test_primary_cell_cap_uses_actual_electric_scaling_and_full_fock_norms():
    member = contract.member(1)
    actual = stress.cell_energy(member, (1, 1, 1))
    a = F(member['source']['spacing']); ell = member['source']['per_link_word_degree']
    assert actual['norm_upper'] == (180+F(1827, 2)*ell**2+F(3, 7))/a
    assert sum(row['kind'] == 'bond' for row in actual['terms']) == 6
    assert sum(row['kind'] == 'plaquette' for row in actual['terms']) == 9
    assert all(sum(abs(x-y) for x, y in zip((1, 1, 1), cell)) <= 2 for cell in actual['owner_support'])


def test_origin_and_ghost_link_terms_are_retained(first):
    origin = stress.cell_energy(first)
    assert len(origin['owner_support']) == 4
    assert sum(row['kind'] == 'electric' for row in origin['terms']) == 3
    corner = stress.cell_energy(first, (1, 1, 1))
    assert sum(row['kind'] == 'electric' for row in corner['terms']) == 3
    assert sum(row['kind'] == 'bond' for row in corner['terms']) == 3
    assert sum(row['kind'] == 'plaquette' for row in corner['terms']) == 0


def test_whole_source_has_no_external_flux_but_regional_pairs_are_complete(first):
    whole = stress.regional_bounds(first, tuple(product(range(2), repeat=3)))
    assert whole['crossing_pairs'] == []
    assert whole['outward_flux_norm_upper'] == 0
    region = stress.regional_bounds(first)
    assert len(region['crossing_pairs']) > 0
    assert len({tuple(row['to']) for row in region['crossing_pairs']}) == len(region['crossing_pairs'])
    assert all(tuple(row['from']) == (0, 0, 0) for row in region['crossing_pairs'])


def test_readout_bound_retains_state_error_and_quadrature_error(first):
    budget = stress.complete_record_budget(first, duration='1/100', half_trace_error='1/1000', quadrature_slices=7)
    k, m, b = (budget[key] for key in ('regional_energy_norm_upper','outward_flux_norm_upper','source_H_norm_upper'))
    assert budget['measured_continuity_residual_upper'] == 4*k/1000+2*m/100000+b*m/70000
    zero = stress.complete_record_budget(first, duration='0', half_trace_error='0')
    assert zero['measured_continuity_residual_upper'] == 0
    fine = stress.complete_record_budget(first, duration='1/100', half_trace_error='1/1000', quadrature_slices=14)
    assert fine['left_rule_quadrature_error_upper']*2 == budget['left_rule_quadrature_error_upper']


def test_original_source_is_regenerated_before_admitting_current_receipt(first):
    receipt = stress.common_contract_bound(first)
    stress.verify_bound(first, receipt)
    changed = deepcopy(first)
    changed['reference']['idle_H'] = 1
    with pytest.raises(ValueError): stress.common_contract_payload(changed)
    changed = deepcopy(receipt)
    changed['result']['stationary_selected_source_regional_net_flux_mean'] = '1'
    changed['result_sha256'] = contract.digest(changed['result'])
    with pytest.raises(ValueError): stress.verify_bound(first, changed)


def test_independent_control_current_has_exact_nonzero_physical_response():
    value = stress.quantum_control()
    h0, h1 = value['h']; total = value['H']; current = value['J']
    assert current == s.I*(h0*h1-h1*h0)
    assert current == current.H
    assert s.I*(total*h0-h0*total)+current == s.zeros(48)
    assert s.I*(total*h1-h1*total)-current == s.zeros(48)
    assert all(op*g == g*op for op in (h0, h1, current) for g in value['Gauss'])
    assert value['physical_dimension'] == 6
    assert value['intervened_outward_flux'] == -s.Rational(3600, 124697)
    assert value['intervened_cell0_energy_rate'] == s.Rational(3600, 124697)


def test_complete_quantum_control_records_retain_both_outcomes_and_reference_domain():
    value = stress.quantum_control()
    assert sum(value['effects'], s.zeros(48)) == s.eye(48)
    assert value['record_probabilities'] == (s.Rational(121697, 249394), s.Rational(127697, 249394))
    assert value['record_cap'] == s.Rational(6, 5)
    assert value['selected_control_mean_current'] == 0
    assert value['idle_purifier_energy_variance'] > 0
    assert value['rho']*value['H'] == value['H']*value['rho']
    assert value['intervened_rho']*value['H'] != value['H']*value['intervened_rho']
    assert not value['same_as_selected_spatial_primary']


def test_supported_entropy_floor_keeps_coin_and_marginal_distinct(first):
    result = stress.regional_first_law_bound(first)
    assert len(result['region_owner_cells']) == 4
    exponent = result['supported_regional_density_floor']['negative_exponent']
    original = result['inherited_selected_support_consumer']['physical_density_floor_factored']['negative_exponent']
    assert exponent == original+4*87
    assert result['regional_modular_norm_upper'] == exponent
    assert result['small_coin_probability_example']['negative_exponent'] == exponent+4
    assert result['small_coin_relative_entropy_remainder_upper']['negative_exponent'] == exponent+8
    assert result['complete_coin_and_quantum_reference_retained']
    assert not result['marginal_first_law_is_full_observer_entropy']
    assert not result['ambient_region_full_rank_claimed']
    assert not result['regional_modular_operator_is_beta_times_cell_energy']


@pytest.mark.parametrize('region', [(), ((0, 0, 0),(0, 0, 0)), ((2, 0, 0),), ((True, 0, 0),)])
def test_invalid_or_duplicate_regions_are_rejected(first, region):
    with pytest.raises(ValueError): stress.regional_bounds(first, region)
