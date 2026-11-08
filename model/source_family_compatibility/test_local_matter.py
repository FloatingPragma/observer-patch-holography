"""Adversarial exact checks for complete local interacting repair histories."""
from fractions import Fraction as F

import pytest
import sympy as s

from model.full_fock_field.bounds import local_boundary_bound
from model.source_family_compatibility import local_matter as m


@pytest.mark.parametrize('law', ['C', 'E'])
def test_disjoint_crossing_seam_fixes_entire_local_car_basis(law):
    identity = s.eye(64)
    a, b = m.annihilator(3, 2), m.annihilator(3, 3)
    basis = [x*y for x in (identity, a, a.H, a.H*a)
             for y in (identity, b, b.H, b.H*b)]
    # The 16 local elements are independent, so checking them is exhaustive.
    assert s.Matrix.hstack(*(s.Matrix(x).reshape(4096, 1) for x in basis)).rank() == 16
    assert all(m.repair_dual(x, sites=3, seam=(0, 2), law=law) == x for x in basis)
    assert m.repair_dual(a, sites=3, seam=(0, 2), law=law, omit_signs=True) != a


def test_complete_archive_countercontrol_is_nonzero_without_local_disturbance():
    result = m.finite_locality_controls()
    assert result['original_remote_collision_probability'] == '1/2'
    assert result['coin_and_raw_slots_per_law'] == 32
    assert result['archive_countercontrol_masks'] == [21, 14]
    assert not result['global_archive_equality_claimed']


def test_centered_same_ground_interaction_derivative_on_whole_fock_space():
    result = m.centered_interaction_control()
    assert result['ground_H_energy'] == '-1/2'
    assert result['full_ground_rank'] == 1
    assert result['one_site_opposite_occupation'] == '1/2'
    assert result['exact_four_point_variance_quadratic_coefficient'] == '4'
    assert not result['finite_dimer_is_infinite_ground_identification']


def test_original_Weyl_action_and_bare_repair_are_directly_local():
    result=m.original_weyl_bare_repair_controls()
    assert [row['axis']for row in result['source_cases']]==[0,1,2]
    assert all(row['charged_second_moment_sum_rule']=='41/8'for row in result['source_cases'])
    assert result['bare_disjoint_physical_seam']==[1,2]
    assert result['local_CAR_basis_elements']==16
    assert not result['repair_Kraus_maps_conjugated']
    assert not result['hidden_alignment_or_chart_covariance_required']


@pytest.mark.parametrize('radius,time,coupling', [(20,F(1,8),1), (80,F(1,2),1), (160,1,1)])
def test_dyadic_boundary_encloses_inherited_rational_lr(radius,time,coupling):
    inherited = local_boundary_bound(radius, 0, 1, time, coupling)
    assert m.boundary_upper(radius, 0, 1, time, coupling) >= inherited


@pytest.mark.parametrize('eta,attempts,interventions', [(F(1,10),8,2),(F(1,10000),100,6),(F(1,2),0,0)])
def test_rational_planner_controls_requested_complete_local_error(eta,attempts,interventions):
    result = m.plan_local_history(1, 4, F(1,8), 1, interventions, attempts, eta)
    assert result['accessible_output_halftrace_upper'] <= eta
    assert result['outer_cube_side'] >= 2*result['radius']+1
    assert result['outer_cube_edges'] == 3*(result['outer_cube_side']-1)*result['outer_cube_side']**2
    assert result['outer_cube_sites'] == result['outer_cube_side']**3
    assert not result['global_raw_archive_compared']
    assert not result['discarded_raw_feedback_allowed']
    assert result['initial_reference_entanglement_allowed']


def test_more_interventions_and_repairs_are_charged_to_bound():
    common = (60,0,1,F(1,8),1)
    base = m.local_history_bound(*common,0,1,10**12)
    doubled = m.local_history_bound(*common,0,2,10**12)
    controls = m.local_history_bound(*common,3,1,10**12)
    assert doubled['accessible_output_halftrace_upper'] == 2*base['accessible_output_halftrace_upper']
    assert controls['suffix_operator_error'] == 7*base['suffix_operator_error']


def test_zero_time_has_only_near_repair_cost():
    result = m.local_history_bound(2,0,1,0,1,4,3,10**6)
    assert result['suffix_operator_error'] == 0
    assert result['accessible_output_halftrace_upper'] == F(2250,10**6)


@pytest.mark.parametrize('args', [(0,1,F(1,8),1,1,1,0), (0,2,1,1,0,1,F(1,2)),
                                  (0,1,1.,1,0,1,F(1,2)), (0,1,1,1,-1,1,F(1,2)),
                                  (0,1,1,1,0,True,F(1,2)), (0,1,1,1,0,1,2)])
def test_invalid_planner_domains_rejected(args):
    with pytest.raises(ValueError):
        m.plan_local_history(*args)


def test_interaction_interval_and_quantitative_connected_lower_bound():
    c,u = F(3,2),F(4)
    maximum = u/(8*(9*c*c+12*c*u+u*u))
    result = m.interaction_witness(c,u,maximum)
    assert result['four_point_variance_lower'] == 4*result['one_connected_four_point_absolute_lower'] > 0
    assert result['variance_quadratic_coefficient'] == 4
    with pytest.raises(ValueError):
        m.interaction_witness(c,u,maximum+F(1,10**9))
    with pytest.raises(ValueError):
        m.interaction_witness(c,0,F(1,1000))


def test_domain_controls_do_not_promote_local_construction_to_particle_assignment():
    result = m.certificate()
    assert result['odd_matter_and_even_spin_share_full_CAR_action']
    assert result['original_Weyl_bare_repair_local_history_instantiated']
    assert not result['bare_repairs_identified_under_site_dependent_frame']
    assert result['same_centered_thermodynamic_ground_subsequence']
    assert not result['effective_ground_subsequence_preparation_modulus']
    assert not result['source_selected_capacity_or_repair_duration']
    assert not result['global_raw_archive_equality']
    assert not result['relativistic_continuum_or_physical_particle_assignment']
