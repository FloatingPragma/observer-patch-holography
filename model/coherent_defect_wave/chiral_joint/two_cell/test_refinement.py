"""Exact source transfer, joint apparatus momentum and dispersion controls."""
import pytest
import sympy as s
from .refinement import (exact_transfer_control,compensated_momentum_control,
    primitive_link_variance_control,momentum_design,shell_feasibility,rational_error_budget)


def test_actual_physical_input_exposes_hidden_blank_link_leakage():
    value=exact_transfer_control()
    assert value['compressed_residual_squared']=='0'
    assert value['full_residual_squared']=='3/16'
    assert s.sympify(value['actual_new_link_yes_probability_lower_bound']).is_positive
    assert value['embedded_open_new_link_yes_probability']=='0'
    assert not value['old_ground_or_spectrum_transferred_to_periodic']


def test_joint_compensated_momentum_is_complete_and_Gauss_physical():
    value=compensated_momentum_control()
    assert value['complete_probabilities']==['0','1/2','0','1/2']
    assert value['all_24_Gauss_constraints_hold_for_input_and_every_phase_branch']
    assert value['phase_observable']=='joint source/reference translation'
    assert not value['source_only_particle_momentum_identified']
    assert not value['source_input_is_stationary_ground']


def test_nonzero_classical_background_needs_more_than_fixed_link_memory():
    value=primitive_link_variance_control()
    assert not value['nonzero_c_number_transfer_can_have_vanishing_relative_residual']
    assert 'dressed/interacting limits are not excluded' in value['scope']


def test_more_momentum_samples_are_needed_to_test_a_two_parameter_shell():
    assert momentum_design(2)['mass_squared_speed_squared_design_rank']==1
    assert momentum_design(4)['mass_squared_speed_squared_design_rank']==2
    assert momentum_design(4)['independent_mass_shell_tests_after_two_parameter_fit']==0
    assert momentum_design(6)['independent_mass_shell_tests_after_two_parameter_fit']==1
    good=shell_feasibility([(1,5,5),(9,13,13),(25,29,29)])
    assert good['parameter_vertices']==[['4','1']]
    assert not shell_feasibility([(1,5,5),(9,13,13),(25,30,30)])['feasible']
    with pytest.raises(ValueError):shell_feasibility([(1,1.0,2)])


def test_complete_comparison_retains_time_tail_and_detector_errors():
    assert rational_error_budget(s.Rational(1,10),1,s.Rational(15,2))==s.Rational(3,4)
    assert rational_error_budget(0,1,0,tail=s.Rational(1,10))==s.Rational(1,5)
    assert rational_error_budget(2,1,1)==1
    with pytest.raises(ValueError):rational_error_budget(1,0,1)
