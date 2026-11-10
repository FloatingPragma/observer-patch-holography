"""Collision checks bind the small exact block to the full source action."""
import sympy as s
import pytest

from .dark_collision import (actual_matrix,certificate,expected_matrix,
                            singlet_and_doublon,transfer_probability,source_column)


def test_full_source_intertwiner_and_singlet_block():
    matrix=actual_matrix();singlet,doublon=singlet_and_doublon()
    assert matrix==matrix.H
    assert matrix==expected_matrix(s.pi/4)
    assert matrix*singlet==(4+s.pi/2)*singlet-s.pi*doublon/2


def test_complete_collision_records_and_preparation_failures():
    result=certificate()
    conditional=[s.sympify(row['probability']) for row in result['actual_source_records_at_transfer_maximum']]
    original=[s.sympify(value) for value in result['peak_original_trial_occupation_probabilities']]
    assert s.simplify(sum(conditional)-1)==0
    assert s.simplify(sum(original)+s.Rational(3,4)-1)==0
    assert result['polarized_internal_pair_double_occupation_probability']=='0'
    assert not result['asymptotic_scattering_matrix_or_cross_section']
    assert not result['intercell_regional_topological_charge_transfer']
    counter=result['additive_onsite_countercontrol']
    assert counter['at_t_pi_over_four_and_tau_one_free_double_probability']=='1'
    assert s.simplify(1-s.sympify(counter['same_time_interacting_double_probability_upper_bound'])
                      -s.sympify(counter['minimum_complete_detector_contrast']))==0


def test_zero_hopping_and_time_controls():
    tau=s.Symbol('tau',real=True)
    assert transfer_probability(s.Integer(0),tau)==0
    assert transfer_probability(s.pi/4,s.Integer(0))==0
    with pytest.raises(ValueError):source_column((0,0))
    with pytest.raises(ValueError):source_column((0,4))
