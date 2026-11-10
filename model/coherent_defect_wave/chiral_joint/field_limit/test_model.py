"""Exact full-source symmetries and controlled low-energy graph tests."""
import pytest
import sympy as s
from .model import (star_parity_checks,path_fluctuation_checks,quartic_control,
    reduction_bounds,physical_graph_matrices,connected_effective_energy_bound,sqrt_upper,
    bounded_transporter_control)


def test_actual_operator_columns_preserve_the_dressed_local_parities():
    data=star_parity_checks()
    assert data['actual_hopping_columns']==24576
    assert data['blank_link_effective_H_preserves_each_site_matter_parity']
    assert not data['single_fermion_hopping_with_all_links_returned_blank']
    assert not data['dressed_charged_propagation_excluded']


def test_actual_path_words_have_the_declared_block_dependent_fluctuation():
    data=path_fluctuation_checks()
    assert data['native_block_dimensions']==[3,6,1,2,3]
    for row in data['rows']:
        d,r=row['block_dimension'],row['length']
        assert s.Integer(row['fixed_endpoint_path_norm_squared'])==d**(r-1)
    assert not data['raw_path_equals_unitary_transport_on_blank_reference']


def test_low_graph_includes_all_physical_states_and_an_exact_metric():
    q=physical_graph_matrices()
    assert q['H'].shape==(15,15) and q['X'].shape==(9,6)
    assert q['effective_form']==q['effective_form'].T
    assert q['effective_metric']*q['effective_generator']==q['effective_form']
    bound=sqrt_upper(q['residual_squared_Hilbert_Schmidt'])
    assert 0<bound<s.Rational(1,10**12)
    assert q['leading_virtual_interaction']!=s.zeros(6)
    # Dropping the actual coupling gives an exactly invariant blank block,
    # while this source has a strictly nonzero omitted component.
    assert q['B']!=s.zeros(9,6)


def test_controlled_virtual_interaction_has_positive_connected_energy():
    assert quartic_control()['connected_diagonal_coefficient_of_minus_B_adjoint_B']=='1/4'
    data=connected_effective_energy_bound()
    assert s.Rational(data['exact_dressed_connected_energy_lower_bound'])>0
    assert s.Rational(data['canonical_effective_H_leading_error'])>0
    with pytest.raises(ValueError):reduction_bounds(1)
    with pytest.raises(ValueError):reduction_bounds(128.0)
    with pytest.raises(ValueError):reduction_bounds(128,-1)


def test_iteration_and_clock_errors_have_the_required_scaling():
    previous=reduction_bounds(128,2,time=2)
    following=reduction_bounds(128,3,time=2)
    assert s.Rational(following['complete_effect_error'])<s.Rational(previous['complete_effect_error'])
    unit=reduction_bounds(128,3,time=1)
    assert s.Rational(following['complete_effect_error'])==2*s.Rational(unit['complete_effect_error'])


def test_explicit_normalization_controls_path_norm_without_transferring_a_field_limit():
    result=bounded_transporter_control()
    assert result['all_path_products_are_contractions']
    assert result['blank_reference_path_norm_squared']=='1'
    assert not result['actual_raw_periodic_family_modified']
    assert not result['interacting_relativistic_limit_proved']
