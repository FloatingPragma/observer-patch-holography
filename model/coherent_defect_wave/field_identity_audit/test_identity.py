import sympy as s

from .model import (dark_module,face_charge_obstruction,bound_transport,antisymmetric_clouds,
                    mass_freedom,dark_code,current_action,representation,graded_swap)


def test_actual_dark_module_and_all_native_current_columns():
    result=dark_module()
    assert result['factorized_full_exterior_and_actual_control_columns_checked']==180
    assert result['actual_native_current_intertwiner_columns_checked']==180
    assert result['extra_native_degree_current_columns_checked']==15
    assert result['extra_degree_if_gauged_as_all_left_Weyl_anomaly_sums']=={'linear':40,'cubic':400}
    assert result['native_K_on_module']=='0'
    assert result['current_intertwiner_from_empty_native_flux_to_R15_dimension']==0


def test_a_charged_native_cloud_is_not_the_empty_native_flux():
    subsets,g5,currents,_=representation()
    state=dark_code(1,subsets[0])
    image=current_action(state,g5[-1])
    assert s.simplify(image.norm_squared()-currents[-1][0,0]**2)==0
    assert image.norm_squared()>0


def test_arbitrary_pure_face_current_is_scalar_under_actual_coupling():
    result=face_charge_obstruction()
    assert result['full_pure_face_current_constraint_rank']==24
    assert result['full_pure_face_current_commutant_dimension']==1
    assert len(result['all_three_plus_two_face_splits_fail'])==10
    assert all(s.Rational(row['actual_K_charge_commutator_norm_squared'])>0
               for row in result['all_three_plus_two_face_splits_fail'])


def test_complete_transport_keeps_native_charge_and_exchange_sign():
    result=bound_transport()
    assert result['total_native_current_commutator_columns_checked']==3072
    assert sum(map(s.Rational,result['complete_256_outcome_probabilities']))==1
    assert graded_swap({(1,1):s.Integer(1)})=={(1,1):-1}
    assert not result['full_unconstrained_composite_CAR_derived_from_two_slots']


def test_whole_cloud_antisymmetrization_intertwines_actual_currents():
    result=antisymmetric_clouds()
    assert result['actual_sum_of_native_currents_intertwiner_columns_checked']==1260
    assert not result['fixed_cloud_per_internal_CAR_mode_slot_is_a_current_intertwiner']


def test_internal_weights_do_not_select_dirac_mass_values():
    result=mass_freedom()
    assert result['gauge_invariant_Dirac_mass_parameter_dimension']==5
    assert result['constant_chiral_R15_Weyl_mass_dimension']==0
    assert result['all_corresponding_internal_rest_frequency_ratios']=='1'
    assert result['face_count_coupling_controls'][0]['rest_frequency']!=result['face_count_coupling_controls'][1]['rest_frequency']
