from fractions import Fraction as F
import math
import pytest
import sympy as s

from model.native_wilson_maxwell import model as m


@pytest.fixture(scope='module')
def result():return m.certificate()


def test_open_cube_gap_identity_does_not_assume_commuting_derivatives():
    H,_,_,_=m.matter.cube_matter();_,beta=m.current.matrices()
    beta=s.kronecker_product(s.eye(8),beta)
    assert (beta*H+H*beta)/2==s.kronecker_product(s.eye(8)+m.field.cube()['D'].T*m.field.cube()['D']/2,s.eye(4))
    assert beta*H!=H*beta
    assert (m.field.cube()['D'].T*m.field.cube()['D']).rank()==7


def test_no_ground_projection_or_missing_trial_probability(result):
    p=result['actual_preparation'];support=result['literal_trial_support']
    assert p['base_trial_success']=='1' and p['native_trial_copies']==32
    assert p['filter_depth']==0 and p['agreement_filter_count']==0
    assert p['no_initial_P0_or_code_measurement']
    assert F(p['global_code_probability_lower'])>F(999999,10**6)
    assert support['four_copy_all_reference_vacuum_intersection_rank']==1


def test_Fock_charge_frame_uses_the_declared_one_body_orientation(result):
    x=result['full_charge_frame_convention']
    assert x['exact_full_pair_Fock_rotation_checked']
    assert x['wrong_inverse_density_action_explicitly_rejected']
    assert x['density_one_body_action']=='M_site^dagger P_x M_site'
    assert x['exact_32_dimensional_grouped_site_permutation_checked']
    assert x['omitted_column_permutation_explicitly_rejected']


def test_complete_count_ranks_include_every_Fock_and_photon_state(result):
    rank=sum(math.comb(16,a)*math.comb(16,b)*math.comb(n+4,4)
             for a in range(17) for b in range(17) for n in range(4))
    assert rank==2**32*56
    assert result['count_instrument']['physical_count_alphabet_size']==1160
    assert result['count_instrument']['full_pointer_dimension']==1161
    assert result['full_Gauss_reference']['full_longitudinal_dimension']==9**7
    assert result['full_Gauss_reference']['full_static_reference_dimension']==65


def test_two_independent_current_coefficients_and_full_cap_boundary(result):
    coefficient=result['independent_current_coefficient']
    assert coefficient['strict_outward_interval_inclusion']
    exact=coefficient['exact_CRT_method']
    assert exact['complete_CRT_projector_ranks']==[8,16,8]
    assert not exact['binomial_sign_approximation_or_numerical_eigensolver_used']
    assert result['field_source']['photon_dimension']==56
    assert result['field_source']['total_photon_cutoff']==3


def test_full_Dyson_and_actual_complete_connected_probability(result):
    r=result['complete_interaction_records'];x={k:F(v) for k,v in r['exact'].items()}
    assert r['target_count']==[1,1,1] and r['free_target_probability']=='0'
    assert x['actual_unconditional_target_probability_lower']>F(475,10**10)
    assert x['actual_unconditional_binary_connected_response_lower']>F(471,10**10)
    assert x['entire_operator_twice_Duhamel_remainder_upper']<x['leading_target_amplitude_lower']/20
    assert x['independent_depth1_conditional_connected_response_lower']>0
    assert r['all_virtual_matter_and_photon_occupations_in_Duhamel_bound']
    assert not r['source_state_is_actual_coupled_ground']


def test_same_full_clock_capacity_and_global_code_scope(result):
    c=result['clock_and_capacity'];lo,hi=map(F,c['full_native_source_reference_trial_capacity_nats_interval'])
    assert 10000<lo<hi<12000
    assert c['source_constituents']==32
    assert c['global_exact_code_dimension']==2**32
    assert c['global_spectral_code_compression_is_declared_architecture']
    assert 'sqrt(p0^32)' in c['arbitrary_pulse_error_target_probability_lower']
    assert not result['boundaries']['all_state_native_spatial_locality_transferred']


def test_public_custody_replays_both_independent_research_references(result):
    c=result['independent_research_comparison']
    assert c['all_scientific_fields_equal']
    assert c['reference_sha256']=='2fcb312421565591e3a62308986b98d3803519320949185ef72e9c8f19c4495e'
    assert c['shared_field_reference_sha256']=='f646ac992fc589a042b1ddb635f6d7f288943d7bec15aeca846badda4608cb7e'
    assert c['both_independent_current_methods_and_all_source_record_bounds_replayed']
    m.filtered.higgs.validate(result['input_sha256'])
