"""Complete field, filled-sea and mutation controls for one source."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy import linalg as la

from . import modes, backreaction, packet


def test_all_scalar_modes_and_same_graph():
    data=modes.build_scalar_source()
    assert len(data['frequencies'])==343
    assert data['modes'].shape==(343,343)
    assert np.array_equal(data['K'],np.eye(343)+data['laplacian']+4096*np.diag(data['defect_amplitude_squared']))
    assert len(data['core_sites'])==2
    assert max(data['errors'].values())<2e-8
    selected=modes.selected_response(data)
    assert np.allclose(selected['core_profile_values'],[.707106633672487]*2,atol=2e-12)
    assert selected['selected_profile_minimum']>0
    assert abs(selected['strong_to_weak_ratio']-1.000122033061)<2e-10
    assert np.allclose(selected['vacuum_spectral_normalized_strengths'],
                       np.sqrt(2*selected['frequency'])*np.array(selected['unit_ladder_strengths']))


def test_literal_cap_first_moment_and_retarded_metric():
    omega=np.array([[2.,.3],[.3,1.5]])
    h=modes.scalar_hamiltonian(omega)
    fields=[modes.scalar_field(x,0,sites=2) for x in range(2)]
    def comm(a,b):return a@b-b@a
    actual=np.array([[comm(a,comm(h,b.conj().T))[0,0] for b in fields] for a in fields])
    assert np.allclose(actual,2*omega)
    data=modes.build_scalar_source()
    z=.7+1.2j
    response=modes.vacuum_response(z,data)
    metric=data['scalar_field_to_canonical']
    assert la.norm(metric@response@metric-la.inv(z*z*np.eye(343)-data['K']))<2e-10
    with pytest.raises(ValueError):modes.vacuum_response(1.,data)


def test_texture_decoupled_control_is_uniform_without_strength_split():
    data=modes.build_scalar_source(False)
    response=modes.selected_response(data)
    assert np.array_equal(data['K'],np.eye(343)+data['laplacian'])
    assert response['frequency']==1
    assert abs(response['strong_to_weak_ratio']-1)<2e-12
    assert np.allclose(response['unit_ladder_strengths'],1/np.sqrt(343))


def test_native_singlet_coefficients_and_proved_pointer_bound():
    native=backreaction.native_coefficients()
    assert native['global_native_gauge_singlet']
    assert native['all_twelve_generator_traces_zero']
    assert native['native_internal_labels']==15
    assert native['literal_local_matrix_unit_words']==112
    assert native['vacuum_to_scalar_branch_coefficient_absolute_sums']==[14]*4
    assert native['local_interaction_norm_upper_exact']==28
    assert native['complete_Ward_checks']['full_twelve_current_pair_Ward_checks']==72


def test_independent_car_image_and_complete_original_distribution():
    control=backreaction.small_car_control()
    assert control['complete_particle_antiparticle_Gram_error']<2e-12
    data=backreaction.calculate_backreaction()
    assert abs(data['interaction_image_norm_squared']-29.0179791642)<2e-8
    rows=data['complete_scalar_mode_records'];pointer=data['complete_pointer']
    assert len(rows)==343 and [row['mode'] for row in rows]==list(range(343))
    assert abs(sum(row['conditional_first_image_probability'] for row in rows)-1)<2e-12
    assert abs(pointer['failure_probability']+sum(row['original_pointer_and_mode_probability'] for row in rows)-1)<2e-12
    assert abs(rows[1]['conditional_first_image_probability']-.964506586904)<2e-10
    assert abs(rows[0]['conditional_first_image_probability']-.000411882689)<2e-10
    assert pointer['norm_bound']==28*343
    assert pointer['success_probability']<1e-6
    assert data['reference']['native_fermion_number']==5145
    assert data['reference']['combined_source_reference_fermion_parity']=='even'


def test_all_native_first_image_and_strict_ritz_lowering():
    data=backreaction.calculate_backreaction()
    c2=data['interaction_image_norm_squared'];cost=data['normalized_first_image_H0_excitation']
    ritz=la.eigvalsh([[0,np.sqrt(c2)],[np.sqrt(c2),cost]])[0]
    assert abs(ritz-data['two_vector_Ritz_energy_relative_to_sea'])<2e-12
    assert abs(ritz+2.46785318204)<2e-8
    exact=data['strict_variational_bound']
    assert exact['complete_spatial_gap_certificate_regenerated']
    assert F(exact['certified_parent_core_abs_energy_lower_exact'])>F(1,1000)
    assert F(exact['ground_energy_lowering_at_least_exact'])>0
    assert not exact['exact_interacting_ground_or_pole_reconstructed']
    assert data['dominant_pair_mode']['strong_to_weak_ratio']<1+1e-9
    assert la.norm(data['kernel'].imag)<2e-10


def test_archive_mutations_reject_mode_truncation_and_gram_forgery():
    arrays=packet.archive_arrays();packet.validate_archive(arrays)
    bad=dict(arrays);bad['modes']=bad['modes'][:,:2]
    with pytest.raises(ArithmeticError):packet.validate_archive(bad)
    bad=dict(arrays);bad['pair_Gram']=bad['pair_Gram'].copy();bad['pair_Gram'][0,0]+=.01
    with pytest.raises(ArithmeticError):packet.validate_archive(bad)
    bad=dict(arrays);bad.pop('antiparticle_Gram')
    with pytest.raises(ArithmeticError):packet.validate_archive(bad)


def test_receipt_mutations_reject_postselection_and_generation_claim():
    fresh=packet.calculate();bad=deepcopy(fresh)
    bad['native_filled_sea_first_image']['complete_pointer']['success_probability']=1.
    assert not packet.compare(bad,fresh)
    bad=deepcopy(fresh);bad['native_filled_sea_first_image']['complete_pointer']['success_probability']*=1.01
    assert not packet.compare(bad,fresh)
    bad=deepcopy(fresh);bad['scope']['three_generations_derived']=True
    assert not packet.compare(bad,fresh)
    bad=deepcopy(fresh);bad['complete_scalar_frequencies']=bad['complete_scalar_frequencies'][:2]
    assert not packet.compare(bad,fresh)
