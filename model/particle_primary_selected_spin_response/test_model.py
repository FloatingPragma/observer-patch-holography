"""Exact local CAR, actual cube labels, coherent current and complete clocks."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_actual_Gauss_kernel_and_all_local_occupations_supply_complete_effects():
    c=m.local_controls();assert c['active_source_modes']==(6,10,11,21,25,26)
    assert set(c['families'])=={'G1','H'}
    for label,row in c['families'].items():
        d=2 if label=='G1' else 4;assert row['dimension']==d
        assert row['effect'][0,0]==d
        assert all(x>=0 for x in row['deficit_spectrum'])
        for mask in range(64):
            if mask.bit_count()==1:assert row['effect'][mask,mask]==s.Rational(d,2)
        for A in row['operators']:
            assert (A[:,0].H*A[:,0])[0]==1
    E=c['total_success_effect'];assert all(0<=x<=1 for x in E.eigenvals())
    failure=s.eye(64)-E
    assert failure*failure!=failure  # the actual failure must be its positive square root


def test_twelve_mode_residual_has_exact_small_block_rational_certificate():
    for label,upper,den,blocks in [('G1',F(19,3),6,1922),('H',F(12),3,1858)]:
        row=m.residual_control(label)
        assert row['squared_residual_norm_upper']==upper
        assert row['integer_Gram_denominator']==den and row['connected_blocks']==blocks
        assert row['largest_exact_block']==32 and row['smallest_positive_LDL_pivot']>0
        assert row['exact_rational_factorization_checked']
    counter=m.residual_control('G1')['three_times_dimension_shortcut_countercontrol']
    assert F(counter['Rayleigh_numerator'],counter['Rayleigh_denominator'])>6
    assert counter['domain'].startswith('uncompressed twelve-mode')


def test_actual_native_Wilson_images_are_jointly_orthogonal_in_all_six_modes():
    source=m.odd.sectors.branch.response;gamma=source.density.gamma_data();chart=m.spin.beta_chart()
    # Reuse the proved native/coordinate chart permutation, and a nonreal actual link.
    controls=m.odd.exact_controls();permutation=controls['native_to_coordinate_basis_permutation']
    z=(3+4*s.I)/5;u=s.Matrix([[s.Rational(3,5),4*s.I/5],[4*s.I/5,s.Rational(3,5)]])
    fundamental=s.diag(z**-2*s.diag(u,s.ones(1)),z**3*u)
    raw=source.density.native_representation(fundamental);R=raw.extract(permutation,permutation)
    assert m.clean(R.H*R)==s.eye(15)
    active=m.local_controls()['active_source_modes'];original=s.eye(60)[:,list(active)]
    for axis in range(3):
        hop=(-gamma[3]+gamma[3]*gamma[axis])/2
        for S,G in ((hop,R),(hop.H,R.H)):
            A=s.kronecker_product(m.clean(chart.H*S*chart),G)[:,list(active)]
            assert m.clean(A.H*A)==s.eye(6)/2
            left=s.Matrix.vstack(original,s.zeros(60,6));right=s.Matrix.vstack(s.zeros(60,6),A)
            assert left.H*right==s.zeros(6)


def test_actual_binary_cube_twirl_distinguishes_probe_irreps_and_component_laws():
    c=m.local_controls();cube=m.spin.cube_restriction(m.DELTA)['matrices'];reps={}
    for label,row in c['families'].items():
        W=row['coefficient_basis'];reps[label]=[m.clean(W.H*D*W) for D in cube]
        d=W.cols;X=s.Matrix(d,d,lambda i,j: (i+1)*(j+2)+s.I*(i-j))
        average=m.clean(sum((D.H*X*D for D in reps[label]),s.zeros(d))/48)
        assert average==s.trace(X)*s.eye(d)/d
    X=s.ones(2,4)
    cross=m.clean(sum((A.H*X*B for A,B in zip(reps['G1'],reps['H'])),s.zeros(2,4))/48)
    assert cross==s.zeros(2,4)
    control=m.covariance_controls()
    assert not control['initial_selected_state_is_assumed_a_cube_singlet']
    assert not control['irreducible_probe_label_is_final_total_state_label']


def test_native_current_increment_and_complete_three_state_meter_alphabet():
    c=m.current_meter_control()
    assert c['meter_count']==5 and c['complete_pointer_dimension']==243
    assert c['two_coupling_identity_preserves_full_source_reference_histories']
    assert not c['initial_absolute_numbers_measured']
    assert c['extra_transfers_outside_alphabet_could_alias']


def test_coherent_transfer_meter_preserves_entangled_source_reference_history():
    # A literal unknown-number superposition; the pointer's difference is known.
    dimension=4;radix=3
    K=np.zeros((dimension,dimension),complex)
    for n in range(dimension-1):K[n+1,n]=.5
    complement=np.diag(np.sqrt(1-np.diag(K.conj().T@K).real))
    psi=np.diag(np.array([1,2,3,1j])/np.sqrt(15))
    for operator,delta in ((K,1),(complement,0)):
        output=np.zeros((dimension,dimension,radix),complex)
        for initial,reference in product(range(dimension),repeat=2):
            before=(-initial)%radix
            for final in range(dimension):output[final,reference,(before+final)%radix]+=operator[final,initial]*psi[initial,reference]
        expected=np.zeros_like(output);expected[:,:,delta]=operator@psi
        assert np.linalg.norm(output-expected)<1e-14
    assert np.linalg.norm(psi-np.diag(np.diag(psi)))==0
    density=np.outer(psi.reshape(-1),psi.reshape(-1).conj())
    assert abs(density[0,5])>.1  # reference coherence exists and is not initially measured


@pytest.fixture(scope='module')
def bound_members():return [m.bind_member(m.contract.member(n)) for n in (0,1,2,8)]


def test_actual_members_have_disjoint_positive_original_G1_and_H_events(bound_members):
    for value in bound_members:
        row=value['result'];assert F(row['complete_two_clock_probability_error_upper'])<=F(1,4096)
        assert row['complete_address_and_success_failure_outcomes']==7*row['source_site_count']
        total=F()
        for channel in row['probe_channels']:
            simple=F(1,200) if channel['probe_irrep']=='G1' else F(3,200)
            assert F(channel['finite_record_original_weight_lower'])>simple
            total+=F(channel['original_success_weight_lower'])
            assert F(channel['each_component_reported_weight_lower'])*channel['dimension']==F(channel['finite_record_original_weight_lower'])
            assert channel['weights_are_original_and_not_success_normalized']
        assert 0<total<=1
        assert row['source_input_is_original_selected_density'] and not row['source_preparation_replaced_by_anchor_or_trial']
        assert row['all_failures_memory_and_reference_retained'] and row['physical_mass_table_GeV']==[]


def test_member_zero_sharp_band_and_corrected_moment_constants(bound_members):
    rows=bound_members[0]['result']['probe_channels']
    for row,constant,weight in zip(rows,(F(19,12),F(3)),(F(121,512),F(61,128))):
        assert F(row['centered_second_moment_upper'])==constant
        narrow=row['member_zero_sharp_record']
        assert narrow['original_transfer_interval']==[3,15]
        assert list(map(F,narrow['finite_reported_transfer_interval']))==[F(11,4),F(61,4)]
        assert F(narrow['finite_reported_original_weight_lower'])>weight


def test_complete_spatial_pointer_rotation_preserves_aggregate_not_each_row():
    c=m.spatial_analyzer_control()
    assert c['all_spatial_rows_retained'] and c['aggregate_each_internal_spin_channel_weight_unchanged']
    assert not c['a_fixed_spatial_row_positive_weight_proved']
    assert not c['spatial_label_is_conserved_momentum']
    assert not c['prepared_anchor_transport_moment_transferred_to_selected_spectrum']
    rng=np.random.default_rng(381);K=rng.normal(size=(8,2,2))+1j*rng.normal(size=(8,2,2))
    cells=list(product(range(2),repeat=3));W=np.array([[(-1)**sum(a*b for a,b in zip(x,y)) for x in cells] for y in cells])/np.sqrt(8)
    transformed=np.einsum('ab,bij->aij',W,K)
    assert np.linalg.norm(sum(A.conj().T@A for A in K)-sum(A.conj().T@A for A in transformed))<1e-12
    assert abs(np.linalg.norm(K[0])**2-np.linalg.norm(transformed[0])**2)>.01


def test_binding_rejects_changed_clock_and_false_state_spin_assignment(bound_members):
    member=m.contract.member(0);m.verify_bound(member,bound_members[0])
    value=deepcopy(bound_members[0]);value['result']['covariance']['irreducible_probe_label_is_final_total_state_label']=True
    with pytest.raises(ValueError):m.verify_bound(member,value)
    changed=deepcopy(member);changed['source']['spacing']='2'
    with pytest.raises(ValueError):m.bind_member(changed)


def test_certificate_is_a_canonical_JSON_value():
    value=m.certificate()
    assert json.loads(json.dumps(value,sort_keys=True))==value
    assert len(value["members"])==4
