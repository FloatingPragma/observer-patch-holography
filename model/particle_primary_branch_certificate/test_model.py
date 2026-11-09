"""Independent CAR, full quadratic form, symmetry and selected-floor controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_literal_CAR_move_against_independent_exterior_action():
    # Independent full64-state Jordan-Wigner matrices, including blocked moves.
    n=6;ann=[]
    for mode in range(n):
        A=np.zeros((2**n,2**n),dtype=np.int64)
        for mask in range(2**n):
            if mask>>mode&1:
                A[mask^(1<<mode),mask]=(-1)**sum(mask>>j&1 for j in range(mode))
        ann.append(A)
    for i,j in product(range(n),repeat=2):
        if i==j:continue  # Source edges always join distinct modes.
        actual=ann[i].T@ann[j]
        for mask in range(2**n):
            moved=m._move(mask,i,j)
            column=np.zeros(2**n,dtype=np.int64)
            if moved:column[moved[0]]=moved[1]
            assert np.array_equal(column,actual[:,mask])
    assert m._gaussian(7-13j)==(7,-13)
    with pytest.raises(ArithmeticError):m._gaussian(.5+0j)
    with pytest.raises(ArithmeticError):m._gaussian(complex(2**20,0))


def test_all_admitted_two_link_packets_have_correct_Haar_geometry():
    packet=m.source_packets();pairs=packet['two_distinct_link_pairs'];groups=packet['second_action_groups']
    assert len(pairs)==66 and sum(r['share_vertex'] for r in pairs)==24
    assert {r['norm_squared'] for r in pairs}=={F(1785,2),F(900)}
    assert len(groups)==30
    assert sum(r['norm_squared'] for r in groups)==59220
    assert sum(r['norm_squared']*r['relative_H0_energy'] for r in groups)==F(11350341,2)
    assert all(x['imaginary']==0 for r in groups for x in r['one_packet_couplings'])
    assert sum(x['real'] for r in groups for x in r['one_packet_couplings'])==59220
    assert len(packet['source_first_packet_norms'])==7
    assert packet['actual_source_modes']==480
    assert packet['all_first_and_second_actions_within_actual_product_cutoff']
    assert not packet['zero_net_matter_profile_in_second_action']


def _dense_trial():
    packets=m.source_packets();trial=m.trial_certificate()
    first=packets['source_first_packet_norms'];groups=packets['second_action_groups']
    n=np.array([float(x['norm_squared']) for x in first])
    N=np.array([float(x['norm_squared']) for x in groups])
    cas=packets['native_Casimir_order'];k=len(first);q=len(groups)
    g=float(F(6,22050));A=np.zeros((1+k+q+5,)*2)
    v=np.ones(A.shape[0]);v[1:1+k]=-np.sqrt(n)*np.array(trial['one_packet_real_coefficients'],float)
    v[1+k:1+k+q]=-np.sqrt(N)*np.array(trial['second_action_real_coefficients'],float)
    v[-5:]=-np.sqrt(g)*np.array(trial['loop_real_coefficients'],float)
    A[0,1:1+k]=np.sqrt(n);A[1:1+k,0]=np.sqrt(n)
    for i,row in enumerate(first):
        r,d=row['packet'];A[1+i,1+i]=float(F(105,2)*cas[r]+3*d)
    for j,row in enumerate(groups):
        A[1+k+j,1+k+j]=float(row['relative_H0_energy']+F(6,7))
        for i,firstrow in enumerate(first):
            c=next(r['real'] for r in row['one_packet_couplings'] if r['packet']==firstrow['packet'])
            A[1+i,1+k+j]=A[1+k+j,1+i]=float(c)/np.sqrt(n[i]*N[j])
    loop=m.sectors.branch.variational.character_control()['cubic_loop_singlet_matrix']
    A[-5:,-5:]=np.diag([210*float(c) for c in cas])-np.array(loop,float)/420
    A[0,-5:]=np.sqrt(g);A[-5:,0]=np.sqrt(g)
    return A,v,1+k,q


def test_exact_trial_matches_independent_dense_quadratic_form_and_bound_direction():
    A,v,start,count=_dense_trial();trial=m.trial_certificate()
    assert A.shape==(43,43)
    assert abs(v@v-float(trial['trial_norm_squared']))<1e-13
    quotient=v@A@v/(v@v)
    assert abs(quotient-float(trial['source_Rayleigh_shift_upper']))<1e-12
    assert np.linalg.eigvalsh(A)[0]<=quotient+1e-12
    assert -715.10020<float(trial['source_Rayleigh_quotient_upper'])<-715.10018
    # The certified two-link block is an upper form. An arbitrary actual
    # magnetic block bounded above by6/7 cannot reverse its direction.
    rng=np.random.default_rng(117)
    R=rng.normal(size=(count,count));R=(R+R.T)/2
    R=R/max(abs(np.linalg.eigvalsh(R)))*float(F(6,7))
    actual=A.copy();actual[start:start+count,start:start+count]+=R-float(F(6,7))*np.eye(count)
    assert v@actual@v<=v@A@v+1e-12
    assert not trial['trial_is_an_eigenvector']
    assert not trial['full_source_omitted_space_assumed_invariant']


def test_omitted_directions_can_lower_the_ground_below_a_Ritz_value():
    A,v,_,_=_dense_trial()
    enlarged=np.zeros((44,44));enlarged[:43,:43]=A
    enlarged[43,43]=100;enlarged[0,43]=enlarged[43,0]=20
    assert np.linalg.eigvalsh(enlarged)[0]<np.linalg.eigvalsh(A)[0]
    assert np.allclose(enlarged[:43,:43],A)


def test_exact_source_cube_matrix_units_and_odd_degeneracy_scope():
    out=m.symmetry_reduction()
    assert out['full_odd_sector_irrep_dimensions']=={'G1':2,'G2':2,'H':4}
    assert out['matrix_coefficient_orthogonality_pairs_checked']==576
    assert out['complete_odd_sector_diagonal_projector_outcomes']==8
    assert out['actual_local_diagonal_projector_ranks']=={'G1':[1,1],'G2':[0,0],'H':[1,1,1,1]}
    assert out['full_sector_rank_one_complement_bottom_equals_sector_bottom']
    assert not out['rank_one_obstruction_applies_after_diagonal_symmetry_reduction']
    # A multiplicity-two ground survives every full-space rank-one removal.
    H=np.kron(np.eye(2),np.diag([1.,4.,9.]))
    rng=np.random.default_rng(512)
    for _ in range(12):
        v=rng.normal(size=6)+1j*rng.normal(size=6);v/=np.linalg.norm(v)
        Q=np.linalg.qr(np.column_stack((v,np.eye(6))),mode='complete')[0][:,1:]
        assert abs(np.linalg.eigvalsh(Q.conj().T@H@Q)[0]-1)<1e-12
    # One diagonal representation matrix unit leaves the multiplicity space.
    assert np.linalg.eigvalsh(H[:3,:3][1:,1:])[0]==4


def test_polynomial_energy_local_floor_without_exact_Gibbs_substitution():
    # An independent finite positive Taylor-polynomial density has exact
    # rational weights. The source proof applies the same modular-slope bound.
    energies=[F(-2),F(-1),F(0),F(3)];beta=F(2);scale=4;order=3
    weights=[]
    for E in energies:
        y=beta*(max(energies)-E)/(2*scale)
        T=sum((y**j/F(s.factorial(j)) for j in range(order+1)),F())
        weights.append(T**(2*scale))
    total=sum(weights);weights=[w/total for w in weights]
    for U in energies:
        exponent=2+(2*beta*(U-min(energies))).__ceil__()
        assert all(p>=F(1,2**exponent) for E,p in zip(energies,weights) if E<=U)
    # Multiplying the actual positive effect preserves the tiny density
    # floor; subtracting an absolute implementation loss would destroy it.
    p=F(1,2**30);effect=F(511,512)
    assert p*effect>0
    assert p-F(1,512)<0


@pytest.fixture(scope='module')
def bound():
    member=m.contract.member(0)
    return member,m.bind_member(member)


def test_full_source_improvement_and_energy_local_original_event(bound):
    member,value=bound;r=value['result'];floor=r['original_selected_energy_local_density_floor']
    lo,hi=map(F,r['full_odd_sector_bottom_interval']);glo,ghi=map(F,r['full_source_ground_interval'])
    assert ghi<lo<hi<F(r['thirteen_coordinate_Rayleigh_upper'])
    assert F(r['Rayleigh_upper_improvement'])>F(31,100)
    assert F(r['Ritz_upper_minus_largest_method_ceiling'])>F(3,10)
    assert floor['negative_exponent']==552953
    assert floor['ambient_dimension_bits_upper']==696
    assert floor['global_physical_density_floor_negative_exponent']==582407133
    assert not floor['exact_Gibbs_substituted_for_selected_polynomial']
    event=r['complete_original_final_energy_event'];loss=F(event['one_complete_clock_loss_upper'])
    assert loss<F(1,512)
    assert F(event['per_signed_sector_weight_lower_factored']['rational_multiplier'])==2*(1-loss)
    assert F(event['either_signed_sector_weight_lower_factored']['rational_multiplier'])==4*(1-loss)
    assert not event['initial_energy_measurement_performed']
    assert not event['event_is_an_energy_transfer_or_particle_creation_record']
    assert not event['probability_divided_by_occupation_success']
    assert not r['analysis_trial_used_as_original_input']
    assert not r['isolated_branch_or_resolved_mass_proved']
    assert r['physical_mass_table_GeV']==[]


def test_source_or_clock_forgery_cannot_promote_the_bound(bound):
    member,value=bound;m.verify_bound(member,value)
    for field in ('full_odd_sector_bottom_interval','physical_mass_table_GeV'):
        bad=deepcopy(value);bad['result'][field]=['electron',1]
        bad['result_sha256']=m.contract.digest(bad['result'])
        with pytest.raises(ValueError):m.verify_bound(member,bad)
    bad=deepcopy(value)
    bad['result']['original_selected_energy_local_density_floor']['negative_exponent']=1
    bad['result_sha256']=m.contract.digest(bad['result'])
    with pytest.raises(ValueError):m.verify_bound(member,bad)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
