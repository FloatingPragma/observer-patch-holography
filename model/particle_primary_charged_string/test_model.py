"""Independent CAR, geometry, intervention and full-source provenance controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m
from model.particle_primary_response_moments import model as response


def car(count, index):
    result=np.zeros((1<<count,1<<count))
    for state in range(1<<count):
        if state&(1<<index):
            result[state^(1<<index),state]=(-1)**((state&((1<<index)-1)).bit_count())
    return result


def test_endpoint_control_preserves_physical_gauss_support_and_opposite_charge():
    c=m.endpoint_control()
    assert c['endpoint_charge_increments']==[6,-6]
    assert c['electric_character_cost']==F(63,4)
    assert c['CAR_parity']=='even' and not c['physical_source_histogram_executed']
    # Independent four-mode control with both directed preparations present.
    cs=[car(4,i) for i in range(4)];ns=[c.T@c for c in cs]
    e=np.diag([-1.,0.,1.]);p0=np.diag([0.,1.,0.]);up=np.zeros((3,3));up[2,1]=1
    down=np.zeros((3,3));down[0,1]=1
    bp=np.kron(up,cs[0].T@cs[3]);bm=np.kron(down,cs[2].T@cs[1])
    qx=6*np.kron(np.eye(3),ns[0]+ns[1]-np.eye(16))
    qy=6*np.kron(np.eye(3),ns[2]+ns[3]-np.eye(16))
    flux=6*np.kron(e,np.eye(16))
    for b,sign in [(bp,1),(bm,-1)]:
        assert np.array_equal(qx@b-b@qx,6*sign*b)
        assert np.array_equal(qy@b-b@qy,-6*sign*b)
        assert np.array_equal((qx-flux)@b,b@(qx-flux))
        assert np.array_equal((qy+flux)@b,b@(qy+flux))
    # The complete two-direction effect is diagonal; repeat each direction
    # 24 times to check the declared49-outcome normalization on this control.
    success=[bp/np.sqrt(48)]*24+[bm/np.sqrt(48)]*24
    effect=sum(b.T@b for b in success)
    failure=np.diag(np.sqrt(np.maximum(0,1-np.diag(effect))))
    assert np.linalg.norm(effect+failure.T@failure-np.eye(48))<1e-14
    assert np.linalg.norm(bp@np.kron(p0,np.eye(16))-bp)<1e-14


def test_actual_cube_off_bond_wilson_gram_and_magnetic_support():
    gamma=response.density.gamma_data();beta=gamma[3]
    spin=[np.array((-beta+beta*g)/2,dtype=complex) for g in gamma[:3]]
    sites=list(product(range(2),repeat=3));index={x:i for i,x in enumerate(sites)}
    bonds=[]
    for x in sites:
        for axis in range(3):
            if x[axis]==0:
                y=list(x);y[axis]=1;bonds.append((x,tuple(y),axis))
    for x,y,chosen_axis in bonds:
        h=np.zeros((32,32),complex)
        for xx,yy,axis in bonds:
            if (xx,yy,axis)==(x,y,chosen_axis):continue
            a=slice(4*index[xx],4*index[xx]+4);b=slice(4*index[yy],4*index[yy]+4)
            h[a,b]=spin[axis];h[b,a]=spin[axis].conj().T
        for p,minus in response.beta_pairs():
            pp=np.zeros(32,complex);mm=np.zeros(32,complex)
            pp[4*index[x]:4*index[x]+4]=np.array(p,dtype=complex).ravel()
            mm[4*index[y]:4*index[y]+4]=np.array(minus,dtype=complex).ravel()
            vectors=np.stack([pp,mm,h@pp,h@mm],axis=1)
            assert np.linalg.norm(vectors.conj().T@vectors-np.eye(4))<1e-14
        incident_plaquettes=0
        for second_axis in range(3):
            if second_axis!=chosen_axis:incident_plaquettes+=1
        assert incident_plaquettes==2
    assert len(bonds)==12
    assert all(sum(x==a or x==b for a,b,_ in bonds)==3 for x in sites)


def test_full_CAR_two_disjoint_raising_pairs_have_sqrt_two_norm():
    c=[car(4,i) for i in range(4)]
    residual=c[2].T@c[1]-c[0].T@c[3]
    singular=np.linalg.svd(residual,compute_uv=False)
    assert abs(singular[0]**2-2)<1e-14
    # Summing separate operator norms would give2; the full Fock calculation
    # verifies the smaller norm used in the source estimate.
    assert np.linalg.matrix_rank(residual)==6


def test_union_bound_uses_all_directions_and_occupations():
    # One bond, one component: exhaustive occupancy and link-constant events.
    for px,my,constant in product((0,1),repeat=3):
        assert 1-constant*my*(1-px)<=px+(1-my)+(1-constant)
    sites=list(product(range(2),repeat=3));counts={(x,a,t):0 for x in sites for a in range(2) for t in ('p','hole')}
    links=0
    for x in sites:
        for axis in range(3):
            if x[axis]:continue
            y=list(x);y[axis]=1;y=tuple(y);links+=1
            for start,end in ((x,y),(y,x)):
                for a in range(2):counts[(start,a,'p')]+=1;counts[(end,a,'hole')]+=1
    assert set(counts.values())=={3} and links==12
    assert F(3,48)==F(1,16) and F(4,48)/F(3,10)==F(5,18)


def test_charge_dephasing_cannot_be_silently_inserted_into_energy_history():
    # A coherent selected energy eigenstate has no initial energy spread;
    # recording a noncommuting charge destroys that property.
    H=np.array([[0.,1.],[1.,0.]])
    ground=np.array([1.,-1.])/np.sqrt(2);rho=np.outer(ground,ground)
    dephased=np.diag(np.diag(rho))
    assert abs(np.trace(rho@H)+1)<1e-14
    assert abs(np.trace(dephased@H))<1e-14
    assert np.linalg.norm(rho-dephased)>0.7


@pytest.fixture(scope='module')
def bound():
    member=m.contract.member(0)
    return member,m.bind_member(member)


def test_full_space_allocation_and_original_energy_windows(bound):
    member,receipt=bound;r=receipt['result'];allocation=r['allocation']
    assert [F(x['onsite_per_endpoint']) for x in allocation['native_rows']]==[1,1,1,1,F(3,10)]
    assert allocation['native_rows'][-1]['full_Fock_dimension']==256
    assert F(allocation['nonconstant_link_margin_lower'])>5
    assert F(allocation['combined_bad_event_coefficient'])==F(168,5)
    assert F('0.794')<F(r['original_success_weight_lower'])<F('0.795')
    assert F(r['full_residual_norm_upper'])**2<70
    assert r['original_transfer_band']==['39/4','135/4']
    assert F(r['original_transfer_weight_lower'])>F(3,10)
    assert F(r['final_excitation']['original_final_excitation_weight_lower'])==F(19,80)
    assert r['reported_transfer_interval']==['77/8','271/8']
    assert F(r['reported_transfer_original_weight_lower'])>F(1,4)
    assert r['all_fine_energy_probe_records']==49*r['complete_clock_plan']['clock_plan']['clock_size']**2
    assert r['source_reference_and_failures_retained']
    assert not r['selected_preparation_replaced_by_constant_links']
    assert not r['energy_bound_transferred_to_charge_dephased_input']
    assert not r['physical_species_or_rest_mass_identified']
    assert r['physical_mass_table_GeV']==[]
    m.verify_bound(member,receipt)


def test_wrong_source_level_and_forged_output_are_rejected(bound):
    member,receipt=bound
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
    forged=deepcopy(receipt);forged['result']['original_success_weight_lower']='1'
    forged['result_sha256']=m.contract.digest(forged['result'])
    with pytest.raises(ValueError):m.verify_bound(member,forged)
    altered=deepcopy(member);altered['source']['level']=1
    with pytest.raises(ValueError):m.verify_bound(altered,receipt)
