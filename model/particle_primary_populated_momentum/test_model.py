"""Source phase identities, complete physical controls and bound tampering."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import permutations, product

import numpy as np
import pytest
import sympy as s

from . import model as m


def car(count,index):
    result=np.zeros((1<<count,1<<count),complex)
    for mask in range(1<<count):
        if mask&(1<<index):
            result[mask^(1<<index),mask]=(-1)**((mask&((1<<index)-1)).bit_count())
    return result


def test_actual_cube_phase_residual_for_arbitrary_nonconstant_links_and_q():
    beta=m.staggered.odd.response.beta_pairs()
    chart=s.Matrix.hstack(beta[0][0],beta[1][0],beta[0][1],beta[1][1])
    gamma=m.staggered.odd.response.density.gamma_data()
    hops=[np.array(chart.H*((-gamma[3]+gamma[3]*g)/2)*chart,dtype=complex) for g in gamma[:3]]
    cells=m.CELLS;index={x:i for i,x in enumerate(cells)}
    h=np.zeros((32,32),complex)
    for x in cells:
        for axis in range(3):
            if x[axis]:continue
            y=list(x);y[axis]=1;y=tuple(y);i=4*index[x];j=4*index[y]
            phase=np.exp(1j*(index[x]+2*axis+1)/3)
            h[i:i+4,j:j+4]=phase*hops[axis]
            h[j:j+4,i:i+4]=phase.conjugate()*hops[axis].conj().T
    plus=[4*x+j for x in range(8) for j in (0,1)]
    minus=[4*x+j for x in range(8) for j in (2,3)]
    for q in [(0,0,0),(np.pi,)*3,(np.pi/16,)*3,(.21,.8,2.1),(np.pi-.25,np.pi+.25,np.pi)]:
        raising=np.zeros((32,32),complex)
        for x in cells:
            i=4*index[x];raising[i:i+2,i+2:i+4]=np.exp(-1j*np.dot(q,x))*np.eye(2)
        residual=h@raising-raising@h
        creation=residual[np.ix_(plus,minus)]
        assert abs(np.linalg.norm(creation)**2-16*sum(np.cos(t/2)**2 for t in q))<2e-13
        assert np.linalg.norm(residual[np.ix_(minus,plus)])<1e-14
        number=np.block([[residual[np.ix_(plus,plus)],np.zeros((16,16))],
                         [np.zeros((16,16)),-residual[np.ix_(minus,minus)].T]])
        assert abs(np.trace(number))<1e-14 and np.linalg.norm(number,2)<=1.5+1e-13
    assert len(m.source_phase_control()['actual_axis_identities'])==3


def test_complete_root_of_unity_analyzer_has_every_internal_and_spatial_row():
    W,labels=m.analyzer();base,_=m.pairs.coefficient_transform('walsh')
    assert W*W.H==s.eye(32) and len(set(labels))==32
    assert W!=base
    scalar_rows=[j for j,(_,a) in enumerate(labels) if a==0]
    scalar=W[scalar_rows,:]
    assert scalar.H*scalar==s.diag(*(int(a==0) for x in m.CELLS for a in range(4)))
    for j,(k,a) in enumerate(labels):
        for i,x in enumerate(m.CELLS):
            for b in range(4):
                expected=((-1)**sum(u*v for u,v in zip(k,x))*s.exp(-s.I*s.pi*sum(x)/16)/s.sqrt(8)) if a==b else 0
                assert W[j,4*i+b]==expected


def test_full_fock_pair_creation_and_trace_zero_number_bounds():
    # All64 occupation states include blocked and fully occupied pair sectors.
    cs=[car(6,j) for j in range(6)];D=sum((c.conj().T@c for c in cs),np.zeros((64,64),complex))
    coefficients=np.array([[1,.3j,0],[.2,.7,-.4j],[0,.1,1.2]],complex)
    L=sum((coefficients[i,j]*cs[i].conj().T@cs[3+j].conj().T for i in range(3) for j in range(3)),np.zeros((64,64),complex))
    bound=np.linalg.norm(coefficients)**2*(np.eye(64)+D/2)-L.conj().T@L
    assert np.linalg.eigvalsh(bound).min()>-1e-12
    n=np.zeros((6,6),complex)
    for i in range(6):n[i,(i+1)%6]=(.2+.1j)*(i+1)
    N=sum((n[i,j]*cs[i].conj().T@cs[j] for i in range(6) for j in range(6)),np.zeros((64,64),complex))
    for occupation in range(7):
        indices=[j for j in range(64) if j.bit_count()==occupation]
        sector=N[np.ix_(indices,indices)]
        assert np.linalg.norm(sector,2)<=min(occupation,6-occupation)*np.linalg.norm(n,2)+1e-12
    assert all(min(q,32-q)**2<=16*q for q in range(33))


def test_centered_spatial_coefficients_have_actual_cube_covariance():
    q=np.array([.17,.43,1.2]);centered=[np.array(x)-.5 for x in m.CELLS]
    for perm in permutations(range(3)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        for signs in product((-1,1),repeat=3):
            if sign*np.prod(signs)!=1:continue
            R=np.zeros((3,3))
            for i in range(3):R[i,perm[i]]=signs[i]
            for x in centered:
                assert abs(np.exp(-1j*q@x)-np.exp(-1j*(R@q)@(R@x)))<1e-14
    for x in centered:assert abs(np.exp(-1j*q@x)-np.exp(-1j*(-q)@(-x)))<1e-14


def test_bond_control_retains_the_full_physical_Gauss_sector_and_parity():
    c=m.physical_bond_control();V=c['physical_embedding']
    assert V.shape==(48,6) and V.H*V==s.eye(6)
    assert c['Gauss_x']*V==c['Gauss_y']*V==s.zeros(48,6)
    assert (s.eye(48)-V*V.H)*c['full_H']*V==s.zeros(48,6)
    assert c['H']*c['rho']==c['rho']*c['H'] and s.trace(c['rho'])==1
    assert c['parity'].H*c['H']*c['parity']==c['H']
    assert c['parity'].H*c['probes'][0]*c['parity']==-c['probes'][1]


def test_complete_energy_histories_keep_the_original_quantum_reference():
    c=m.physical_bond_control();H=np.array(c['H'],complex);rho=np.array(c['rho'],complex)
    root=np.array(c['root'],complex);purification=root.reshape(-1)
    _,U=np.linalg.eigh(H);P=[np.outer(U[:,i],U[:,i].conj()) for i in range(6)]
    operators=m.bond_instrument(np.pi/16)
    assert np.linalg.norm(sum(K.conj().T@K for K in operators)-np.eye(6))<1e-14
    total=0;density=np.zeros((6,6),complex)
    for initial in P:
        for K in operators:
            for final in P:
                L=final@K@initial;v=np.kron(L,np.eye(6))@purification
                probability=np.vdot(v,v).real;total+=probability
                assert abs(probability-np.trace(L@rho@L.conj().T).real)<1e-14
                amplitude=v.reshape(6,6)
                density+=amplitude.T@amplitude.conj()
    assert abs(total-1)<1e-12
    assert np.linalg.norm(density-rho.T)<1e-12


def test_analyzer_is_not_a_conserved_PVM_and_coherent_choices_matter():
    r=m.control_receipt()
    assert r['probability_normalization_verified'] and r['inversion_energy_covariance_verified']
    assert F(r['numerical_probability_tolerance'])==F(1,10**12)
    assert r['success_effect_commutator_H_norm']>.8
    assert r['success_effect_idempotence_defect']>.3
    assert r['coherent_incoherent_probability_difference']>3e-5


def test_all_scan_bins_have_positive_source_bounds_and_the_union_is_not_renormalized():
    rows=m.phase_scan_rows(F(46,25))
    expected={0:F(7057,1317120),1:F(30027,3920000),2:F(80081,8467200),3:F(6883,705600)}
    for row in rows:
        assert row['original_transfer_weight_lower']==expected[row['shifted_axis_count']]
        assert all((a-b)%2==0 for a,b in zip(row['wavevector_over_pi'],row['principal_wavevector_over_pi']))
    union=sum(row['original_transfer_weight_lower'] for row in rows)
    assert union==F(1367869,20580000) and union-F(1,8192)>F(1,16)
    assert union-F(1,128)>F(15,256)
    assert F(177,6400)-m.centered_moment(F(46,25),F(3,64),7)/9==F(391,44800)
    assert m.analyzer_error(F(1,1024))==F(1,1024)
    with pytest.raises(ValueError):m.analyzer_error(-1)
    with pytest.raises(ValueError):m.centered_moment(33,3,1)


@pytest.fixture(scope='module')
def selected():
    member=m.contract.member(0)
    return member,m.bind_member(member)


def test_same_source_complete_spatial_and_energy_records_have_actual_positive_bounds(selected):
    member,value=selected;r=value['result']
    assert r['complete_probe_outcomes']==64 and r['phase_gate_order']==32
    assert len(r['scalar_spatial_rows'])==8
    assert F(r['selected_scan_union_reported_original_weight_lower'])>F(1,16)
    assert F(r['union_final_excitation_original_weight_lower'])>F(15,256)
    assert F(r['union_final_excitation_tail_upper'])==F(1,128)
    for row in r['scalar_spatial_rows']:
        assert F(row['reported_joint_original_weight_lower'])>0
        assert F(row['final_excitation']['original_final_excitation_weight_lower'])>0
    assert not r['selected_preparation_changed'] and not r['source_clock_or_energy_zero_changed']
    assert r['all_failures_and_quantum_reference_retained']
    assert not r['spatial_analyzer_is_conserved_momentum_PVM']
    assert not r['populated_relative_GNS_momentum_component_proved']
    assert not r['wavevector_zero_establishes_a_rest_frame'] and not r['refinement_survival_proved']
    assert r['physical_mass_table_GeV']==[]
    m.verify_bound(member,value)


def test_source_result_and_phase_error_changes_require_recomputation(selected):
    member,value=selected
    wrong=deepcopy(value);wrong['result']['physical_mass_table_GeV']=[1]
    wrong['result_sha256']=m.contract.digest(wrong['result'])
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
    with pytest.raises(ValueError):m.verify_bound(member,value,phase_diagonal_error='1/1024')
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
    noisy=m.bind_member(member,'1/1024')['result']
    assert F(noisy['selected_scan_union_reported_original_weight_lower'])==F(value['result']['selected_scan_union_reported_original_weight_lower'])-F(1,1024)
    for old,row in zip(value['result']['scalar_spatial_rows'],noisy['scalar_spatial_rows']):
        assert F(row['reported_joint_original_weight_lower'])==F(old['reported_joint_original_weight_lower'])-F(1,1024)
