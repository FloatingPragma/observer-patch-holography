"""Exact graph/Haar checks and an independent charge-six gauge compression."""
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_actual_cube_Gauss_support_and_private_edges():
    g=m.cube_support();counts=g['supports_with_no_degree_one_vertex_counts']
    assert counts[0]==1 and counts[4]==6 and counts[6]==16
    assert all(n not in counts for n in (1,2,3,5))
    assert g['private_edge_patterns_checked_by_total_degree']=={2:15,3:50,4:105}
    assert len(g['interior_edges'])==12
    assert len({f['edges'] for f in g['oriented_outward_faces']})==6


def test_actual_degree_one_native_Haar_irreps_and_moments():
    data=m.native_Haar_control()
    assert data['native_character_orthogonality_pairs_checked']==100
    assert data['actual_degree_one_link_dimension']==119
    assert data['one_link_electric_energies']==(F(42),F(189,8),F(441,8),F(147,4),F(63,4))
    assert data['first_electric_Gauss_level_multiplicity']==12
    assert data['next_electric_Gauss_level_lower']==F(189,2)
    assert data['ground_anchor_residual_squared']==F(1,735)
    assert data['odd_loop_residual_squared_upper']==F(59,44100)
    assert not data['odd_residual_bound_equals_compressed_second_moment_claimed']


def test_exact_cube_inversion_separates_spatial_and_fermion_parity():
    c=m.cube_channels()
    assert c['first_cluster_decomposition']=={
        'real':{'A1g':1,'Eg':1,'T1u':1},
        'imaginary':{'A1u':1,'Eu':1,'T1g':1}}
    assert c['spatial_odd_A1u_multiplicity_in_cluster']==1
    assert c['A1u_CAR_parity']==0
    assert not c['A1u_is_Lorentz_or_fermion_spin']


@lru_cache(maxsize=1)
def q6_compression():
    """Actual q6-link subspace control; other native links are deliberately absent."""
    graph=m.cube_support();edges=graph['interior_edges'];cells=graph['cells']
    incidence=s.zeros(7,12)
    for i,(x,y) in enumerate(edges):
        if x<7:incidence[x,i]=-1
        if y<7:incidence[y,i]=1
    tree=[];parent=list(range(8))
    def find(x):
        while parent[x]!=x:x=parent[x]
        return x
    for i,(x,y) in enumerate(edges):
        a,b=find(x),find(y)
        if a!=b:parent[a]=b;tree.append(i)
    chords=[i for i in range(12) if i not in tree]
    transform=-incidence[:,tree].inv()*incidence[:,chords]
    states=[]
    for free in product((-1,0,1),repeat=5):
        assigned=transform*s.Matrix(free)
        if any(v not in (-1,0,1) for v in assigned):continue
        q=[0]*12
        for i,v in zip(tree,assigned):q[i]=int(v)
        for i,v in zip(chords,free):q[i]=v
        assert incidence*s.Matrix(q)==s.zeros(7,1)
        states.append(tuple(q))
    where={q:i for i,q in enumerate(states)}
    fluxes=[]
    for face in graph['oriented_outward_faces']:
        v0=min(face['vertices']);neighbors=[y if x==v0 else x for x,y in (edges[i] for i in face['edges']) if v0 in (x,y)]
        v1,v3=sorted(neighbors);v2=next(v for v in face['vertices'] if v not in (v0,v1,v3))
        u=np.array(cells[v1])-cells[v0];v=np.array(cells[v3])-cells[v0]
        outward=np.zeros(3,dtype=int);outward[face['axis']]=2*face['side']-1
        tour=(v0,v1,v2,v3)
        if np.cross(u,v)@outward<0:tour=(v0,v3,v2,v1)
        q=[0]*12
        for a,b in zip(tour,tour[1:]+tour[:1]):
            e=(min(a,b),max(a,b));q[edges.index(e)]=1 if a<b else -1
        assert tuple(q) in where
        fluxes.append(tuple(q))
    n=len(states);T=np.zeros((n,n));A=np.zeros((n,n),complex)
    for j,q in enumerate(states):
        for face in fluxes:
            for sign in (-1,1):
                target=tuple(a+sign*b for a,b in zip(q,face))
                if target in where:
                    T[where[target],j]-=1/210
                    A[where[target],j]+=sign/(12j)
    electric=np.array([float(F(63,4))*sum(x*x for x in q) for q in states])
    H=np.diag(electric)+T
    omega=np.zeros(n);omega[where[(0,)*12]]=1
    phi=np.zeros(n,complex)
    for face in fluxes:
        phi[where[face]]+=1/(1j*np.sqrt(12))
        phi[where[tuple(-x for x in face)]]-=1/(1j*np.sqrt(12))
    return H,A,omega,phi,electric


def test_independent_exact_Gauss_flow_carrier_and_actual_native_omission():
    H,A,omega,phi,electric=q6_compression()
    assert len(H)==69
    assert np.count_nonzero(electric==0)==1
    assert np.count_nonzero(electric==63)==12
    assert min(electric[electric>63])==94.5
    assert np.allclose(A,A.conj().T)
    assert np.linalg.norm(A@omega-phi/np.sqrt(12))<1e-15
    assert abs(np.vdot(phi,phi)-1)<1e-15
    assert np.linalg.norm(H@omega)**2==pytest.approx(float(F(1,3675)))
    # The full source has five native character pairs. The q6 compression
    # cannot replace it: its anchor residual omits four native contributions.
    assert F(1,735)==5*F(1,3675)
    assert abs(np.vdot(phi,(H-63*np.eye(len(H)))@phi))<1e-13


def test_isolated_cluster_and_nonzero_residue_in_independent_compression():
    H,A,omega,phi,electric=q6_compression()
    E,V=np.linalg.eigh(H);cluster=(E>62)&(E<64)
    assert np.count_nonzero(cluster)==12
    assert E[13]>93
    u=V[:,cluster]@(V[:,cluster].conj().T@phi);u/=np.linalg.norm(u)
    Eu=float(np.vdot(u,H@u).real)
    assert np.linalg.norm(H@u-Eu*u)<1e-10
    bounds=m.spectral_bounds();lo,hi=bounds['gauge_ground_to_A1u_energy_gap_interval']
    assert float(lo)<Eu-E[0]<float(hi)
    assert abs(np.vdot(u,A@V[:,0]))**2>float(F(4,49))
    assert not bounds['raw_loop_probe_has_no_other_transition_frequencies_claimed']


def test_Temple_and_phase_bounds_are_rational_outward_bounds():
    d=m.spectral_bounds()
    assert d['interacting_gauge_ground_energy_interval']==[-F(1,45675),F()]
    assert d['interacting_A1u_energy_interval']==[63-F(59,1351350),F(63)]
    assert d['A1u_gap_to_same_symmetry_complement_lower']==F(429,14)
    assert d['phase_aligned_eigenvector_distance_sum_upper']<F(1,364)
    assert F(15,52)**2<F(1,12)
    assert F(15,52)-F(2,7)==F(1,364)
    assert not d['A1u_energy_is_globally_nondegenerate_across_other_cube_irreps']


@pytest.fixture(scope='module')
def bound():
    member=m.contract.member(0)
    return member,m.bind_member(member)


def test_actual_original_atom_and_complete_clock_weight_keep_extreme_rarity(bound):
    member,value=bound;r=value['result']
    assert r['original_energy_local_population_floor']['negative_exponent']==32406376
    assert F(r['original_resolved_positive_transfer_atom_lower_factored']['rational_multiplier'])==F(4,49)
    assert F(r['original_resolved_retarded_positive_residue_lower_factored']['rational_multiplier'])==F(2,49)
    event=r['complete_clocked_original_event']
    loss=F(event['one_clock_loss_upper']);factor=F(event['original_weight_lower_factored']['rational_multiplier'])
    assert loss<F(1,512)
    assert factor==F(4,49)*(1-loss)**2>F(1,13)
    assert not event['empty_sector_success_divided_out']
    assert not event['joint_reference_Loewner_bound_asserted']
    assert not event['all_other_initial_state_clock_leakage_excluded']
    assert r['extreme_sector_supported_contraction_weight_vanishes_on_selected_refinement']
    assert not r['CAR_odd_trion_branch_isolated']
    assert not r['ambient_primary_vacuum_gap_proved']
    assert r['physical_mass_table_GeV']==[]
    upper=F(r['original_extreme_sector_suppression']['combined_original_sector_probability_upper'])
    assert 0<upper<F(1,10000)


def test_rehashed_source_or_mass_promotions_are_rejected(bound):
    member,value=bound;m.verify_bound(member,value)
    for key,replacement in (('physical_mass_table_GeV',[63]),('CAR_odd_trion_branch_isolated',True)):
        bad=deepcopy(value);bad['result'][key]=replacement;bad['result_sha256']=m.contract.digest(bad['result'])
        with pytest.raises(ValueError):m.verify_bound(member,bad)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
