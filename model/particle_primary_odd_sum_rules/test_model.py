"""Independent CAR, original spectral, complete-clock and parity controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest
import sympy as s
from scipy.sparse import csr_matrix

from . import model as m


def creation(mode,count):
    A=np.zeros((1<<count,1<<count),complex)
    for mask in range(1<<count):
        if not mask>>mode&1:
            A[mask|(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
    return A


def test_positive_original_fermionic_measure_has_correct_reflected_adjoint_moments():
    C=[creation(i,6) for i in range(6)];N=[A@A.conj().T for A in C]
    H=sum((i+1)*N[i]/7 for i in range(6))
    H+=sum((C[i]@C[3+i].conj().T+C[3+i]@C[i].conj().T)/3 for i in range(3))
    H+=N[0]@N[3]/5
    energies,U=np.linalg.eigh(H);p=(10-energies)**6;p/=sum(p);rho=(U*p)@U.conj().T
    fields=(C[0]@C[1]@C[2],C[3]@C[4]@C[5]);omega=energies[:,None]-energies[None,:]
    matrices=[U.conj().T@A@U for A in fields]
    for x,y,r in product(range(2),range(2),range(3)):
        Fx,Fy=fields[x],fields[y];Ax,Ay=matrices[x],matrices[y]
        plus=np.sum(p[None,:]*Ax.conj()*Ay*omega**r)
        # Indices transpose on the annihilation kernel; energy reverses too.
        minus=np.sum(p[None,:]*(Ay.conj().T).conj()*(Ax.conj().T)*(-omega)**r)
        D=Fy.copy()
        for _ in range(r):D=H@D-D@H
        expected=np.trace(rho@(Fx.conj().T@D+D@Fx.conj().T))
        assert abs(plus+minus-expected)<2e-13
    for sign in (-1,1):
        A=(fields[0]+sign*fields[1])/np.sqrt(2)
        actual=A.conj().T@A+A@A.conj().T
        exact=sum(Fx.conj().T@Fx+Fx@Fx.conj().T for Fx in fields)/2
        assert np.linalg.norm(actual-exact)<1e-14
    plus=(fields[0]+fields[1])/np.sqrt(2);minus=(fields[0]-fields[1])/np.sqrt(2)
    assert np.linalg.norm(plus.conj().T@plus-minus.conj().T@minus)>1


def test_full_twelve_mode_literal_hopping_commutator_reproduces_Haar_norms():
    spin=m.coherent.spin;local=spin.local_controls();active=local['active_source_modes'];blocks=dict(m.composite.sectors.native_blocks())
    controls=m.neighbor_graded_control();actual=[]
    C=[spin._creation_sparse((i,),12) for i in range(12)];empty=np.zeros(4096,dtype=np.int64);empty[0]=1
    for label,family in local['families'].items():
        for polynomial in family['polynomials']:
            square,coefficients=spin._integer_polynomial(polynomial)
            Fx=sum((v*spin._creation_sparse(w,12) for w,v in coefficients.items()),csr_matrix((4096,4096),dtype=np.int64))
            Fy=sum((v*spin._creation_sparse(tuple(i+6 for i in w),12) for w,v in coefficients.items()),csr_matrix((4096,4096),dtype=np.int64))
            witness=Fx@empty;norm=F()
            for charge in (-3,6):
                dimension=len(blocks[charge])
                for left,right in product(range(dimension),repeat=2):
                    T=csr_matrix((4096,4096),dtype=np.int64)
                    for spinrow in range(2):
                        a=active.index(15*spinrow+blocks[charge][left]);b=active.index(15*spinrow+blocks[charge][right])
                        T+=C[a]@C[6+b].T
                    D=T@Fy-Fy@T;G=Fx.T@D+D@Fx.T;G.eliminate_zeros()
                    row,col=G.nonzero()
                    assert all((int(i)&63).bit_count()-(int(j)&63).bit_count()==-2 for i,j in zip(row,col))
                    assert all((int(i)>>6).bit_count()-(int(j)>>6).bit_count()==2 for i,j in zip(row,col))
                    output=G@witness;norm+=square**3*int(output@output)/(4*dimension)
            actual.append(norm)
    assert actual==[F(row['first_graded_response_witness_norm_squared']) for row in controls['all_normalized_physical_trion_witnesses']]
    assert actual==[F(3,8),F(3,8),F(3,4),F(1,4),F(1,4),F(3,4)]


def test_distinct_nonneighbor_first_graded_response_vanishes_by_actual_support():
    graph=m.coherent.cube_stars();pairs=list(product(graph['cells'],repeat=2));seen=[0,0]
    for x,y in pairs:
        if x==y:continue
        neighbor=sum(a!=b for a,b in zip(x,y))==1
        if neighbor:seen[0]+=1
        else:
            assert x not in graph['stars'][y]['matter'];seen[1]+=1
    assert seen==[24,32]
    controls=m.sum_rule_control();assert controls['first_summed_shells']==['m0','m1','0','0']
    assert controls['first_Walsh_moments']==['m0+3*m1','m0+m1','m0-m1','m0-3*m1']
    assert not controls['m1_original_value_or_sign_computed']


def test_Cayley_graded_clicks_are_complete_and_show_controlled_spatial_contrast():
    C=[creation(i,3) for i in range(3)];H=sum(C[i]@C[i].conj().T for i in range(3))
    H-=sum((C[i]@C[i+1].conj().T+C[i+1]@C[i].conj().T)/2 for i in range(2))
    R=F(4);N=4;c=2.**-N;t=c/(64*float(R)**2);plan=m.finite_time_plan(R,N)
    U=np.linalg.solve(np.eye(8)+1j*t*H/2,np.eye(8)-1j*t*H/2);rho=np.eye(8)/8
    assert np.linalg.norm(U.conj().T@U-np.eye(8))<1e-14
    probabilities=[]
    for y in (1,2):
        Y=U.conj().T@C[y]@U;K=(C[0].conj().T@Y+Y@C[0].conj().T)/2
        values,V=np.linalg.eigh(np.eye(8)-K.conj().T@K);failure=(V*np.sqrt(np.maximum(values,0)))@V.conj().T
        assert min(values)>0 and np.linalg.norm(K.conj().T@K+failure.conj().T@failure-np.eye(8))<1e-14
        p=np.trace(rho@K.conj().T@K).real;probabilities.append(p)
        assert abs(p+np.trace(rho@failure.conj().T@failure).real-1)<1e-14
    assert probabilities[0]>=c*t*t/64
    assert probabilities[1]<=9*float(R)**4*t**4
    assert probabilities[0]-probabilities[1]>=c*t*t/128
    assert plan['nonneighbor_to_neighbor_bound_ratio_upper']==F(9,64)


def test_native_chirality_has_same_gauge_representation_and_separate_parity():
    row=m.chirality_parity_control();g5=s.Matrix(row['native_gamma5']);plus,minus=map(s.Matrix,row['native_chiral_projectors'])
    beta=m.composite.sectors.branch.response.density.gamma_data()[3]
    assert beta*plus*beta==minus and plus+minus==s.eye(4)
    assert g5.H==g5 and g5**2==s.eye(4)
    assert row['native_gauge_generators_checked']==12 and row['proper_spin_lifts_checked']==48
    assert row['same_internal_representation_in_both_projectors']
    assert not row['onsite_H_conserves_native_gamma5']
    assert not row['native_gamma5_is_already_proved_physical_Lorentz_chirality']
    assert not row['parity_conjugates_internal_gauge_representation']


def test_complete_chiral_and_clock_records_have_original_reflected_laws():
    gamma=[np.array(A,dtype=complex) for A in m.composite.sectors.branch.response.density.gamma_data()];beta=gamma[3]
    hop=(-beta+beta@gamma[0])/2;H=np.block([[3*beta,hop],[hop.conj().T,3*beta]])
    P=np.kron(np.array([[0,1],[1,0]]),beta)
    assert np.linalg.norm(P@H@P.conj().T-H)<1e-14
    raw=10*np.eye(8)-H;rho=raw@raw;rho/=np.trace(rho)
    plus,minus=[np.array(s.Matrix(A),complex) for A in m.chirality_parity_control()['native_chiral_projectors']]
    probes={}
    for sign,G in ((1,plus),(-1,minus)):
        for x,y in ((0,1),(1,0)):
            B=np.zeros((8,8),complex);B[4*x:4*x+4,4*y:4*y+4]=G/4;probes[sign,x,y]=B
    effect=sum(A.conj().T@A for A in probes.values());values,V=np.linalg.eigh(np.eye(8)-effect);failure=(V*np.sqrt(values))@V.conj().T
    assert min(values)>0 and np.linalg.norm(P@failure@P.conj().T-failure)<1e-14
    U=np.linalg.solve(np.eye(8)+1j*H/16,np.eye(8)-1j*H/16);clocks=((np.eye(8)+U)/2,(np.eye(8)-U)/2)
    for key,A in probes.items():
        sign,x,y=key;paired=probes[-sign,1-x,1-y]
        assert np.linalg.norm(P@A@P.conj().T-paired)<1e-14
        for first,last in product(clocks,repeat=2):
            K=last@A@first;L=last@paired@first
            assert abs(np.trace(K@rho@K.conj().T)-np.trace(L@rho@L.conj().T))<1e-14
    assert np.linalg.norm(effect+failure.conj().T@failure-np.eye(8))<1e-14


def test_inversion_on_interior_graph_is_uniform_but_not_an_ambient_ghost_rotation():
    for L in (2,3,4):
        cells=tuple(product(range(L),repeat=3));edges={(x,y) for x,y in product(cells,repeat=2) if x<y and sum(abs(a-b) for a,b in zip(x,y))==1}
        reflect=lambda x:tuple(L-1-a for a in x)
        assert {tuple(sorted((reflect(x),reflect(y)))) for x,y in edges}==edges
    for level in (0,1,2,8):
        row=m.bind_chirality_member(m.contract.member(level))['result']
        assert row['source_link_degree']==level+1
        assert row['complete_success_records_per_native_block']==12*row['source_side']**2*(row['source_side']-1)
        assert row['one_link_ghost_Gauss_constraints_reduce_every_dangling_link_to_constant']
        assert row['native_Clifford_controls']['same_native_gauge_currents_in_both_projectors']
        assert not row['one_handed_physical_low_energy_sector_proved']


def test_trions_are_actual_gauge_singlets_not_gauge_charged_lepton_labels():
    row=m.gauge_charge_control();assert row['trion_total_Q6_increment']==0
    assert row['actual_local_gauge_kernel_dimension']==6
    assert row['all_twelve_local_gauge_generators_annihilate_trion_kernel']
    assert row['known_physical_string_endpoint_Q6_increments']==[-6,6]
    assert row['known_physical_string_probe_CAR_parity']==0
    assert not row['trion_native_increment_is_nonzero_electric_charge']
    assert not row['single_gauge_covariant_creator_preserves_original_physical_carrier']
    assert not row['gauge_singlet_alone_identifies_or_excludes_a_dressed_neutrino']


@pytest.fixture(scope='module')
def actual():return m.bind_member(m.contract.member(0))


def test_actual_original_source_record_bound_uses_tiny_factored_time_and_weights(actual):
    r=actual['result'];plan=r['finite_complete_graded_probe'];N=r['physical_density_floor_factored']['negative_exponent']
    assert N==582407133
    assert plan['time']['negative_exponent']==N
    assert plan['neighbor_original_click_probability_lower']['negative_exponent']==3*N
    assert plan['nonneighbor_original_click_probability_upper']['negative_exponent']==4*N
    assert plan['complete_outcome_count']==2 and plan['all_failures_source_and_reference_retained']
    assert r['nonzero_original_four_point_spatial_contrast_proved']
    assert not r['nonzero_original_odd_two_point_spatial_shell_proved']
    assert not r['isolated_particle_pole_or_relativistic_dispersion_proved']
    for row in r['original_signed_success_channels']:
        assert F(row['each_complete_Walsh_signed_success_lower'])>0
        assert row['same_original_success_probability_for_all_eight_rows']
    assert r['physical_mass_table_GeV']==[]


def test_changed_source_and_physical_promotions_are_rejected(actual):
    m.verify_bound(m.contract.member(0),actual)
    wrong=deepcopy(actual);wrong['result']['nonzero_original_odd_two_point_spatial_shell_proved']=True
    with pytest.raises(ValueError):m.verify_bound(m.contract.member(0),wrong)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
    member=m.contract.member(1);bound=m.bind_chirality_member(member);m.verify_chirality_bound(member,bound)
    bad=deepcopy(bound);bad['result']['one_handed_physical_low_energy_sector_proved']=True
    with pytest.raises(ValueError):m.verify_chirality_bound(member,bad)


def test_canonical_receipt_reuses_exact_source_and_finite_chirality_contract():
    c=m.certificate();assert json.loads(json.dumps(c,sort_keys=True))==c
    assert len(c['members'])==1 and len(c['chirality_members'])==4
    assert 'Lean/Screen/SMStructureAdequacySurface.lean' in c['source_sha256']
    assert 'model/particle_primary_composite_transport/model.py' in c['source_sha256']
