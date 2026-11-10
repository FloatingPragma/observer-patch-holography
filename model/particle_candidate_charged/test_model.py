"""Independent exact coefficient and full-occupation charged apparatus tests."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import model as m


def creator(j,n):
    out=np.zeros((1<<n,1<<n),complex)
    for bits in range(1<<n):
        if not bits>>j&1:out[bits|1<<j,bits]=(-1)**((bits&((1<<j)-1)).bit_count())
    return out


def test_exact_native_coefficient_ranks_and_absent_flavor_factor():
    r=m.coefficient_control()
    assert [x['exact_constraint_rank'] for x in r['native_blocks']]==[8,3,35,8,0]
    assert r['gauge_and_each_Wilson_coefficient_commutant_dimension']==5
    assert r['generated_local_complex_algebra_dimension']==944
    T=m.native.native_currents();X=s.zeros(15);X[0,11]=1
    assert X*T[-1]!=T[-1]*X
    # A declared additional family factor does have non-scalar commuting coefficients.
    family=s.zeros(3);family[0,1]=1
    A=s.kronecker_product(s.eye(15),family)
    assert all(A*s.kronecker_product(t,s.eye(3))==s.kronecker_product(t,s.eye(3))*A for t in T)
    assert A!=s.trace(A)*s.eye(45)/45
    assert not r['full_H_commutant_computed'] and not r['composite_or_spatial_generations_excluded']


def test_literal_clifford_mass_space_and_parity_rejection():
    r=m.mass_control();g=m.density.gamma_data();b=g[3];g5=m.clean(g[0]*g[1]*g[2]*b)
    alpha=[s.I*b*x for x in g[:3]];odd=s.I*b*g5
    assert all(m.clean(b*x+x*b)==s.zeros(4) for x in alpha)
    assert all(m.clean(odd*x+x*odd)==s.zeros(4) for x in alpha)
    assert m.clean(b*odd*b)==-odd and m.clean(b*b*b)==b
    assert r['unrestricted_mass_space_dimension']==2 and r['parity_even_mass_space_dimension']==1
    assert not r['onsite_coefficient_is_a_rest_mass']


def test_three_same_hypercharge_scalar_vectors_are_distinct_native_sectors():
    rows=m.ward_control()['bare_Yukawa_channels']
    deltas=[tuple(x['required_H_dagger_native_increment']) for x in rows]
    assert deltas==[(0,-1,0,0,-1),(0,0,-1,-1,0),(1,0,1,0,0)]
    assert {sum(a*b for a,b in zip(d,m.CHARGES)) for d in deltas}=={-3}
    assert len(set(deltas))==3
    # Exact character phases distinguish the first two despite equal Q6.
    assert deltas[0][1]!=deltas[1][1]


def test_full64_CAR_scalar_and_dressed_effect_spectra_independently():
    C=[creator(j,6) for j in range(6)];a=[x.conj().T for x in C]
    phi=[(a[2*j]@a[5]-a[2*j+1]@a[4])/np.sqrt(2) for j in range(2)]
    D=[sum(C[2*j+t]@phi[j] for j in range(2)) for t in range(2)]
    actual=m.dressing_control()
    for i in range(2):
        assert np.max(np.abs(phi[i]-np.array(actual['phi'][i],complex)))<1e-14
        assert np.max(np.abs(D[i]-np.array(actual['D'][i],complex)))<1e-14
        assert abs(np.linalg.eigvalsh(phi[i].conj().T@phi[i])[-1]-1)<1e-13
        assert abs(np.linalg.eigvalsh(D[i].conj().T@D[i])[-1]-3)<1e-13
    E=sum(x.conj().T@x for x in D)
    assert abs(np.linalg.eigvalsh(E)[-1]-4.5)<1e-13
    # Anchor normalization cannot replace the full-Fock norm.
    assert np.linalg.eigvalsh(np.eye(64)-E/4)[0]<-.12


def test_dressed_anchor_overlap_is_not_a_full_operator_identity():
    r=m.dressing_control();C=r['creators'];a=[x.H for x in C]
    eps=s.Matrix([[0,1],[-1,0]])
    for j in range(2):
        bare=s.sqrt(2)*sum((eps[j,k]*a[4+k] for k in range(2)),s.zeros(64))
        assert r['D'][j]*r['anchor']==bare*r['anchor']
        empty_weak=s.eye(64)[:,(1<<4)+(1<<5)]
        assert r['D'][j]*empty_weak==s.zeros(64,1)
        assert bare*empty_weak!=s.zeros(64,1)


def test_full256_graded_compensated_menu_normalization_and_anchor():
    C=[creator(j,8) for j in range(8)];a=[x.conj().T for x in C]
    phi=[(a[2*j]@a[5]-a[2*j+1]@a[4])/np.sqrt(2) for j in range(2)]
    D=[sum(C[2*j+t]@phi[j] for j in range(2)) for t in range(2)]
    rows=[C[p]@D[t]/3 for p in (6,7) for t in range(2)]
    E=sum(x.conj().T@x for x in rows);eigs,U=np.linalg.eigh(E)
    assert eigs[0]>-1e-13 and abs(eigs[-1]-1)<1e-13
    failure=(U*np.sqrt(np.maximum(0,1-eigs)))@U.conj().T
    assert np.linalg.norm(failure.conj().T@failure+E-np.eye(256))<1e-12
    actual=m.compensated_control()
    for x,y in zip(rows,actual['rows']):assert np.max(np.abs(x-np.array(y,complex)))<1e-14
    anchor=np.eye(256)[:,63];targets=np.column_stack([x@anchor for x in rows])
    assert np.linalg.norm(targets.conj().T@targets-(2/9)*np.eye(4))<1e-13
    parity=np.diag([(-1)**i.bit_count() for i in range(256)])
    assert all(np.linalg.norm(parity@x-x@parity)==0 for x in rows)


def test_compensating_link_closes_each_endpoint_and_bare_D_is_charged():
    r=m.compensated_control();Qx=6*r['Np'];Qy=r['Qy']
    for B in r['rows']:
        assert m.clean(Qx*B-B*Qx)==6*B
        assert m.clean(Qy*B-B*Qy)==-6*B
        assert m.clean((Qx+Qy)*B-B*(Qx+Qy))==s.zeros(256)
    L=s.zeros(3);L[2,1]=1;E=s.diag(-1,0,1)
    assert E*L-L*E==L
    # The nonzero charge commutator is precisely what the bare regional factor lacks.
    d=m.dressing_control();Q=-3*d['NL']+6*d['NE']
    assert all(m.clean(Q*D-D*Q)==-6*D for D in d['D'])


def test_scalar_Ward_zero_does_not_erase_compensated_or_reference_records():
    # Independent finite two-mode model with exact integer counts and stationary diagonal density.
    C=[creator(j,2) for j in range(2)];N=[x@x.conj().T for x in C]
    rho=np.diag([1,2,3,4])/10
    assert np.trace(rho@C[0].conj().T@C[1])==0
    assert np.trace(rho@C[0]@C[1])==0
    compensated=C[0]@C[1].conj().T
    assert np.trace(rho@compensated.conj().T@compensated).real>0
    assert np.linalg.norm((N[0]+N[1])@compensated-compensated@(N[0]+N[1]))==0
    bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    source_raise=np.kron(np.array([[0,0],[1,0]]),np.eye(2))
    X=(source_raise@np.outer(bell,bell.conj())).reshape(2,2,2,2)
    reference=X[0,:,0,:]+X[1,:,1,:]
    assert abs(np.trace(reference))<1e-14 and np.linalg.norm(reference)>0


def test_exact_positive_mixed_scalar_identity_on_all64_occupations():
    r=m.dressing_control();h=(r['creators'][5].H,-r['creators'][4].H)
    overlap=m.clean(sum((h[j].H*r['D'][j] for j in range(2)),s.zeros(64)))
    effect=m.clean(sum((x.H*x for x in r['phi']),s.zeros(64)))
    assert overlap==m.clean(s.sqrt(2)*effect)
    assert m.clean(r['anchor'].H*overlap*r['anchor'])[0]==2*s.sqrt(2)
    # A charged intermediate Gauss projection is not the product in this identity.
    Q=-3*r['NL']+6*r['NE'];P=s.diag(*(int(Q[j,j]==0) for j in range(64)))
    assert all(P*D*P==s.zeros(64) for D in r['D'])
    assert P*overlap*P!=s.zeros(64)


def test_complete_nine_record_original_overlap_and_phase_countercontrol():
    r=m.coherent_overlap_control();comp=m.compensated_control()
    plus=[np.array(x,complex) for x in r['plus_rows']];minus=[np.array(x,complex) for x in r['minus_rows']]
    E=sum(x.conj().T@x for x in plus+minus);v,U=np.linalg.eigh(E)
    assert v[0]>-1e-13 and v[-1]<=1+1e-13
    failure=(U*np.sqrt(np.maximum(0,1-v)))@U.conj().T
    assert np.linalg.norm(E+failure.conj().T@failure-np.eye(256))<1e-12
    difference=sum(x.conj().T@x for x in plus)-sum(x.conj().T@x for x in minus)
    assert np.linalg.eigvalsh(difference)[0]>-1e-13
    anchor=np.array(comp['anchor'],complex).ravel()
    assert abs(np.vdot(anchor,difference@anchor)-2*np.sqrt(2)/3)<1e-13
    # Independent mixed finite density, not the primary state: the positive-effect bound is state universal.
    rho=(np.eye(256)+np.outer(anchor,anchor.conj()))/257
    probs=[np.trace(rho@x.conj().T@x).real for x in plus+minus+[failure]]
    assert abs(sum(probs)-1)<1e-13
    assert sum(probs[:4])-sum(probs[4:8])>=2*np.sqrt(2)/(3*257)-1e-14
    A=[np.array(x,complex) for x in r['bare_rows']];B=[np.array(x,complex) for x in comp['rows']]
    imaginary=sum(a.conj().T@(1j*b) for a,b in zip(A,B))
    assert abs(np.trace(rho@imaginary).real)<1e-14


def test_optimal_weak_hole_residual_on_all_occupation_blocks():
    C=[creator(j,6) for j in range(6)];a=[x.conj().T for x in C]
    phi=[(a[2*j]@a[5]-a[2*j+1]@a[4])/np.sqrt(2) for j in range(2)]
    D=[sum(C[2*j+t]@phi[j] for j in range(2)) for t in range(2)]
    h=[a[5],-a[4]];R=[D[j]-np.sqrt(2)*h[j] for j in range(2)]
    G=sum(x.conj().T@x for x in R);maxima=[]
    for holes in range(5):
        indices=[b for b in range(64) if 4-(b&15).bit_count()==holes]
        block=G[np.ix_(indices,indices)];maxima.append(np.linalg.eigvalsh(block)[-1])
    assert np.linalg.norm(np.array(maxima)-[0,1,3,3,4])<1e-13
    deficit=np.diag([4-(b&15).bit_count() for b in range(64)])
    assert np.linalg.eigvalsh(1.5*deficit-G)[0]>-1e-13
    assert np.linalg.eigvalsh(1.49*deficit-G)[0]<-.019
    exact=m.weak_hole_residual_control()
    assert [x['largest_eigenvalue'] for x in exact['blocks']]==['0','1','3','3','4']


def test_success_history_norm_with_same_unsharp_clocks_and_different_failures():
    comp=m.compensated_control();overlap=m.coherent_overlap_control()
    A=[np.array(x,complex) for x in overlap['bare_rows']];B=[np.array(x,complex) for x in comp['rows']]
    r=2*np.sqrt(2)/3;R=[b-r*a for a,b in zip(A,B)]
    # Independent finite control. Clocks are nonprojective functions of one diagonal H.
    energies=np.linspace(-2,3,256);rho=(np.arange(256)+1).astype(float);rho/=rho.sum()
    root=np.diag(np.sqrt(rho));clock=[np.diag(np.cos(energies)),np.diag(np.sin(energies))]
    initial=[f@root for f in clock]
    residual=sum(np.linalg.norm(f@x@v,'fro')**2 for f in clock for x in R for v in initial)
    direct=sum(np.trace(np.diag(rho)@x.conj().T@x).real for x in R)
    assert abs(residual-direct)<1e-13
    # The same retained final bin and subset of success labels defines both events.
    pa=sum(np.linalg.norm(clock[1]@x@v,'fro')**2 for x in A[:2] for v in initial)
    pb=sum(np.linalg.norm(clock[1]@x@v,'fro')**2 for x in B[:2] for v in initial)
    assert np.sqrt(pb)+1e-13>=r*np.sqrt(pa)-np.sqrt(residual)
    EA=sum(x.conj().T@x for x in A);EB=sum(x.conj().T@x for x in B)
    assert np.linalg.norm(EA-EB)>0.1


def test_actual_members_preserve_identity_scope_and_tiny_original_floor():
    for n in (0,1,2,8):
        member=m.contract.member(n);bound=m.bind_member(member);r=bound['result']
        assert all(r[k]==v for k,v in member['ids'].items())
        assert r['physical_mass_table_GeV']==[] and r['physical_species_confirmed']==[]
        assert F(r['dressed_success_original_weight_lower_factored']['prefactor'])==F(8,9)
        assert not r['dressed_positive_narrow_energy_band_derived']
        assert not r['source_state_reference_clock_changed']
        if n==0:
            assert r['primary_member_zero_regional_energy_attachment']['raw_spin_per_joint_row_original_energy_weight_simple_lower']=='1/16'
            w=r['original_compensated_dressed_energy_window']
            assert F(w['original_dressed_recorded_weight_strict_lower'])==F(109561,4000000)>F(1,40)
            assert w['reported_transfer_interval']==['77/8','271/8']
            assert w['additional_clock_loss_subtracted']==0 and w['complete_probe_outcomes']==97
            assert r['dressed_positive_reported_energy_window_derived']
        else:
            assert r['primary_member_zero_regional_energy_attachment'] is None
            assert r['original_compensated_dressed_energy_window'] is None
    member=m.contract.member(0);value=m.bind_member(member);bad=deepcopy(value)
    bad['result']['candidate_discriminators'][0]['mass_GeV']=1
    bad['result_sha256']=m.contract.digest(bad['result'])
    with pytest.raises(ValueError):m.verify_bound(member,bad)
