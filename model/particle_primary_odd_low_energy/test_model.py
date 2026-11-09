"""Independent complete-filter, clock, source-bound and depth controls."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from . import model as m


def _toy():
    rng=np.random.default_rng(9041);H=np.zeros((16,16),complex)
    for number in range(5):
        rows=[x for x in range(16) if x.bit_count()==number]
        X=rng.normal(size=(len(rows),len(rows)))
        H[np.ix_(rows,rows)]=(X+X.T)/5+number*np.eye(len(rows))
    H-=np.linalg.eigvalsh(H)[0]*np.eye(16)
    energies,vec=np.linalg.eigh(H)
    rho=vec@np.diag((1+energies[-1]-energies)**2)@vec.conj().T;rho/=np.trace(rho)
    A=np.array(m.odd.creation_word((0,1,2),4),complex)
    D=np.eye(16)-A.conj().T@A
    return H,energies,vec,rho,A,D


def _selector(H,A,D,*,cayley=False,carrier=.7):
    grid=m.apparatus.quadrature(F(1,3),6)
    weights=np.array([float(r['weight']) for r in grid['rows']]);v=np.sqrt(weights)
    e=np.zeros(len(v));e[0]=1;u=e-v
    R=np.eye(len(v))-2*np.outer(u,u)/(u@u)
    operators=[[],[]]
    for r in range(len(v)):
        out=[np.zeros_like(A,dtype=complex),np.zeros_like(A,dtype=complex)]
        for j,row in enumerate(grid['rows']):
            t=float(row['time'])
            if cayley:
                steps=3;I=np.eye(len(H))
                C=np.linalg.matrix_power(np.linalg.solve(I+1j*t*H/(2*steps),I-1j*t*H/(2*steps)),steps)
                z=((1-1j*t*carrier/(2*steps))/(1+1j*t*carrier/(2*steps)))**steps
            else:C=expm(-1j*t*H);z=np.exp(-1j*carrier*t)
            for p,M in enumerate((A,D)):out[p]+=R[j,r]*v[j]*z*C.conj().T@M@C
        for p in (0,1):operators[p].append(out[p])
    return operators


@pytest.mark.parametrize('cayley',[False,True])
def test_every_selector_preserves_each_transition_and_complete_reference(cayley):
    H,E,V,rho,A,D=_toy();ops=_selector(H,A,D,cayley=cayley)
    complete=sum(K.conj().T@K for group in ops for K in group)
    assert np.max(abs(complete-np.eye(16)))<2e-13
    for group,M in zip(ops,(A,D)):
        resolved=sum(abs(V.conj().T@K@V)**2 for K in group)
        assert np.max(abs(resolved-abs(V.conj().T@M@V)**2))<2e-13
    root=V@np.diag(np.sqrt(np.diag(V.conj().T@rho@V).real))@V.conj().T
    branches=[K@root for group in ops for K in group]
    assert abs(sum(np.vdot(x,x).real for x in branches)-1)<1e-13
    reference=sum(x.T@x.conj() for x in branches)
    assert np.max(abs(reference-rho.T))<2e-13
    original=np.diag(V.conj().T@A@rho@A.conj().T@V).real
    filtered=sum(np.diag(V.conj().T@K@rho@K.conj().T@V).real for K in ops[0])
    assert np.max(abs(original-filtered))<1e-13


def test_stationarity_and_coherent_selector_boundaries_are_real():
    H,E,V,rho,A,D=_toy();ops=_selector(H,A,D)
    rng=np.random.default_rng(294);psi=rng.normal(size=16)+1j*rng.normal(size=16);psi/=np.linalg.norm(psi)
    nonstationary=np.outer(psi,psi.conj())
    original=np.diag(V.conj().T@A@nonstationary@A.conj().T@V).real
    after=sum(np.diag(V.conj().T@K@nonstationary@K.conj().T@V).real for K in ops[0])
    assert np.max(abs(original-after))>1e-6
    # Keeping only the blank changes the law even for the stationary input.
    blank=np.diag(V.conj().T@ops[0][0]@rho@ops[0][0].conj().T@V).real
    stationary=np.diag(V.conj().T@A@rho@A.conj().T@V).real
    assert np.max(abs(blank-stationary))>1e-4


def test_carrier_modulation_selects_frequency_without_changing_energy():
    H=np.diag([0.,9.]);A=np.array([[0.,0.],[1.,0.]])
    D=np.diag([0.,1.]);ops=_selector(H,A,D,carrier=9.)
    assert np.max(abs(ops[0][0]-A))<1e-13
    state=np.array([1.,0.]);out=ops[0][0]@state
    assert np.vdot(out,H@out).real==pytest.approx(9.)
    assert np.linalg.eigvalsh(H).tolist()==[0.,9.]


def test_arbitrary_static_number_phase_only_changes_each_local_amplitude():
    A=np.array(m.odd.creation_word((0,1,2),6),complex)
    B=np.array(m.odd.creation_word((3,4,5),6),complex)
    theta=(.37,-.21)
    phase=np.array([theta[0]*(mask&7).bit_count()+theta[1]*(mask>>3).bit_count() for mask in range(64)])
    Z=np.diag(np.exp(1j*phase))
    for C,t in zip((A,B),theta):
        assert np.max(abs(Z@C@Z.conj().T-np.exp(3j*t)*C))<1e-13
        failure=np.eye(64)-C.conj().T@C
        assert np.max(abs(Z@failure@Z.conj().T-failure))<1e-13
    # Relative phases change a coherently combined spatial operator.
    coherent=(A+B)/np.sqrt(2)
    transformed=Z@coherent@Z.conj().T
    assert np.linalg.norm(transformed-np.exp(3j*theta[0])*coherent,2)>.5


def test_selected_moment_allows_both_present_and_absent_fixed_energy_weight():
    # These comparison measures diagnose the moment's information; neither
    # is substituted for the actual selected source state.
    for n in (1,4,8):
        a=F(1,2**n);center=9/a
        for low in (F(),F(1)):
            second=F(1,9)*(center-low)**2
            assert second<=9/a**2
            assert F(8,9)>F(1,50)
        # Entire mass at the centered frequency satisfies the same cap.
        assert F()==(center-center)**2
    bound=m.low_pass_bound(m.contract.member(8))
    assert F(bound['ideal_zero_centered_blank_success_upper'])<F(113,1000)
    assert bound['fixed_T_lim_sup_blank_success_upper']=='1/9'
    assert m.low_pass_bound(m.contract.member(0))['ideal_zero_centered_blank_success_upper']=='97/400'
    for n in (0,1,2,8):
        assert F(m.low_pass_bound(m.contract.member(n))['ideal_zero_centered_blank_success_upper'])<F(1,4)


def test_projector_purification_tail_bound_uses_original_probability():
    H=np.diag([0.,.2,2.,7.]);Fmat=np.array([[0.,.1,.03,.02],[.1,0,.1,.04],[.03,.1,0,.2],[.02,.04,.2,0.]])
    double=Fmat@(H@Fmat-Fmat@H)-(H@Fmat-Fmat@H)@Fmat
    drift=np.sqrt(np.linalg.norm(double,2)/2);U=expm(1j*Fmat)
    sigma=np.diag([0.,0.,.1,.3]);p=.4;L=2.;R=.5
    P=np.diag([1.,1.,0.,0.]);Q=np.eye(4)-P
    assert np.linalg.norm(Q@U.conj().T@P,2)<=(np.sqrt(R)+drift)/np.sqrt(L)+1e-13
    tail=np.trace(U@sigma@U.conj().T@Q).real
    assert tail>=max(0,np.sqrt(p)-(np.sqrt(R)+drift)/np.sqrt(L))**2-1e-13
    assert np.trace(sigma)==pytest.approx(.4)


@pytest.fixture(scope='module')
def receipts():
    return {n:m.bind_member(m.contract.member(n)) for n in (0,1,2,8)}


def test_actual_complete_filter_weights_and_clock_budgets(receipts):
    for n,bound in receipts.items():
        r=bound['result'];f=r['populated_centered_filter'];a=F(1,2**n)
        assert f['selector_dimension']==1024
        assert F(f['half_width'])==a/6 and F(f['analyzer_carrier_frequency'])==9/a
        assert f['quadrature_operator_error_upper']=='1/2048'
        assert F(f['exact_node_instrument_operator_error_upper'])<=F(1,16384)
        assert F(f['blank_operator_error_to_continuous_triangle_upper'])<=F(9,16384)
        assert F(f['complete_energy_clock_event_error_upper'])<=F(1,4096)
        assert F(f['reported_original_blank_band_weight_lower'])>=F(1,80)
        assert r['complete_selector_summed_creation_energy_measure_exactly_original']
        assert not f['analyzer_carrier_subtracted_from_source_energy']
        assert not f['division_by_success_probability']
        assert r['physical_mass_table_GeV']==[]
    sharp=receipts[0]['result']['populated_centered_filter']['member_zero_sharper_reported_band']
    assert F(sharp['reported_original_blank_success_weight_lower'])>F(3,5)
    assert list(map(F,sharp['reported_transfer_band']))==[F(11,4),F(61,4)]


def test_required_depth_and_invalid_domains():
    assert m.dressing_depth(20)['critical_Weyl_word_length_necessary_lower']==3
    assert m.dressing_depth(40)['critical_Weyl_word_length_necessary_lower']==2331
    for level in (True,-1,1.5):
        with pytest.raises(ValueError):m.dressing_depth(level)
    for target in ('1/50','-1',.001):
        with pytest.raises(ValueError):m.dressing_depth(8,tail_target=target)
    with pytest.raises(ValueError):m.low_pass_bound(m.contract.member(0),'0')


def test_source_and_forged_low_energy_claims_are_rejected(receipts):
    member=m.contract.member(0);m.verify_bound(member,receipts[0])
    wrong=deepcopy(receipts[0]);wrong['result']['fixed_positive_energy_odd_residue_proved']=True
    wrong['result_sha256']=m.contract.digest(wrong['result'])
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
    wrong=deepcopy(member);wrong['clock']['Schrodinger']='exp(-2itH)'
    with pytest.raises(ValueError):m.bind_member(wrong)
