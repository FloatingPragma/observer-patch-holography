"""Independent CAR, coherent-control and reference-history checks."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm

from . import model as m
from model.particle_primary_current_tail.test_model import _lift
from model.particle_local_wilson_density import density


@pytest.mark.parametrize('R',[1,2])
def test_actual_open_tent_gradient_and_normalization(R):
    L=2*R*R;a=F(1,R)
    t=np.array([max(0,1-min(k,L-k)/R) for k in range(L)])
    profile=t[:,None,None]*t[None,:,None]*t[None,None,:]
    assert abs(profile.sum()-R**3)<1e-13
    cap=F(1,10) if R==1 else (2*R*R+1)**2/F(45*R**8)
    for X in product(range(L),repeat=3):
        q=np.roll(profile,X,axis=(0,1,2))/(30*R**3)
        exact=30/float(a)*sum(np.sum(np.diff(q,axis=j)**2) for j in range(3))
        assert exact<=float(cap)+1e-15
    assert cap<=F(9,80)*a**4


def test_full_CAR_hopping_and_balanced_meter_identity():
    g=density.gamma_data();S=np.array((-g[3]+g[3]*g[0])/2,complex)
    forward=np.zeros((8,8),complex);forward[:4,4:]=S
    hop=_lift(forward)
    shift=sparse.csr_matrix([[0,0,0],[1,0,0],[0,1,0]],dtype=complex)
    T=(sparse.kron(shift,hop)+sparse.kron(shift.conj().T,hop.conj().T)).toarray()
    q=np.tile(np.array([((x&15).bit_count()-2)/30 for x in range(256)]),3)
    dq=1/30
    double=-(q[:,None]-q[None,:])**2*T
    assert np.max(abs(double+dq*dq*T))<1e-14
    r=F(3,2);phase=np.pi*float(r)*q
    c=np.cos(phase);s=np.sin(phase)
    changed=c[:,None]*T*c[None,:]+s[:,None]*T*s[None,:]-T
    assert np.max(abs(changed-(np.cos(np.pi*float(r)*dq)-1)*T))<1e-13


def test_member_zero_work_uses_all_eight_origins_and_twelve_real_bonds():
    r=7/3;total=0.
    for X in product(range(2),repeat=3):
        phase=np.zeros((2,2,2));phase[X]=np.pi*r/30
        total+=sum(np.sum(1-np.cos(np.diff(phase,axis=j))) for j in range(3))/8
    k=2*(1-np.cos(np.pi*r/30))
    assert abs(total-12*k/8)<1e-14
    assert 4*r*r/900<=k<=10*r*r/900


def test_arbitrary_entangled_control_uses_the_same_double_commutator():
    rng=np.random.default_rng(811)
    X=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));H=X.conj().T@X
    H-=np.linalg.eigvalsh(H)[0]*np.eye(4)
    Q=np.diag([-1.,-.2,.3,1.]);A=np.array([[0.,1.],[1.,0.]])
    Hjoint=np.kron(H,np.eye(2));Fjoint=np.pi*.7*np.kron(Q,A)
    local=Q@(H@Q-Q@H)-(H@Q-Q@H)@Q
    double=Fjoint@(Hjoint@Fjoint-Fjoint@Hjoint)-(Hjoint@Fjoint-Fjoint@Hjoint)@Fjoint
    assert np.max(abs(double-(np.pi*.7)**2*np.kron(local,A@A)))<1e-12
    psi=rng.normal(size=8)+1j*rng.normal(size=8);psi/=np.linalg.norm(psi)
    U=expm(1j*Fjoint);before=np.vdot(psi,Hjoint@psi).real
    after=np.vdot(U@psi,Hjoint@(U@psi)).real
    assert after<=(np.sqrt(before)+np.sqrt(np.linalg.norm(double,2)/2))**2+1e-11


def test_complete_pointer_and_waiting_outputs_retain_entangled_reference():
    H=np.diag([0.,.4,2.]);Q=np.array([[0.,.2,0.],[.2,.1,.3],[0.,.3,-.1]])
    rho=np.diag([.97,.02,.01]);root=np.sqrt(rho)
    # Source x reference x pointer; coherent signs keep all correlations.
    psi=np.kron(root.reshape(-1),np.array([1.,1.])/np.sqrt(2))
    K=np.kron(np.kron(H,np.eye(3)),np.eye(2));r=.6
    Fjoint=np.pi*r*np.kron(np.kron(Q,np.eye(3)),np.diag([1.,-1.]))
    U=expm(1j*Fjoint);state=U@psi
    e=np.vdot(state,K@state).real;t=.37;waited=expm(-1j*t*K)@state
    distance=np.sqrt(max(0,1-abs(np.vdot(state,waited))**2))
    assert distance<=np.sqrt(2*t*e)+1e-14
    # Measuring a complete pointer retains total probability and reference.
    tensor=state.reshape(3,3,2)
    assert abs(sum(np.vdot(tensor[:,:,j],tensor[:,:,j]).real for j in range(2))-1)<1e-14
    reference=np.einsum('sri,sqi->rq',tensor,tensor.conj())
    assert np.max(abs(reference-rho))<1e-13
    final_prob=np.sum(abs(tensor[2])**2)
    assert final_prob<=e/2+1e-14


def test_controlled_ground_phase_is_not_an_uncontrolled_wait():
    E0=-3.;time=.2;pointer=np.array([1.,1.])/np.sqrt(2)
    changed=expm(-1j*time*E0*np.diag([1.,-1.]))@pointer
    distance=np.sqrt(1-abs(np.vdot(pointer,changed))**2)
    assert distance>.5
    # A one-dimensional source has zero ground-relative energy, while
    # its controlled absolute-energy phase is visible in the pointer.


@pytest.fixture(scope='module')
def receipts():
    return {n:m.bind_member(m.contract.member(n)) for n in (0,1,2,8)}


def test_actual_source_fixed_strength_obstruction_and_critical_boundary(receipts):
    last=1
    for n,bound in receipts.items():
        r=bound['result'];h=r['fixed_history'];critical=r['critical_Weyl'];a=F(1,2**n)
        cap=F(h['original_positive_transfer_probability_upper'])
        assert cap<=last;last=cap
        assert F(h['normalization']['double_commutator_norm_upper'])<=F(9,80)*a**4
        assert F(critical['original_Weyl_strength'])*a*a==30
        assert critical['energy_root_drift_upper']=='45/2'
        assert not critical['subcritical_obstruction_applies']
        assert not critical['fixed_positive_energy_survival_decided']
        assert h['complete_reference_and_pointer_outputs_retained']
        assert not h['controlled_waiting_time_indistinguishability_claimed']
        assert h['physical_mass_table_GeV']==[]
        low,high=map(F,r['balanced_finite_work']['original_signed_energy_mean_interval'])
        assert 0<low<=high<=F(9,16)*a**4
    assert last<F(1,4000)
    assert receipts[8]['result']['necessary_strength']['unit_bounded_phase_calls_necessary_lower']>7000


def test_relative_root_bounds_do_not_leave_a_fixed_rounding_floor():
    for n in (0,1,2,63,64,129):
        root=F(1,2**((n+1)//2+2))*(F(3,2) if n%2 else 1)
        eps=F(1,2**(n+4))
        assert eps<=root*root<=F(9,8)*eps


def test_domains_and_forged_full_source_no_go_are_rejected(receipts):
    member=m.contract.member(0);m.verify_bound(member,receipts[0])
    for strengths,threshold in [((), '1'),((True,), '1'),((.2,), '1'),(('1',),'0')]:
        with pytest.raises(ValueError):m.fixed_history(member,strengths,threshold)
    with pytest.raises(ValueError):m.minimum_strength(member,'2')
    with pytest.raises(ValueError):m.finite_work(member,'31')
    wrong=deepcopy(receipts[0]);wrong['result']['critical_Weyl']['fixed_positive_energy_survival_decided']=True
    wrong['result_sha256']=m.contract.digest(wrong['result'])
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
    changed=deepcopy(member);changed['clock']['Schrodinger']='exp(-2itH)'
    with pytest.raises(ValueError):m.bind_member(changed)
