"""Independent source CAR, form-domain and complete-history controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
from scipy.linalg import expm
from scipy import sparse

from . import model as m
from model.particle_local_wilson_density import density
from model.particle_primary_response_moments import model as response


def _car(modes,i,j):
    rows=[];cols=[];values=[]
    for state in range(1<<modes):
        if not state&(1<<j):continue
        sign=(-1)**((state&((1<<j)-1)).bit_count());after=state^(1<<j)
        if after&(1<<i):continue
        sign*=(-1)**((after&((1<<i)-1)).bit_count())
        rows.append(after|(1<<i));cols.append(state);values.append(sign)
    return sparse.csr_matrix((values,(rows,cols)),shape=(1<<modes,1<<modes),dtype=complex)


def _lift(onebody):
    total=sparse.csr_matrix((1<<len(onebody),1<<len(onebody)),dtype=complex)
    for i,j in zip(*np.nonzero(onebody)):
        total+=onebody[i,j]*_car(len(onebody),int(i),int(j))
    return total


def test_literal_full_CAR_link_double_commutator_and_sign_average():
    g=density.gamma_data();beta=np.array(g[3],complex)
    block=np.array((-g[3]+g[3]*g[0])/2,complex)
    forward=np.zeros((8,8),complex);forward[:4,4:]=block
    onsite=np.zeros((8,8),complex);onsite[:4,:4]=3*beta;onsite[4:,4:]=3*beta
    hop=_lift(forward);h0=_lift(onsite)
    shift=sparse.csr_matrix([[0,0,0],[1,0,0],[0,1,0]],dtype=complex)
    J=(sparse.kron(shift,hop)+sparse.kron(shift.conj().T,hop.conj().T)).toarray()
    electric=sparse.diags([15.75,0,15.75],dtype=complex)
    H=J+(sparse.kron(sparse.eye(3),h0)+sparse.kron(electric,sparse.eye(256))).toarray()
    n0=np.array([(x&15).bit_count()-2 for x in range(256)])
    n1=np.array([(x>>4).bit_count()-2 for x in range(256)])
    phase=np.tile(np.pi*(n0/3-n1/7),3)
    delta=np.pi*(1/3+1/7)
    HF=H*(phase[None,:]-phase[:,None])
    double=(phase[:,None]-phase[None,:])*HF
    assert np.max(abs(double+delta**2*J))<1e-13
    u=np.exp(1j*phase)
    plus=u.conj()[:,None]*H*u[None,:]
    minus=u[:,None]*H*u.conj()[None,:]
    assert np.max(abs((plus+minus)/2-H-(np.cos(delta)-1)*J))<1e-13


def test_semibounded_energy_form_lemma_on_noncommuting_full_matrices():
    rng=np.random.default_rng(703)
    A=rng.normal(size=(7,7))+1j*rng.normal(size=(7,7));K=A.conj().T@A
    K-=np.linalg.eigvalsh(K)[0]*np.eye(7)
    F=rng.normal(size=(7,7))+1j*rng.normal(size=(7,7));F=(F+F.conj().T)/2
    double=F@(K@F-F@K)-(K@F-F@K)@F
    cap=np.linalg.norm(double,2);drift=np.sqrt(cap/2)
    for time in (-3.,-.25,.2,1.3):
        U=expm(1j*time*F)
        for _ in range(4):
            psi=rng.normal(size=7)+1j*rng.normal(size=7);psi/=np.linalg.norm(psi)
            before=max(0,np.vdot(psi,K@psi).real)
            after=np.vdot(U@psi,K@(U@psi)).real
            assert after<=(np.sqrt(before)+abs(time)*drift)**2+1e-10


def test_filtered_contraction_root_commutator_and_mixed_word_bound():
    rng=np.random.default_rng(704);e=np.array([0.,.1,1.,7.,30.]);K=np.diag(e)
    A=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5));A/=np.linalg.norm(A,2)
    T=.7;C=np.sinc((e[:,None]-e[None,:])*T/(2*np.pi))**2*A
    d=2/T
    assert np.linalg.norm(C,2)<=1+1e-14
    assert np.linalg.norm(K@C-C@K,2)<=d
    assert np.linalg.norm(np.diag(np.sqrt(e))@C-C@np.diag(np.sqrt(e)),2)<=2*np.sqrt(d)
    F=rng.normal(size=(5,5));F=(F+F.T)/2
    double=F@(K@F-F@K)-(K@F-F@K)@F
    r=np.sqrt(np.linalg.norm(double,2)/2);U=expm(1j*F)
    p=np.exp(-20*e);p/=p.sum();rho=np.diag(p);epsilon=np.trace(rho@K).real
    W=U@C@expm(.3j*K)@U.conj().T@expm(-.3j*K)
    after=np.trace(rho@W.conj().T@K@W).real
    assert after<=(np.sqrt(epsilon)+2*r+2*np.sqrt(d))**2
    # The response uses the stationary state and original contraction mass.
    weights=p[None,:]*abs(W)**2;frequencies=e[:,None]-e[None,:]
    absolute=np.sum(weights*abs(frequencies))
    assert absolute<=after+epsilon+1e-12
    for time in (.01,.1,.7):
        translated=np.exp(1j*time*frequencies)*W
        distance=np.trace(rho@(translated-W).conj().T@(translated-W)).real
        assert distance<=2*abs(time)*absolute+1e-12


def test_complete_random_sign_histories_preserve_reference_and_energy_bound():
    H=np.array([[0.,.3],[.3,2.]]);E0=np.linalg.eigvalsh(H)[0];K=H-E0*np.eye(2)
    F=np.diag([.4,-.2]);U=[expm(-1j*F),expm(1j*F)]
    even=sum(V.conj().T@K@V for V in U)/2-K
    cap=np.linalg.norm(even,2)
    rng=np.random.default_rng(705);psi=rng.normal(size=(2,3))+1j*rng.normal(size=(2,3));psi/=np.linalg.norm(psi)
    initial=np.trace(psi.conj().T@K@psi).real
    rows=[psi]
    for step in range(3):
        time=expm(-1j*(step+1)*.17*H)
        rows=[V@time@row/np.sqrt(2) for row in rows for V in U]
    assert len(rows)==8
    assert max(abs(np.vdot(row,row).real-1/8) for row in rows)<1e-14
    reference=sum((row.conj().T@row for row in rows),np.zeros((3,3),complex))
    assert np.max(abs(reference-psi.conj().T@psi))<1e-13
    energy=sum(np.trace(row.conj().T@K@row).real for row in rows)
    assert energy<=initial+3*cap+1e-13


@pytest.mark.parametrize('radius',[1,2])
def test_all_actual_open_bonds_and_origins_give_anchor_formula(radius):
    R=radius;L=2*R*R;a=1/R
    t=np.array([max(0,1-min(k,L-k)/R) for k in range(L)])
    f=t[:,None,None]*t[None,:,None]*t[None,None,:]
    k=np.sum(1-np.cos(np.pi*a*(np.roll(f,1,axis=0)-f)))
    exact=0
    for X in product(range(L),repeat=3):
        phase=np.pi*a*np.roll(f,X,axis=(0,1,2))
        exact+=30/a**2*sum(np.sum(1-np.cos(np.diff(phase,axis=axis))) for axis in range(3))/L**3
    target=30*(3*L*L*(L-1))*k/(a*a*L**3)
    assert abs(exact-target)<1e-9
    assert exact>=160/a
    if R==1:assert abs(exact-180)<1e-12


def test_actual_corner_CAR_action_retains_both_orientations():
    pair=response.beta_pairs()
    chart=np.column_stack([np.array(pair[0][0],complex).ravel(),np.array(pair[1][0],complex).ravel(),
                           np.array(pair[0][1],complex).ravel(),np.array(pair[1][1],complex).ravel()])
    g=density.gamma_data();negative=np.diag([0.,0.,1.,1.])
    for axis in range(3):
        S=chart.conj().T@np.array((-g[3]+g[3]*g[axis])/2,complex)@chart
        forward=np.zeros((8,8),complex);forward[:4,4:]=S
        h=_lift(forward)+_lift(forward.conj().T)
        omega=np.zeros(256);omega[(1<<2)|(1<<3)|(1<<6)|(1<<7)]=1
        packet=h@omega
        assert abs(np.vdot(packet,packet)-1)<1e-13
        # A pi phase at either endpoint changes every hopping packet sign.
        phase=np.array([np.pi*((mask&15).bit_count()-2) for mask in range(256)])
        result=packet-np.exp(1j*phase)*packet
        assert abs(np.vdot(result,result)-4)<1e-13


def test_uniform_double_commutator_can_have_an_elastic_limit():
    for n in (3,5,8):
        a=4.**(-n);H=np.diag([0.,1/a]);F=np.pi*np.sqrt(a)*np.array([[0.,1.],[1.,0.]])
        double=F@(H@F-F@H)-(H@F-F@H)@F
        assert abs(np.linalg.norm(double,2)-2*np.pi**2)<1e-12
        U=expm(1j*F);weight=abs(U[1,0])**2
        assert weight<=np.pi**2*a
        assert weight/a>8
        assert weight/a**2>8/a
        if n>=5:assert 1/a>100


@pytest.fixture(scope='module')
def receipts():
    return {n:m.bind_member(m.contract.member(n)) for n in (0,1,2,8)}


def test_actual_source_mixed_words_and_operator_obstruction_are_scoped(receipts):
    for n,value in receipts.items():
        r=value['result'];word=r['mixed_word'];pure=r['pure_Weyl_word'];anchor=r['physical_anchor']
        assert F(word['source_form_bound']['double_commutator_norm_upper'])<=F(2025,2)
        assert F(word['fixed_word_final_ground_relative_energy_mean_upper'])<=F(word['uniform_final_ground_relative_energy_mean_upper'])
        assert word['common_unfiltered_Weyl_and_filtered_probe_time_GNS_constructed']
        assert word['canonical_GNS_energy_nonnegative']
        assert not word['original_spectral_mass_is_exactly_one']
        assert pure['original_spectral_mass_is_exactly_one']
        assert not word['point_norm_continuity_of_unfiltered_sections_proved']
        assert not word['spatial_continuity_or_momentum_for_unfiltered_Weyls_proved']
        assert not word['uniform_fixed_positive_energy_weight_proved']
        assert F(anchor['averaged_anchor_commutator_square_interval'][0])>=F(160)*2**n
        assert not anchor['anchor_is_selected_state']
        assert not anchor['selected_state_weighted_second_moment_lower_inferred']
        assert word['physical_mass_table_GeV']==[]
    assert F(receipts[0]['result']['physical_anchor']['averaged_anchor_commutator_square_interval'][0])==180


def test_actual_source_complete_modular_response_is_uniform_and_marginal(receipts):
    for value in receipts.values():
        row=value['result']['source_marginal_modular_response']
        delta=F(row['relative_modular_slope_defect_upper'])
        low,high=map(F,row['normalized_complete_CQ_relative_entropy_interval'])
        uniform_low,uniform_high=map(F,row['uniform_normalized_complete_CQ_relative_entropy_interval'])
        assert 0<delta<=F(1,3)
        assert 0<uniform_low<=low<=high<=uniform_high
        assert row['all_sign_origin_outcomes_included']
        assert not row['source_state_replaced_by_exact_Gibbs']
        assert not row['source_reference_joint_relative_entropy_identified']
        assert not row['regional_diamond_modular_kernel_identified']
        assert not row['response_mean_passes_to_common_GNS_limit_proved']


def test_polynomial_CQ_entropy_bracket_with_all_original_sign_records():
    energies=np.array([0.,.3,1.]);beta=2.;order=2;scale=4
    y=beta*(1-energies)/(2*scale)
    q=(1+y+y*y/2)**scale;p=q*q/np.sum(q*q);K=-np.log(p)
    slope=m.modular.polynomial_slope('2','0','1',order,scale)
    delta=float(slope['relative_modular_slope_defect_upper'])
    F=np.array([[0.,.4,.1],[.4,0.,.7],[.1,.7,0.]])
    work=0.;relative=0.
    for sign in (-1,1):
        U=expm(sign*1j*F);changed=U@np.diag(p)@U.conj().T
        difference=np.diag(changed).real-p
        work+=np.dot(difference,energies)/2
        # Source entropy is unchanged by the unitary. Identical original
        # sign probabilities cancel the classical relative-entropy terms.
        relative+=np.dot(difference,K)/2
    assert work>0
    assert beta*(1-delta)*work<=relative<=beta*work


def test_source_word_domains_and_forged_survival_are_rejected(receipts):
    member=m.contract.member(0);m.verify_bound(member,receipts[0])
    for k,widths in [(True,()),(-1,()),(0,()),(1,('0',)),(1,(.5,))]:
        with pytest.raises(ValueError):m.mixed_word_bound(member,k,widths)
    wrong=deepcopy(receipts[0]);wrong['result']['mixed_word']['uniform_fixed_positive_energy_weight_proved']=True
    wrong['result_sha256']=m.contract.digest(wrong['result'])
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
    changed=deepcopy(member);changed['state']['target_ground_excess']='1'
    with pytest.raises(ValueError):m.bind_member(changed)
