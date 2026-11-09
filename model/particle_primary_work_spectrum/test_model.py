"""Independent finite transition, CAR and literal polynomial controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json
import math

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm

from . import model as m


def _transition_control():
    energies=np.array([-1.,-.2,.7,2.])
    beta=1.3;order=4;scale=8
    y=beta*(energies[-1]-energies)/(2*scale)
    Q=sum(y**k/math.factorial(k) for k in range(order+1))**scale
    p=Q**2/np.sum(Q**2)
    F0=np.array([[.1,.4,.2j,0.],[.4,-.3,.1,.2],[-.2j,.1,.2,.3],[0.,.2,.3,.4]])
    units=[expm(1j*F0),expm(-1j*F0)]
    # T has incoming state on the first index, as in Lean push.
    T=sum(abs(U.T)**2 for U in units)/2
    paths=p[:,None]*T
    omega=energies[None,:]-energies[:,None]
    kappa=np.log(p[:,None]/p[None,:])
    return energies,p,units,T,paths,omega,kappa,beta,order,scale


def test_actual_finite_fluctuation_theorem_premises_and_distinct_entropy_variable():
    E,p,units,T,P,omega,kappa,beta,order,scale=_transition_control()
    assert np.max(abs(T-T.T))<1e-15
    assert np.max(abs(T.sum(axis=1)-1))<1e-15
    pi=np.full(4,.25);q=p@T
    assert np.min(q)>0 and np.max(abs(pi@T-pi))<1e-15
    sigma=np.log(p[:,None]/q[None,:])
    assert np.max(abs(sigma-kappa))>.01
    assert np.sum(P*np.exp(-sigma))==pytest.approx(1.,abs=2e-15)
    assert np.sum(P*np.exp(-kappa))==pytest.approx(1.,abs=2e-15)
    assert np.max(abs(P-np.exp(sigma)*(q[:,None]*T).T))<2e-16
    H=lambda x: -np.sum(x*np.log(x))
    kl=np.sum(q*np.log(q/p))
    assert np.sum(P*sigma)==pytest.approx(H(q)-H(p),abs=3e-16)
    assert np.sum(P*kappa)==pytest.approx(H(q)-H(p)+kl,abs=3e-16)
    slope=m.modular.polynomial_slope(str(beta),str(E[0]),str(E[-1]),order,scale)
    W=np.sum(P*omega);K=np.sum(P*kappa)
    assert float(slope['modular_slope_lower'])*W<=K<=beta*W


def test_retained_signed_instrument_reference_and_affinity_domain():
    E,p,units,T,P,omega,kappa,*_=_transition_control()
    root=np.diag(np.sqrt(p));rho=np.diag(p)
    reference=np.zeros((4,4),complex);effect=np.zeros((4,4),complex)
    affinity=0.;distance=0.
    for U in units:
        K=U/np.sqrt(2);effect+=K.conj().T@K
        branch=K@root;reference+=branch.T@branch.conj()
        affinity+=np.trace(root@U@root@U.conj().T).real/2
        distance+=np.sum(abs(np.linalg.eigvalsh(U@rho@U.conj().T-rho)))/4
    assert np.max(abs(effect-np.eye(4)))<1e-15
    assert np.max(abs(reference-rho.T))<1e-15
    assert affinity==pytest.approx(np.sum(P*np.exp(-kappa/2)),abs=2e-15)
    assert 0<=affinity<=1
    W=np.sum(P*omega);width=E[-1]-E[0]
    slope=m.modular.polynomial_slope('1.3','-1','2',4,8)
    low=float(slope['modular_slope_lower'])
    assert 1-affinity>=low*W/(4+low*width)
    assert distance>=W/width
    # A diagonal phase is invisible on the source but changes its retained
    # canonical purification. Its overlap is not source quantum affinity.
    U=np.diag([1.,-1.,1.,-1.])
    A=np.trace(root@U@root@U.conj().T).real
    joint_overlap_squared=abs(np.trace(rho@U))**2
    assert A==pytest.approx(1.) and joint_overlap_squared<.8


def test_negative_tail_pairing_and_outward_nonzero_factored_bound():
    E,p,units,T,P,omega,kappa,*_=_transition_control()
    low=float(m.modular.polynomial_slope('1.3','-1','2',4,8)['modular_slope_lower'])
    for cutoff in (.1,.9,1.4,2.5):
        assert P[omega<=-cutoff].sum()<=1/(1+np.exp(low*cutoff))
    negative=np.sum(P*np.maximum(-omega,0))
    assert negative<=1/(3*low)
    absolute=np.sum(P*abs(omega));W=np.sum(P*omega)
    assert absolute==pytest.approx(W+2*negative)
    row=m.downward_tail(F(22271),F(1,2))
    assert row['rational_upper_factored']['exponent']==11135
    assert row['rational_upper_factored']['numerator']==1
    with pytest.raises(ValueError):m.downward_tail(0,1)
    with pytest.raises(ValueError):m.downward_tail(1,0)


def test_full_CAR_trace_norm_for_literal_Wilson_bonds():
    gamma=[np.array(x.tolist(),dtype=complex) for x in m.density.gamma_data()]
    beta=gamma[3];theta=.37;a=.5
    for g in gamma[:3]:
        S=(-beta+beta@g)/2
        assert np.linalg.svd(S,compute_uv=False)==pytest.approx([1,1,0,0],abs=2e-15)
        z=np.exp(1j*theta)-1
        t=np.block([[np.zeros((4,4)),z*S/a],[z.conjugate()*S.conj().T/a,np.zeros((4,4))]])
        ev=np.linalg.eigvalsh(t)
        fock=np.array([sum(v for v,bit in zip(ev,bits) if bit) for bits in product((0,1),repeat=8)])
        assert np.max(abs(fock))==pytest.approx(np.sum(abs(ev))/2)
        assert np.trace(t@t).real==pytest.approx(8*(1-np.cos(theta))/a**2)
        # The fifteen total internal rows contribute their dimensions,
        # independently of the link's unitary gauge chart.
        native_phases=np.exp(1j*np.arange(15)/7)
        native=np.diag(native_phases)
        block=np.kron(z*S/a,native)
        full=np.block([[np.zeros((60,60)),block],[block.conj().T,np.zeros((60,60))]])
        assert np.trace(full@full).real==pytest.approx(120*(1-np.cos(theta))/a**2)
        assert np.max(abs(fock))**2<=np.linalg.matrix_rank(t)*np.trace(t@t).real/4+1e-12


def test_parent_origin_support_and_gradient_compression_bound():
    R=2;L=2*R*R;a=1/R;coords=np.indices((L,L,L));f=np.ones((L,L,L))
    for i in range(3):f*=np.maximum(1-np.minimum(coords[i],L-coords[i])/R,0)
    active=np.zeros_like(f,dtype=bool);cosine_sum=0.;d2=0.
    for i in range(3):
        delta=np.roll(f,-1,axis=i)-f
        changes=abs(delta)>0
        active|=changes|np.roll(changes,1,axis=i)
        cosine_sum+=np.sum(2*np.sin(np.pi*a*delta/2)**2)
        d2+=np.sum(delta**2)
    assert active.sum()==(2*R-1)**2*(2*R+5)
    row=m.current_second_moment(m.contract.member(1))
    assert row['active_source_sites_upper']==active.sum()
    assert d2==pytest.approx(float(m.critical.profile_geometry(1)['three_axis_gradient_square_sum']))
    assert cosine_sum/3<=5*a*a*d2/3
    # A direct full-Fock norm control uses the one-body signed spectrum.
    # Arbitrary orthogonal compression can only reduce its operator norm.
    rng=np.random.default_rng(71);Z=rng.normal(size=(10,10))+1j*rng.normal(size=(10,10))
    A=Z+Z.conj().T;A-=np.trace(A)*np.eye(10)/10
    Q=np.linalg.qr(rng.normal(size=(10,4))+1j*rng.normal(size=(10,4)))[0]
    assert np.linalg.norm(Q.conj().T@A@Q,2)<=np.linalg.norm(A,2)+1e-12


def test_second_moment_improves_all_integer_refinements_but_is_not_uniform():
    R=sp.Symbol('R',positive=True);t=sp.Symbol('t',nonnegative=True)
    sites=(2*R-1)**2*(2*R+5)
    numerator=108*R**9-sites*(2*R**2-1)*(2*R**2+1)**2
    shifted=sp.Poly(sp.expand(numerator.subs(R,t+2)),t)
    assert all(x>0 for x in shifted.all_coeffs())
    limit=sp.limit(3000*sites*(2*R**2-1)*(2*R**2+1)**2/R**9,R,sp.oo)
    assert limit==192000
    for n in (1,2,8):
        row=m.current_second_moment(m.contract.member(n))
        assert row['selected_second_moment_upper']<row['inherited_second_moment_upper']
        assert row['selected_second_moment_upper']>F(190000*2**(4*n))
        assert not row['uniform_second_moment_proved']
    assert m.current_second_moment(m.contract.member(0))['selected_second_moment_upper']==19620


def test_literal_primary_polynomial_control_is_faithful_and_has_elastic_limit():
    works=[];affinities=[]
    with mp.workdps(80):
        for n in (0,1,2,8,12):
            member=m.contract.member(n);row=m.polynomial_control(member)
            data=member['state'];beta=mp.mpf(str(data['beta']));scale=int(data['dyadic_scale'])
            a=mp.mpf(2)**(-n);y=beta/(2*scale*a);order=data['Taylor_order']
            increment=sum((y**k/mp.factorial(k) for k in range(1,order+1)),mp.mpf(0))
            ell=2*scale*mp.log1p(increment)
            lower=F(row['negative_log_excited_to_ground_probability_ratio_interval'][0])
            lower_mp=mp.mpf(lower.numerator)/lower.denominator
            assert ell>=lower_mp*(1-mp.mpf('1e-70'))
            assert ell<=beta/a*(1+mp.mpf('1e-70'))
            ratio=mp.exp(-ell);pg=1/(1+ratio);pe=ratio/(1+ratio)
            assert ratio>0 and pe>0
            assert pe/a<=1/lower_mp/a
            phase_probability=mp.sin(mp.sqrt(a))**2
            W=(pg-pe)*phase_probability/a;M2=phase_probability/a**2
            A=1-phase_probability*(mp.sqrt(pg)-mp.sqrt(pe))**2
            works.append(float(W));affinities.append(float(A))
            assert 0<W<1 and M2>=W/a
            assert 0<phase_probability<=a
            assert M2*a>=1-a/3
            assert not row['same_carrier_or_source_as_primary']
            # Subtracting the declared scalar ground energy changes no commutator.
            H=np.diag([0.,2.**n]);F0=np.sqrt(2.**(-n))*np.array([[0.,1.],[1.,0.]])
            comm=H@F0-F0@H;double=F0@comm-comm@F0
            assert np.linalg.norm(double,2)==pytest.approx(2.)
            for time in (.1,1.,10.):
                U=expm(1j*F0);flow=np.diag(np.exp(1j*time*np.diag(H)))
                assert np.linalg.norm(flow@U@flow.conj().T-np.eye(2),2)<=np.sqrt(2.**(-n))+1e-14
        assert all(x<y for x,y in zip(works,works[1:]))
        assert all(x<y for x,y in zip(affinities,affinities[1:]))
        assert works[-1]>.9999 and affinities[-1]>.9997


def test_actual_member_bindings_clock_budget_and_unpromoted_physical_scope():
    for n in (0,1,2,8):
        member=m.contract.member(n);bound=m.bind_member(member);m.verify_bound(member,bound)
        row=bound['result'];slope=row['actual_modular_slope'];b=F(slope['modular_slope_lower'])
        assert F(row['negative_absolute_first_moment_upper'])==1/(3*b)
        assert F(row['negative_absolute_first_moment_upper'])<F(member['state']['target_ground_excess'])
        assert F(row['source_CQ_half_trace_distance_from_no_operation_lower'])>0
        assert F(row['source_CQ_affinity_deficit_lower'])>0
        assert F(row['finite_complete_clock_record_total_variation_lower'])>0
        ideal=row['ideal_original_positive_window'];reported=row['finite_reported_positive_window']
        eta=F(reported['energy_resolution']);loss=F(reported['two_complete_clock_probability_loss_upper'])
        assert F(reported['reported_transfer_interval'][0])==F(ideal['original_transfer_interval'][0])-2*eta>0
        assert F(reported['reported_transfer_interval'][1])==F(ideal['original_transfer_interval'][1])+2*eta
        assert F(reported['original_probability_lower'])==F(ideal['original_weight_lower'])-loss>F(ideal['original_weight_lower'])/2
        assert not reported['conditioned_or_postselected_probability']
        assert not row['affinity_is_untouched_joint_purification_overlap']
        assert not row['modular_work_and_classical_entropy_production_are_equal']
        assert not row['nonelastic_fixed_energy_limit_proved']
        assert not row['uniform_selected_energy_weighted_tail_proved']
        assert row['physical_mass_table_GeV']==[]
        assert json.loads(json.dumps(bound,sort_keys=True,allow_nan=False))==bound
    bad=deepcopy(bound);bad['result']['original_signed_work_interval'][0]='1'
    with pytest.raises(ValueError):m.verify_bound(member,bad)
    bad_member=deepcopy(member);bad_member['state']['beta']='1'
    with pytest.raises(ValueError):m.bind_member(bad_member)


def test_signed_mean_window_and_complete_clock_error_have_original_normalization():
    values=np.array([-10.,-.5,0.,.1,1.,3.,20.]);prob=np.array([.001,.05,.60,.04,.2,.1,.009])
    assert prob.sum()==pytest.approx(1.)
    mean=float(prob@values);second=float(prob@(values**2))
    w=F(str(mean));cap=F(str(second));window=m.local.positive_window(w,cap)
    lo,hi=map(float,window['original_transfer_interval'])
    actual=prob[(values>=lo)&(values<=hi)].sum()
    assert actual>=float(window['original_weight_lower'])
    # A complete noisy record is a stochastic kernel: retain every failure,
    # then bound a widened event by the bad-clock probability, without
    # renormalizing on the positive records.
    eta=.01;tail=.001;outcomes=[];weights=[]
    for energy,p in zip(values,prob):
        outcomes.extend([energy-2*eta,energy+2*eta,-100.])
        weights.extend([p*(1-2*tail)/2,p*(1-2*tail)/2,p*2*tail])
    outcomes=np.array(outcomes);weights=np.array(weights)
    assert weights.sum()==pytest.approx(1.)
    reported=weights[(outcomes>=lo-2*eta)&(outcomes<=hi+2*eta)].sum()
    assert reported>=actual-2*tail
