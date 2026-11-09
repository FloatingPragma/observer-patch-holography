"""Independent complete-channel, joint-kernel and authenticated-source tests."""
from copy import deepcopy
from fractions import Fraction as F
import math

import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from . import model as m


def test_nonhermitian_weighted_selector_has_all_records_and_original_coherence():
    value=m.exact_control();instrument=value['instrument']
    assert instrument['complete_effect']==s.eye(2)
    assert instrument['fine_outcomes']==8
    assert value['reference_difference']==s.zeros(2)
    assert sum(value['probabilities'])==1
    assert value['original_success']!=value['incoherent_success']
    assert value['binary_complete']==s.eye(2)
    W=np.array(instrument['analyzer'],dtype=complex)
    np.testing.assert_allclose(W.conj().T@W,np.eye(4),atol=1e-14)
    Ks=[np.array(row['K'],dtype=complex) for row in instrument['rows']]
    rho=np.diag([.75,.25])
    assert sum(np.trace(rho@K.conj().T@K).real for K in Ks)==pytest.approx(1)
    assert np.trace(rho@Ks[0].conj().T@Ks[0]).real==pytest.approx(float(value['original_success']))


def test_complete_binary_failure_is_not_the_fine_lcu_failure_channel():
    value=m.exact_control();rows=value['instrument']['rows'];B=value['instrument']['success']
    root=m.binary_instrument(B)[1];rho=s.Matrix([[s.Rational(1,2),s.Rational(1,4)],[s.Rational(1,4),s.Rational(1,2)]])
    fine=sum((r['K']*rho*r['K'].H for r in rows[1:]),s.zeros(2))
    binary=root*rho*root.H
    assert s.simplify(s.trace(fine-binary))==0
    assert m.apparatus.clean(fine-binary)!=s.zeros(2)


def test_cartesian_fourier_analyzer_keeps_every_bin_and_matches_the_coherent_kernel():
    B=s.Matrix([[0,s.Rational(1,2)],[0,0]])
    groups=((B,B.T),(B.T,s.I*B),(-B,B.T),(s.I*B,-B.T))
    weights=(F(3,4),F(1,4));value=m.factored_fourier_instrument(groups,weights,(2,2))
    assert value['fine_outcomes']==16 and value['complete_effect']==s.eye(2)
    selected=next(row['K'] for row in value['rows']
        if (row['branch'],row['inner_analyzer'],row['Fourier_index'])==(0,0,(1,0)))
    expected=sum((sign*(s.Rational(3,4)*group[0]+group[1]/4)
        for sign,group in zip((1,1,-1,-1),groups)),s.zeros(2))/4
    assert m.apparatus.clean(selected-expected)==s.zeros(2)
    rho=s.diag(s.Rational(2,3),s.Rational(1,3))
    assert s.simplify(sum(s.trace(rho*row['K'].H*row['K']) for row in value['rows']))==1


def test_exact_weight_analyzer_rejects_hidden_success_normalization():
    B=s.Matrix([[0,1],[0,0]])
    with pytest.raises(ValueError):m.weighted_instrument((B,B),(F(1,4),F(1,4)))
    with pytest.raises(ValueError):m.weighted_instrument((B,),(1,),(2,))
    value=m.weighted_instrument((B,),(1,))
    assert value['success']==B and value['complete_effect']==s.eye(2)


def test_joint_bochner_identity_uses_both_frequency_signs_and_original_mass():
    H=np.diag([0.,2.,5.]);P=np.diag([0.,-1.,3.]);Omega=np.array([1.,0.,0.])
    A=np.array([[0.,0.,0.],[.3,0.,0.],[.4,0.,0.]])
    Q=.9;S=.7;q=.8;nu=1.2;N=800
    spatial=np.zeros_like(A,dtype=complex);temporal=np.zeros_like(A,dtype=complex)
    for y in -Q+(np.arange(N)+.5)*2*Q/N:
        V=expm(1j*y*P);spatial+=np.exp(-1j*q*y)*V@A@V.conj().T/N
    for t in -S+(np.arange(N)+.5)*2*S/N:
        U=expm(1j*t*H);temporal+=np.exp(-1j*nu*t)*U@spatial@U.conj().T/N
    predicted=np.sinc(S*(np.diag(H)-nu)/np.pi)*np.sinc(Q*(np.diag(P)-q)/np.pi)*(A@Omega)
    np.testing.assert_allclose(temporal@Omega,predicted,atol=2e-7)
    probability=np.vdot(temporal@Omega,temporal@Omega).real
    # The midpoint rule error for exp(i omega x) is omega² Q²/(6N²).
    component_error=((Q*(np.diag(P)-q))**2+(S*(np.diag(H)-nu))**2)/(6*N*N)
    amplitude_error=np.linalg.norm(component_error*(A@Omega))
    assert np.linalg.norm(temporal@Omega-predicted)<=amplitude_error
    assert abs(probability-np.vdot(predicted,predicted).real)<=2*amplitude_error
    assert np.linalg.norm(A@Omega)**2==pytest.approx(.25)
    wrong=np.sinc(S*(np.diag(H)+nu)/np.pi)*np.sinc(Q*(np.diag(P)+q)/np.pi)*(A@Omega)
    assert np.linalg.norm(predicted-wrong)>.01


def test_noninvariant_input_cannot_use_the_vacuum_fixed_identity():
    # A diagonal operator is fixed by beta, while this vector is not fixed
    # by its spatial unitary. The invariant relative-state hypothesis matters.
    P=np.diag([0.,1.]);Omega=np.ones(2)/np.sqrt(2);A=np.diag([0.,1.])
    Q=1.;actual=A@Omega;false_prediction=np.sinc(Q*np.diag(P)/np.pi)*(A@Omega)
    assert np.linalg.norm(actual-false_prediction)>.1
    rho=np.outer(Omega,Omega);V=expm(1j*P)
    assert np.linalg.norm(V@rho@V.conj().T-rho)>.5


def test_joint_kernel_rectangle_bound_retains_zero_and_unnormalized_mass():
    Q=3.;S=2.;q=np.array([.5,-.2,.8]);nu=4.
    E=np.array([3.8,4.1,9.]);P=np.array([[.4,-.1,.7],[.6,-.3,.9],[5.,2.,7.]])
    weights=np.array([.1,.2,.15]);dE=2.;dp=np.ones(3)
    factors=np.sinc(S*(E-nu)/np.pi)**2*np.prod(np.sinc(Q*(P-q)/np.pi)**2,axis=1)
    p=float(weights@factors);leak=m.kernel_leakage('3','2','2',('1','1','1'))
    inside=(abs(E-nu)<=dE)&np.all(abs(P-q)<=dp,axis=1)
    lower=m.original_rectangle_lower(str(p),leak)
    assert float(lower)<=sum(weights[inside])
    assert m.original_rectangle_lower(0,leak)==0
    assert m.original_rectangle_lower('1/2','1/4','1/16')==F(3,16)


def test_full_clock_error_has_no_selector_or_reference_factor():
    H=np.diag([-1.,2.]);B=np.array([[0.,.5],[0.,0.]])
    root=np.diag([1.,np.sqrt(.75)]);times=[-.3,.1,.4];p=np.array([.5,.25,.25])
    v=np.sqrt(p);e=np.array([1.,0.,0.]);difference=e-v
    W=np.eye(3)-np.outer(difference,difference)/(1-v[0])
    exact=[];approx=[];eta=0.;r=8
    for t in times:
        U=expm(-1j*t*H);z=1j*t*H/(2*r)
        V=np.linalg.matrix_power((np.eye(2)-z)@np.linalg.inv(np.eye(2)+z),r)
        eta=max(eta,np.linalg.norm(U-V,2))
        exact.append([U.conj().T@K@U for K in (B,root)])
        approx.append([V.conj().T@K@V for K in (B,root)])
    bell=np.array([1,0,0,1])/np.sqrt(2);distance=0.
    for b in (0,1):
        for k in range(3):
            a=sum(W[k,j]*v[j]*exact[j][b] for j in range(3))
            c=sum(W[k,j]*v[j]*approx[j][b] for j in range(3))
            va=np.kron(a,np.eye(2))@bell;vc=np.kron(c,np.eye(2))@bell
            delta=np.outer(va,va.conj())-np.outer(vc,vc.conj())
            distance+=np.abs(np.linalg.eigvalsh(delta)).sum()/2
    assert distance<=2*eta+1e-14


def test_same_origin_average_differs_from_coherent_origin_superposition():
    rho=np.diag([.75,.25]);B=np.array([[0.,.5],[0.,0.]])
    pair=(B,-B)
    retained=sum(np.trace(rho@A.T@A) for A in pair)/2
    wrong=(pair[0]+pair[1])/2
    assert retained>0 and np.trace(rho@wrong.T@wrong)==0


@pytest.fixture(scope='module')
def plans():
    return {n:m.bind_member(m.contract.member(n)) for n in (0,4,8,12)}


def test_source_bound_precision_uses_actual_growing_spectrum(plans):
    for n,receipt in plans.items():
        row=receipt['result'];target=F(1,2**(n+12))
        assert F(row['inner_quadrature_operator_error_upper'])<=target
        assert F(row['one_source_evolution_operator_error_upper'])<=target
        assert row['inner_triangle_cells']%2==0
        assert row['fine_outcomes_per_origin']==2*row['selector_dimension']
        assert row['complete_fine_outcomes_with_retained_origin']==row['fine_outcomes_per_origin']*row['origin_geometry']['origin_count']
        assert row['finite_state_to_limit_convergence_modulus_supplied'] is False
    assert plans[0]['result']['active_source_aperture'] is False
    assert plans[0]['result']['fixed_regular_target_comparison_applicable'] is False
    assert plans[0]['result']['success_amplitude_comparison_upper'] is None
    assert plans[0]['result']['original_probability_comparison_upper'] is None
    assert plans[4]['result']['active_source_aperture']
    assert plans[4]['result']['fixed_regular_target_comparison_applicable']
    assert plans[8]['result']['same_source_bulk_bound']['operator_error']['negative_exponent']>100
    assert F(plans[12]['result']['original_probability_comparison_upper'])<F(plans[8]['result']['original_probability_comparison_upper'])


def test_source_and_scope_mutations_fail_recomputation(plans):
    member=m.contract.member(8);m.verify_bound(member,plans[8])
    bad=deepcopy(plans[8]);bad['result']['nonzero_joint_spectral_weight_proved']=True
    bad['result_sha256']=m.contract.digest(bad['result'])
    with pytest.raises(ValueError):m.verify_bound(member,bad)
    changed=deepcopy(member);changed['clock']['GeV_or_proper_time_conversion']=1
    with pytest.raises(ValueError):m.bind_member(changed)
    with pytest.raises(ValueError):m.bind_member(member,spatial_index=(1.0,0,0))
    with pytest.raises(ValueError):m.bind_member(member,spatial_aperture='1/3')


def test_signed_frequency_and_aliasing_do_not_create_a_finite_momentum_pvm():
    member=m.contract.member(4);row=m.bind_member(member,spatial_index=(-1,2,0),time_index=-2)['result']
    assert row['frequency_centers']=={'q_over_pi':['-1','2','0'],'nu_over_pi':'-2000'}
    assert row['finite_Fourier_label_is_exact_momentum_projection'] is False
    assert row['fixed_origin_state_replaced_or_assumed_invariant'] is False
    assert row['same_relative_observation_canonical_H_P_used'] is True
