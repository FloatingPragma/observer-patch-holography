"""Independent source resolvents, coarea calculus and moving-scattering checks."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import math

import mpmath as mp
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from ..band_edge_contact.limit import scaled_green_diagnostic
from . import bounds as b, diagnostic as d


def test_exact_majorants_and_common_moving_witness():
    proof = b.coarea_certificate()
    assert proof['angular_amplitude_derivative_upper'] == 48
    assert proof['scaled_green_energy_error_coefficient'] == F(25385, 96)
    joined = b.joined_certificate(F(1, 16384), 1, 1, 1, 1)
    assert joined['scattering']['uniform_scattering_operator_error_upper'] == F(737, 3584) < F(1, 4)
    assert joined['bound_clock_operator_error_upper'] < F(28, 1000)
    assert joined['scattering_free_clock_operator_error_upper'] == F(1, 1073741824)
    assert not joined['relativistic_or_many_body_QFT_constructed']


def test_actual_four_spin_source_resolvent_interior_and_exterior():
    pauli = (np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]), np.diag([1.,-1.]))
    singlet = np.array([0, 1, -1, 0])/math.sqrt(2)
    rng = np.random.default_rng(400113)
    for _ in range(30):
        q, k = rng.uniform(-math.pi, math.pi, 3), rng.uniform(-.25,.25,3)
        speed = .73
        h = sum(-speed*(np.sin(q[j]+k[j]/2)*np.kron(pauli[j],np.eye(2))+
                        np.sin(-q[j]+k[j]/2)*np.kron(np.eye(2),pauli[j])) for j in range(3))
        for energy in (2.1+.037j, 3.7+.019j):
            actual = np.vdot(singlet,np.linalg.solve(energy*np.eye(4)-h,singlet))
            assert d.kernel(q,k,energy,speed) == pytest.approx(actual,abs=5e-13)


def test_all_eight_edge_charts_and_coarea_source_normalization():
    k, r = np.array([.23,-.17,.11]), np.array([.17,-.32,.41])
    t, speed = 1.71, .8
    u = np.sin(r)
    e2 = np.sum(np.sin(k/2)**2*u*u)
    de = np.sum(np.cos(k/2)*np.sin(k/2)*u*np.sqrt(1-u*u))
    R2 = np.sum(np.cos(k/2)**2)
    shape = np.sum(np.cos(k)*u*u)+de*de/t**2
    expected = t/(2*speed)*(1-e2/t**2)/(shape+t*t-R2)
    for signs in product((-1,1),repeat=3):
        q=np.asarray(signs)*math.pi/2+r
        assert d.kernel(q,k,2*speed*t,speed) == pytest.approx(expected,abs=2e-12)
    assert 8*F(1,8) == 1


def test_uniform_coarea_jacobian_derivatives_against_high_precision_calculus():
    rng=np.random.default_rng(1115)
    with mp.workdps(70):
        for _ in range(60):
            k=rng.uniform(-.25,.25,3)
            omega=rng.normal(size=3); omega/=np.linalg.norm(omega)
            x=float(rng.uniform(0,.25))
            h=float(rng.uniform(0,1/128))
            t=math.sqrt(float(np.sum(np.cos(k/2)**2)))-h/(2*math.sqrt(3))
            row=d.angular_amplitude(x,omega,k,t)
            assert 31/32-1e-14<=row['f']<=1+1e-14
            assert abs(row['df_dx'])<=max(abs(k))**2+1e-14
            assert abs(row['d2f_dx2'])<=2*max(abs(k))**2+1e-14
            assert row['D']>=61/64
            assert abs(row['dA_dy2'])<48
            kk,ww,tt=tuple(map(mp.mpf,k)),tuple(map(mp.mpf,omega)),mp.mpf(t)
            def f(xx):
                w=sum(mp.cos(ki/2)*mp.sin(ki/2)*wi*mp.sqrt(1-xx*wi*wi) for ki,wi in zip(kk,ww))
                return sum(mp.cos(ki)*wi*wi for ki,wi in zip(kk,ww))+w*w/(tt*tt)
            def A(xx):
                J=mp.fprod((1-xx*wi*wi)**(-mp.mpf('.5')) for wi in ww)
                N=1-xx*sum(mp.sin(ki/2)**2*wi*wi for ki,wi in zip(kk,ww))/(tt*tt)
                return J*N/(mp.sqrt(f(xx))*(f(xx)+xx*mp.diff(f,xx)))
            expected=mp.diff(A,mp.mpf(x))/(f(x)+x*mp.diff(f,mp.mpf(x)))
            assert row['dA_dy2']==pytest.approx(float(expected),abs=4e-14)


def test_PV_subtraction_against_independent_weighted_cauchy_integral():
    k=np.array([.12,-.07,.09]); vector=np.array([1.,.37,-.51])
    R=math.sqrt(float(np.sum(np.cos(k/2)**2)))
    lam=.004; t=math.sqrt(R*R-lam)
    root=brentq(lambda l:d.radial(l,vector,k,t)[0]-lam,0,.25,xtol=1e-15)
    _,n0,fp=d.radial(root,vector,k,t)
    def smooth(l):
        f,n,_=d.radial(l,vector,k,t)
        if abs(l-root)<1e-10:return root*root*n0/fp
        return l*l*n*(l-root)/(f-lam)
    real,error=quad(smooth,0,math.pi/2,weight='cauchy',wvar=root,epsabs=1e-11,epsrel=1e-11)
    assert error<1e-9
    expected=complex(real,-math.pi*root*root*n0/fp)
    assert d.radial_pv(vector,k,t,lam,order=48)==pytest.approx(expected,abs=2e-11)


def test_zero_total_momentum_agrees_with_independent_closed_axis_integration():
    actual=d.scaled_green(1/128,[0,0,0],1,order=24)
    expected=scaled_green_diagnostic(1/128,-1)
    assert actual==pytest.approx(expected,abs=2e-8)


def test_moving_on_shell_scattering_converges_with_correct_time_orientation():
    errors=[]
    for a in (1/64,1/128):
        value=d.scattering(a,1,[1,.5,-.75],1,order=16)
        assert abs(value['reflected_phase'])==pytest.approx(1,abs=2e-15)
        assert value['target']==-1j
        assert value['absolute_error']<float(b.green_certificate(F(a),1,1)['scaled_green_error_upper'])
        assert abs(value['reflected_phase'].conjugate()-value['target'])>1.9
        errors.append(value['absolute_error'])
    assert errors[1]<.6*errors[0]


def test_critical_countertuning_and_nonidentity_are_required():
    for a in (F(1,100),F(1,1000)):
        denominator=1/float(a)+1j
        source=denominator.conjugate()/denominator
        assert abs(source-1)<=2*float(a)
    joined=b.joined_certificate(F(1,16384),1,1,1,1)
    error=joined['scattering']['uniform_scattering_operator_error_upper']
    assert F(5,4)-error>1


@pytest.mark.parametrize('g,p,s,eps',[(1,1,1,F(1,100)),(F(1,2),3,4,F(1,2)),(2,0,F(1,3),1)])
def test_rational_schedule_admits_arbitrary_smaller_spacings(g,p,s,eps):
    result=b.spacing_certificate(g,p,s,eps)
    for factor in (1,2,7):
        smaller=b.joined_certificate(result['spacing_upper']/factor,g,p,s,1)
        assert smaller['scattering']['uniform_scattering_operator_error_upper']<=eps


@pytest.mark.parametrize('call',[lambda:b.green_certificate(F(1,2),1,1),
    lambda:b.scattering_certificate(F(1,128),1,1,1),
    lambda:b.scattering_certificate(F(1,16384),0,1,1),
    lambda:b.green_certificate(True,1,1),lambda:b.green_certificate(F(1,100),1,False),
    lambda:b.spacing_certificate(1,1,1,3)])
def test_outside_domain_and_boolean_parameters_fail(call):
    with pytest.raises(ValueError):call()


def test_receipt_scope_mutations_fail():
    from . import packet
    result=packet.build()
    bad=deepcopy(result);bad['scope']['relativistic_many_body_QFT_constructed']=True
    with pytest.raises(ValueError):packet.validate(bad,result)
    bad=deepcopy(result);bad['coarea']['all_edge_valleys_retained']=1
    with pytest.raises(ValueError):packet.validate(bad,result)
    bad=deepcopy(result);bad['scope']['particle_number_sector']=True
    with pytest.raises(ValueError):packet.validate(bad,result)
