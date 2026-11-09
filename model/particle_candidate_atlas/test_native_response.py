"""Original-probability, unitary-distance and common-kinetic controls."""
import json

import numpy as np
import pytest
import scipy.linalg as la
import sympy as s

from . import native_response as n
from . import memory_diffusion as source


def clean(a):
    return s.Matrix(a).applyfunc(lambda x:s.factor(s.radsimp(s.expand(x))))


def matrices(p,q):
    a,b=p*(1-q),p*q
    t=s.Matrix([[1-a,1-b,1-a],[a,0,0],[0,b,a]])
    e=s.Matrix([[b,a,b],[a,0,0],[0,b,a]])
    eq=s.eye(3)-t;eq[2,:]=s.ones(1,3)
    pi=eq.inv()*s.Matrix([0,0,1])
    return t,e,pi


@pytest.mark.parametrize('p,q',[(s.Rational(1,12),s.Rational(1,6)),
                               (s.Rational(1,24),s.Rational(1,12)),
                               (s.Rational(2,3),s.Rational(3,4))])
def test_original_weights_from_independent_stationary_linear_system(p,q):
    t,e,pi=matrices(p,q)
    w=e*pi
    assert tuple(n.c_duty_weights(p,q)) == (1-sum(w),*w)
    assert t*pi == pi and sum(pi)==1


def test_native_weights_equal_all_raw_records_with_disabled_probability():
    row=source.exact_response('C');w=row['weight']
    assert tuple(n.c_duty_weights(s.Rational(1,12),s.Rational(1,6))) == (
        s.Rational(11,12),w[2],w[3],w[9]+w[10])
    assert row['raw_records_over_all_old_inputs']==1536


@pytest.mark.parametrize('p,q,tangent',[(s.Rational(1,12),s.Rational(1,6),s.Rational(1,10)),
                                       (s.Rational(1,24),s.Rational(1,12),s.Rational(-1,20))])
def test_complete_finite_outputs_and_lower_bound(p,q,tangent):
    value=n.finite_control(p,q,tangent)
    histories=[s.Matrix(x).applyfunc(s.sympify) for x in value['subnormalized_history_densities']]
    weights=list(map(s.sympify,value['original_weights']))
    assert [s.simplify(s.trace(x)) for x in histories] == weights
    out=clean(sum(histories,s.zeros(3)))
    loss=s.simplify(1-s.trace(out*out))
    floor=s.sympify(value['proved_purity_loss_floor'])
    assert loss>=floor>0 and s.trace(out)==1
    # The four pair contributions are all kept; the proof isolates one pair.
    states=[clean(x/w) for x,w in zip(histories,weights,strict=True)]
    pair_loss=sum(weights[i]*weights[j]*(1-s.trace(states[i]*states[j]))
                  for i in range(4) for j in range(4))
    assert s.simplify(pair_loss-loss)==0
    isolated=2*weights[0]*weights[2]*(1-s.trace(states[0]*states[2]))
    assert s.simplify(isolated-floor)==0


def test_lower_bound_vanishes_quadratically_without_asserting_a_finite_limit():
    p,q=s.Rational(1,12),s.Rational(1,6)
    coarse=n.purity_floor(p,q,s.Rational(1,10))
    fine=n.purity_floor(p,q,s.Rational(1,20))
    assert coarse/fine == s.Rational(401,101)
    assert coarse>fine>0


@pytest.mark.parametrize('args',[(0,s.Rational(1,6),1),(1,s.Rational(1,6),1),
                                (s.Rational(1,12),0,1),(s.Rational(1,12),1,1),
                                (s.Rational(1,12),s.Rational(1,6),0),
                                (.1,s.Rational(1,6),1),
                                (True,s.Rational(1,6),1)])
def test_invalid_rate_or_clock_domain(args):
    with pytest.raises(ValueError):n.purity_floor(*args)


def test_independent_noncommuting_hamiltonian_trace_identity_and_zero_control():
    x=s.Matrix([[0,1],[1,0]]);z=s.diag(1,-1)
    y=s.Matrix([[1,1,0],[1,3,1],[0,1,2]])
    a=s.I*s.Matrix([[0,1,2],[-1,0,3],[-2,-3,0]])
    assert y*a-a*y != s.zeros(3)
    h,k=s.Rational(1,5),s.Rational(1,7)
    for potential in (a,s.zeros(3)):
        b=s.kronecker_product(x,y)+h*s.kronecker_product(s.eye(2),potential)
        full=k*s.kronecker_product(z,s.eye(3))+b
        defect=s.expand(s.trace(full**4)-6*k**4-2*k*k*s.trace(b*b)-s.trace(b**4))
        assert defect==8*h*h*k*k*s.trace(potential*potential)


def test_native_ordering_term_has_strictly_positive_trace_defect():
    value=n.hamiltonian_trace_control()
    assert s.sympify(value['trace_A_squared'])>0
    assert s.sympify(value['same_speed_mass_shell_trace_defect'])>0
    assert value['ordering_only_dynamics_is_the_complete_source'] is False


@pytest.mark.parametrize('tangent',[s.Rational(1,10),s.Integer(1)])
def test_full_CAR_pairing_against_independent_two_pair_tensor_product(tangent):
    value=n.four_mode_pairing(tangent)
    c=(1-tangent*tangent)/(1+tangent*tangent);v=2*tangent/(1+tangent*tangent)
    plus=s.eye(4);plus[0,0]=plus[3,3]=c;plus[0,3]=plus[3,0]=-s.I*v
    minus=plus.conjugate()
    independent=s.kronecker_product(minus,plus)
    assert clean(value['unitary']-independent)==s.zeros(16)
    assert value['generator'].eigenvals()=={-2:1,-1:4,0:6,1:4,2:1}
    assert value['vacuum_amplitude']==c*c


def test_full_CAR_source_vacuum_floor_retains_all_records_and_clock_boundary():
    value=n.full_car_vacuum_control()
    weights=list(map(s.sympify,value['original_weights']))
    overlap=s.sympify(value['disabled_P3_overlap_squared'])
    floor=s.sympify(value['full_CAR_purity_loss_floor'])
    assert sum(weights)==1 and floor==2*weights[0]*weights[2]*(1-overlap)>0
    assert value['other_records_retained'] is True
    assert value['full_ninety_mode_four_record_density_replayed'] is False
    assert value['physical_vacuum_selection_claimed'] is False
    with pytest.raises(ValueError):n.full_car_vacuum_control(tangent=0)


def test_common_H0_cancels_out_of_actual_second_generator_coefficient():
    # Direct slow-block Taylor expansion of the 108-dimensional transfer,
    # independent of the algebraic H0-cancellation proof.
    row=source.exact_response('C')
    p=[np.asarray(row['projectors'][i].evalf(30),complex) for i in (2,3,10)]
    t,e,pi=matrices(s.Rational(1,12),s.Rational(1,6))
    t,e,pi=(np.asarray(a,complex) for a in (t,e,pi))
    x=np.array([[0,1],[1,0]],complex);z=np.diag([1.,-1.])
    identity=np.eye(6)
    def liouvillian(h):return -1j*(np.kron(identity,h)-np.kron(h.T,identity))
    l0=liouvillian(np.kron(z,np.eye(3))/7)
    ls=[liouvillian(np.kron(x,a)) for a in p]
    injection=np.kron(pi,np.eye(36));summation=np.kron(np.ones((1,3)),np.eye(36))
    m0=np.kron(t,np.eye(36));projection=injection@summation
    fast=la.inv(np.eye(108)-m0+projection)-projection
    m1=np.block([[t[i,j]*l0+e[i,j]*ls[i] for j in range(3)] for i in range(3)])
    m2=np.block([[ (t[i,j]*(l0@l0)+e[i,j]*(l0@ls[i]+ls[i]@l0+ls[i]@ls[i]))/2
                   for j in range(3)] for i in range(3)])
    g0=summation@m1@injection
    g1=summation@m2@injection+summation@m1@fast@m1@injection-g0@g0/2
    stationary=pi@np.ones((1,3));group=la.inv(np.eye(3)-t+stationary)-stationary
    w=(e@pi)[:,0];k=np.diag(w)/2-np.outer(w,w)/2+e@group@np.diag(w)
    expected=sum((k[i,j]*ls[i]@ls[j] for i in range(3) for j in range(3)),np.zeros((36,36),complex))
    assert la.norm(g0-l0-sum((w[i]*ls[i] for i in range(3)),np.zeros((36,36),complex)),2)<2e-13
    assert la.norm(g1-expected,2)<2e-12


def test_json_certificate_has_no_physical_promotion_and_is_detached():
    value=n.certificate();expected=json.dumps(value,sort_keys=True)
    assert value['physical_mass_assigned'] is None
    assert value['physical_selfenergy_computed'] is False
    assert value['unitary_dilation_or_complete_archive_excluded'] is False
    value['native_finite_step']['original_weights'][0]='forged'
    assert json.dumps(n.certificate(),sort_keys=True)==expected
