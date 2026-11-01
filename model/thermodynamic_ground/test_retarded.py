from fractions import Fraction as F
import math

import pytest

from model.thermodynamic_ground import retarded as r


def test_rational_exponential_and_retarded_tails_enclose_independent_values():
    for value in (F(0),F(1,10),F(1),F(8),F(32)):
        assert float(r.exp_negative_upper(value))>=math.exp(-float(value))*(1-2e-15)
    assert r.exp_negative_upper(0)==1
    assert r.tail_bound(1,32)<r.tail_bound(1,16)<r.tail_bound(1,8)<F(1,1000)


def test_exact_single_spectral_line_obeys_tail_and_derivative_bounds():
    # C(t)=exp(-i omega t), G(z)=1/(z-omega), a legitimate CAR line.
    eta,time,energy=F(1,2),F(3),F(2)
    z=1+complex(0,float(eta))
    exact_tail=math.exp(-float(eta*time))/abs(z-float(energy))
    assert exact_tail<float(r.tail_bound(eta,time))
    for n in range(5):
        assert math.factorial(n)/abs(z-float(energy))**(n+1)<=float(r.derivative_bound(eta,n))


def test_state_discrepancy_is_counted_and_never_inferred_from_boundary_bound():
    args=(64,16,0,1,F(1,8),F(1))
    first=r.comparison_bound(*args,state_discrepancy=F(1,10))
    second=r.comparison_bound(*args,state_discrepancy=F(2,10))
    assert second-first==F(1,80)
    with pytest.raises(TypeError):r.comparison_bound(*args)


@pytest.mark.parametrize('call',[
    lambda:r.tail_bound(0,1),lambda:r.tail_bound(-1,1),lambda:r.tail_bound(1,-1),
    lambda:r.derivative_bound(1,True),lambda:r.derivative_bound(1,-1),
    lambda:r.exp_negative_upper(1.0),lambda:r.exp_negative_upper(1,0),
])
def test_no_real_axis_or_inexact_domain_admission(call):
    with pytest.raises(ValueError):call()
