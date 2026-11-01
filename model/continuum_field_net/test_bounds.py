from fractions import Fraction as F
import math

import pytest

from . import bounds


def test_uniform_time_and_packet_bounds_do_not_contain_inverse_spacing():
    assert bounds.time_lipschitz(3,2,4,5)==22
    value=bounds.window_response(1,2,3,2,1)
    assert value['inverse_norm_upper']==F(35,4)
    assert value['negative_imaginary_lower']==F(1,81)
    assert value['full_source_second_moment_upper']==49
    assert bounds.valley_overlap(F(1,8),2,3,4,5)==F(11,8)


def test_central_difference_remainder_on_independent_plane_waves():
    for a in (.5,.2,.05):
        for p in (-2.3,-.2,.7,3.):
            error=abs(math.sin(a*p)/a-p)
            assert error<=a*a*abs(p)**3/6+1e-15
    assert bounds.free_dynamics_error(F(1,10),2,6,1,3,F(1,4))==F(14,25)
    assert bounds.free_dynamics_error(F(1,20),2,6,1,3,0)==F(3,200)


def test_scaled_cone_bound_dominates_original_full_CAR_estimate():
    # Compare the uncapped exponential exponents in logarithms, so no
    # underflow can create a false zero at small spacing.
    for a in (F(1,8),F(1,16),F(1,32)):
        r,R,t,v=F(1),F(10),F(1,100),F(1)
        value=bounds.cone_bound(a,r,R,t,v)
        ri,Ro=value['inner_lattice_radius'],value['outer_lattice_radius']
        speed=float(v/a);sites=(2*ri+1)**3
        original_log=math.log(12*speed*sites*(2*Ro+1)**2*float(t))-(Ro-ri)*math.log(2)+48*speed*float(t)
        majorant_log=math.log(float(value['prefactor']))-6*math.log(float(a))-float(value['exponential_margin']/a)
        assert original_log<=majorant_log+1e-12
        assert Ro-ri>=float((R-r)/a-4)


def test_cone_remainder_goes_to_zero_at_explicit_linear_rate():
    first=bounds.cone_bound(F(1,32),1,100,F(1,100))
    second=bounds.cone_bound(F(1,64),1,100,F(1,100))
    assert 0<second['linear_in_spacing_upper']<2
    assert first['linear_in_spacing_upper']==2*second['linear_in_spacing_upper']
    assert bounds.cone_bound(F(1,32),1,10,0)['linear_in_spacing_upper']==0


@pytest.mark.parametrize('call',[
    lambda:bounds.cone_bound(F(1,2),1,2,0),
    lambda:bounds.cone_bound(F(1,32),1,2,1),
    lambda:bounds.window_response(2,1,1),
    lambda:bounds.time_lipschitz(-1),
    lambda:bounds.free_dynamics_error(.1,1,1),
])
def test_unproved_domains_and_nonexact_input_rejected(call):
    with pytest.raises(ValueError):call()
