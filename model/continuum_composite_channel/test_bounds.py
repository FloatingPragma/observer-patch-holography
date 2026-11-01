from fractions import Fraction as F
import pytest
from . import bounds


def test_dyson_remainder_against_integer_product():
    for m in (1,2,3,7):
        for n in (1,2,5,10):
            expected=F(1)
            for k in range(n):expected*=F(m+2*k,8*(k+1))
            assert bounds.dyson_remainder(m,n,1,F(1,8))==expected


def test_geometric_short_time_error_bound():
    for m in (1,3,8,40):
        for ell in (1,4,12):
            assert bounds.dyson_remainder(m,m+ell,1,F(1,8))<=F(1,2**ell)


def test_composed_degree_and_error_recurrences():
    m=3;ell=9;K=5;degree=m;norm=F(1);error=F(0)
    for _ in range(K):
        degree=3*degree+2*ell-2
        error+=F(1,2**ell)*norm
        norm*=1+F(1,2**ell)
    answer=bounds.approximation(m,1,F(K,8),ell)
    assert answer=={'steps':K,'degree':degree,'norm_upper':norm,'error_upper':error}


def test_fixed_time_centrality_with_chosen_precision():
    previous=None
    for power in (8,12,16,20):
        a=F(1,2**power);ell=bounds.precision_for_spacing(a,8)
        value=bounds.commutator_upper(3,1,1,ell,2*a**3,3)
        if previous is not None:assert value<previous/1000
        previous=value
    assert previous<F(1,10**10)


def test_ground_bound_retains_nonzero_local_variance_and_decays_for_averages():
    assert bounds.mixed_defect_ground_norm_squared(0,1,1)==4
    assert bounds.mixed_defect_ground_norm_squared(0,1,F(1,1000))==F(1,250)
    values=[bounds.mixed_defect_ground_norm_squared(F(1,2**(4*k)),1,F(1,2**(3*k))) for k in (4,8,12)]
    assert values[0]>values[1]>values[2]>0


def test_smearing_uses_physical_a_cubed():
    result=bounds.smearing_bounds(F(1,8),2,3,4,5)
    assert result['quadratic_coefficient_sup']==F(5,128)
    assert result['I4_upper']==F(25,128)
    assert result['l2_product']==6


def test_invalid_or_inexact_input_is_rejected():
    for value in (True,1.0,-1):
        with pytest.raises(ValueError):bounds.rational(value)
    with pytest.raises(ValueError):bounds.dyson_remainder(0,1,1,1)
    with pytest.raises(ValueError):bounds.approximation(1,1,1,0)
    with pytest.raises(ValueError):bounds.mixed_defect_ground_norm_squared(3,1,1)


def test_free_dynamics_does_not_inflate_degree():
    assert bounds.approximation(3,0,100,1)=={'steps':0,'degree':3,'norm_upper':F(1),'error_upper':F(0)}
