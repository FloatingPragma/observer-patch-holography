from fractions import Fraction as F
import pytest
from .source import small_control, noise_bounds, information_bounds, Q

def test_complete_basis_controls():
    for n in (3,5,7):
        row=small_control(n)
        assert abs(sum(row['all_coordinate_probabilities'])-1)<1e-12

def test_covariance_positivity_and_zero_profile():
    for n in (1,9,10**45):
        lo,hi=noise_bounds(n,F(2,3),F(1,5))
        assert 0<lo<=hi
        assert noise_bounds(n,0,1)==(0,0)

def test_no_signal_precision_erasure():
    n=10**45
    lo,hi=information_bounds(n)
    assert 0<lo<hi and hi/lo<1+F(1,10**40)
    assert lo<4*Q*Q/n<hi

@pytest.mark.parametrize('n',[0,-1,True,1.5])
def test_invalid_carrier_rejected(n):
    with pytest.raises(ValueError): information_bounds(n)
