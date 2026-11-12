from fractions import Fraction as F
import numpy as np
import pytest
from . import mixing as m


def test_exact_square_root_enclosures():
    for value in (F(0),F(1,4),F(3,4),F(7,13),F(1,2)**100):
        lo,hi=m.sqrt_box(value)
        assert lo*lo<=value<=hi*hi and hi-lo<=F(1,2**180)


def test_independent_two_state_gap():
    bound=m.eigensystem_enclosure(np.array([[.75,.25],[.25,.75]]),F(0))
    a,b=map(F,bound['gap_interval'])
    assert a<=F(1,2)<=b and b-a<F(1,10**10)
    for bad,error in [(np.zeros((2,3)),F(0)),(np.array([[1.,.2],[.3,1.]]),F(0)),
                      (np.eye(2)*np.nan,F(0)),(np.eye(2),F(-1))]:
        with pytest.raises(ValueError):m.eigensystem_enclosure(bad,error)


def test_state_lag_uses_original_error_without_renormalization():
    target=[F(1,3),F(2,3)];p=[.4,.6];error=F(1,1000)
    lo,hi=map(F,m.lag_interval(p,target,error))
    assert lo<F(1,15)<hi and hi-lo==error


def test_attempt_budget_and_all_components():
    assert len(m.bath.source()['components'])==26
    assert sum(len(c)==1 for c in m.bath.source()['components'])==12
    assert m.mixing_steps(F(1,2),F(1,2),F(1,1000))>=18
