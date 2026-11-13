from fractions import Fraction as F
import math
import numpy as np
from scipy.integrate import solve_ivp
from . import source as s

def test_derivative_bounds_and_independent_numerical_solution():
    assert s.derivative_certificate()['absolute_third_derivative_upper_exact']=='739/100'
    om,orr=F(63,200),F(23,250000)
    rows,_=s.integrate(om,orr,[0,1,10,100],F(1,64))
    x0=float(1-om-orr)
    def rhs(t,y):
        x=1/(1+np.exp(-y[0]));r=float(orr)*np.exp(t)/(float(om)+float(orr)*np.exp(t))
        return [-1-r-2*np.sqrt(x)]
    oracle=solve_ivp(rhs,[0,math.log(101)],[math.log(x0/(1-x0))],rtol=2e-12,atol=1e-13,dense_output=True)
    assert oracle.success
    for row in rows:
        value=1/(1+np.exp(-oracle.sol(math.log(1+row['redshift']))[0]))
        lo,hi=map(F,row['Omega_dark_interval']);assert float(lo)<=value<=float(hi)

def test_order_and_radiation_term_are_not_optional():
    om,orr=F(63,200),F(23,250000)
    rows,_=s.integrate(om,orr,[0,1,10],F(1,32))
    intervals=[list(map(F,row['Omega_dark_interval'])) for row in rows]
    assert all(b[1]<a[0] for a,b in zip(intervals,intervals[1:]))
    assert s.radiation_ratio(F(0),om,orr).lo>0
    import pytest
    with pytest.raises(ValueError):s.integrate(om,orr,[0,10,1])
    with pytest.raises(ValueError):s.taylor_step(F(0),s.Interval.point(0),F(1),om,orr)
