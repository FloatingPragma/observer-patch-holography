from fractions import Fraction as F
import math
import pytest
from .source import h2,log_interval,event_primitive,event_omega,event_compare,run


def test_horizon_branches_and_normalization():
    for branch in ('fixed_capacity','comoving_area','hubble_area_conserved'):
        assert h2(0,branch,67,F(3,10))==67**2
    for z in (F(1,4),F(1),F(3)):
        assert h2(z,'fixed_capacity',67,F(3,10))<h2(z,'comoving_area',67,F(3,10))<h2(z,'hubble_area_conserved',67,F(3,10))
    with pytest.raises(ValueError):h2(1,'unknown',67,F(3,10))


def test_rational_logarithms_and_primitive_against_independent_arithmetic():
    import mpmath as mp
    with mp.workdps(100):
        m=lambda x:mp.mpf(x.numerator)/x.denominator
        for x in (F(1,10000),F(1),F(1,2),F(7,3),F(10000)):
            lo,hi=log_interval(x);assert m(lo)<=mp.log(m(x))<=m(hi)
        for x in (F(1,10),F(1,2),F(9,10)):
            lo,hi=event_primitive(x);v=m(x)
            direct=2*mp.log(v)+mp.log(1+v)-mp.log(1-v)/3-8*mp.log(1+2*v)/3
            assert m(lo)<=direct<=m(hi)


def test_implicit_event_solution_against_independent_ode():
    from scipy.integrate import solve_ivp
    om=F('0.315')
    for z in (F(1,10),F(1),F(3)):
        lo,hi=event_omega(z,om)
        sol=solve_ivp(lambda t,y:y*(1-y)*(1+2*y**.5),(0,-math.log1p(float(z))),[float(1-om)],rtol=2e-12,atol=1e-14)
        assert float(hi-lo)<1e-15
        assert abs(sol.y[0,-1]-float((lo+hi)/2))<2e-11
    far=event_omega(F(10**12),om)
    assert 0<far[0]<far[1]<F(1,10**10)
    for extreme_om in (F(1,10**60),1-F(1,10**60)):
        lo,hi=event_omega(F(1),extreme_om)
        assert 0<lo<hi<1


def test_all_rows_covariance_and_calibration_boundaries():
    r=run();assert len(r['comparisons'])==8
    assert all(x['rows']==15 and len(x['records'])==15 for x in r['comparisons'])
    assert not r['same_data_alternatives_are_independent']
    assert not r['physical_horizon_or_Lambda_amplitude_derived']
    c={x['branch']:x for x in r['comparisons'] if x['table']=='HzTable_MM_BC03.dat'}
    assert F(c['hubble_area_conserved']['chi_square_interval'][0])>70
    assert F(c['future_event_area_gamma_one']['chi_square_interval'][1])<7
