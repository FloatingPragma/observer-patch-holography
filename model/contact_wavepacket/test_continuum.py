from fractions import Fraction as F
import math
import numpy as np
import pytest
from scipy.integrate import quad
from model.contact_wavepacket import continuum as c
from model.contact_continuum.continuum import broadened_green
from model.contact_scattering.onshell import threshold_enclosure, pi_interval


def test_exact_slope_margin_and_density_majorant():
    cert=c.slope_certificate()
    assert cert['positive_derivative_correction_upper']==F(6112,535815)
    assert cert['negative_slope_magnitude_lower']-F(1,16)==F(80813,8573040)
    for x in [F(1,64),F(1,8),F(1,4)]:
        density=threshold_enclosure(x,speed=F(1))
        # Keep correlated pi factors, then compare the actual exact enclosure.
        p_lo,p_hi=pi_interval()
        assert p_hi*density['density_upper'] < 16*x*x/(63*p_lo)


@pytest.mark.parametrize('n',[0,1,2,5,10])
@pytest.mark.parametrize('x',[.0625,.125,.25])
def test_derivative_recurrence_against_independent_principal_value_quadrature(n,x):
    def value(z):
        integral,_=quad(lambda t:-t**(2*n)/(z+t),0.,1.,weight='cauchy',wvar=z,
                        epsabs=1e-11,epsrel=1e-11)
        return z**3*integral
    h=1e-5
    derivative=(value(x+h)-value(x-h))/(2*h)
    exact_formula=c.pv_correction_derivative_monomial(n,x)
    assert abs(derivative-exact_formula)<2e-9
    assert exact_formula <= (2*n+3)*x**(2*n+2)/(1-x*x)+1e-15


@pytest.mark.parametrize('s',[.1,.5,.9,1.1,1.5,1.9])
def test_two_axis_marginal_against_original_arcsine_convolution(s):
    x=min(s,2-s)
    actual,_=quad(lambda t:2/math.pi**2/math.sqrt(
        (1-x*math.sin(t)**2)*(1-x*math.cos(t)**2)),0.,math.pi/2)
    assert abs(actual-c.two_sine_density(s))<2e-13


def test_marginal_resolvent_matches_original_torus_source():
    z=.125+.3j
    def f(s):
        value=z*c.two_sine_density(s)/(np.sqrt(z*z-4*s)*np.sqrt(z*z-4*(s+1)))
        return value
    value=complex(quad(lambda s:f(s).real,0,2,points=[1],epsabs=1e-11)[0],
                  quad(lambda s:f(s).imag,0,2,points=[1],epsabs=1e-11)[0])
    assert abs(value-broadened_green(z.real,z.imag,order=256))<1e-10


@pytest.mark.parametrize('x',[.0625,.125,.25])
def test_actual_on_shell_green_is_odd_and_has_stronger_negative_slope(x):
    a=c.principal_green_diagnostic(x)
    assert a < -x/16
    assert c.principal_green_diagnostic(-x)==-a
    h=1e-5
    derivative=(c.principal_green_diagnostic(x-h)-c.principal_green_diagnostic(x-2*h))/h
    assert derivative < -1/16
    for speed in [.125,4.]:
        assert abs(c.principal_green_diagnostic(x*speed,speed=speed)-a/speed)<1e-12


@pytest.mark.parametrize('speed',[4.,16.,64.])
@pytest.mark.parametrize('coupling_ratio',[-100.,-8.,-.1,.1,8.,100.])
def test_integrated_bound_for_independent_monotone_green_envelope(speed,coupling_ratio):
    # Synthetic admissible envelope; not a claim this is the exact lattice g.
    u=coupling_ratio*speed**2
    root=-4/coupling_ratio
    def value(e):
        if e==0:return 0.
        s=c.signed_scattering(-e/(4*speed**2),e*e/(4*math.pi**2*speed**3),u)
        return abs(s-1)**2
    points=[0.]+([root] if -1<root<1 else [])
    integral,_=quad(value,-1,1,points=sorted(points),epsabs=1e-10)
    bound=c.window_certificate(1,speed=F(str(speed)))
    assert integral<=float(bound['integrated_squared_scattering_bound'])


def test_actual_tuned_resonance_is_not_operator_norm_convergence():
    rows=c.resonance_diagnostic()
    assert all(row['repulsive_U_diagnostic']>0 for row in rows)
    assert all(abs(row['peak_scattering_distance']-2)<1e-14 for row in rows)
    off=[row['off_resonance_distance'] for row in rows]
    assert all(a>3.9*b for a,b in zip(off,off[1:]))


def test_arbitrarily_localized_packets_do_not_obey_a_fixed_density_cap():
    # A narrowing normalized packet can follow a resonance; D grows with c.
    errors=[]
    for speed in [16.,64.,256.]:
        center=-.5;width=center**2/(math.pi*speed);u=-4*speed**2/center
        def value(z):
            e=center+width*z
            s=c.signed_scattering(-e/(4*speed**2),e*e/(4*math.pi**2*speed**3),u)
            return abs(s-1)**2/2
        errors.append(quad(value,-1,1,epsabs=1e-10)[0])
    assert min(errors)>3.


def test_uniform_rational_packet_schedule():
    e=F(1,8);d=4;v=F(1,8);tol=F(1,1000)
    cert=c.spacing_certificate(e,tol,velocity=v,packet_density_cap=d)
    assert cert['packet_squared_error_bound']<=tol
    for j in range(5):
        a=cert['spacing_upper']/2**j
        row=c.window_certificate(e,speed=v/a,packet_density_cap=d,threshold=F(1,2))
        assert row['packet_squared_error_bound']<=tol/2**j
        assert row['exceptional_energy_measure_upper']<=2*e


def test_signed_extension_conjugation_and_zero_contact():
    for u in [-1e6,-1.,1.,1e6]:
        s=c.signed_scattering(-.1,.01,u)
        assert abs(abs(s)-1)<2e-15
        assert abs(c.signed_scattering(.1,.01,-u)-s.conjugate())<1e-15
    assert c.signed_scattering(.1,.01,0)==1


@pytest.mark.parametrize('bad',[True,.1,float('nan'),-1])
def test_invalid_exact_inputs(bad):
    with pytest.raises(ValueError):c.window_certificate(bad,speed=1)


def test_domains_never_silently_enlarge():
    with pytest.raises(ValueError):c.window_certificate(1,speed=1)
    with pytest.raises(ValueError):c.window_certificate(0,speed=1)
    with pytest.raises(ValueError):c.window_certificate(1,speed=4,threshold=3)
    with pytest.raises(ValueError):c.spacing_certificate(1,1,velocity=1,packet_density_cap=0)
    with pytest.raises(ValueError):c.principal_green_diagnostic(.3)
    with pytest.raises(ValueError):c.signed_scattering(0,0,1)
