from fractions import Fraction as F
import mpmath as mp
import pytest
from . import filters as f, coefficients as c, envelope


def test_exact_compact_filter_and_plateau():
    a,b,d=F(1),F(2),F(1,4)
    assert [f.energy_filter(x,a,b,d) for x in (0,F(3,4),F(7,8),1,F(3,2),2,F(17,8),F(9,4),3)] == [0,0,F(1,2),1,1,1,F(1,2),0,0]
    for n in range(-20,80):
        assert 0 <= f.energy_filter(F(n,20),a,b,d) <= 1


def test_goldstone_filter_retains_the_proved_band():
    m=F(1,1000);x=f.goldstone_filter(m)
    assert x['kernel_L1_upper']==28
    assert x['support_lower']==m/16
    assert x['support_upper'] < 64/m**2
    assert f.energy_filter(m/8,m/8,32/m**2,m/16)==1
    assert f.energy_filter(32/m**2,m/8,32/m**2,m/16)==1
    assert x['kernel_first_absolute_moment_upper']==F(256000,3)


@pytest.mark.parametrize('x',[False,True,0.5,'1'])
def test_rejects_inexact_or_boolean_inputs(x):
    with pytest.raises(ValueError):f.energy_filter(x,1,2,F(1,4))


def mpf(x):return mp.mpf(x.numerator)/x.denominator


@pytest.mark.parametrize('t',[F(0),F(1,7),F(-19,3),F(10**8,3)])
def test_coefficient_rectangle_against_independent_high_precision(t):
    mp.mp.dps=100;a,b,d=F(1),F(2),F(1,4)
    real,imag=c.kernel_rectangle(t,a,b,d,70)
    x=mpf(t);length=mpf(b-a+d)
    z=length/(2*mp.pi) if not t else mp.exp(-1j*mpf((a+b)/2)*x)*mp.sin(length*x/2)/(mp.pi*x)*(mp.sin(mpf(d)*x/4)/(mpf(d)*x/4))**2
    assert mpf(real[0]) <= mp.re(z) <= mpf(real[1])
    assert mpf(imag[0]) <= mp.im(z) <= mpf(imag[1])


def test_constructed_coefficient_has_requested_error_and_phase():
    r=c.midpoint_coefficient(F(1,3),F(1,16),1,2,F(1,4),F(1,10**18))
    assert r['complex_error_upper'] <= F(1,10**18)
    assert r['imag']<0  # opposite sign would select negative excitation energy


def test_explicit_envelope_has_no_sampling_alias_and_finite_certified_tails():
    value=envelope.certificate()
    assert F(value['Fourier_support_radius_lower'])>F(1,2)
    assert F(value['Fourier_support_radius_squared_upper'])<1
    assert F(value['sinc16_integral_over_pi'])>0
    assert 0<F(value['sampling_LCU_bound_Lambda'])<F(value['weighted_sampling_LCU_bound_Lambda4'])
