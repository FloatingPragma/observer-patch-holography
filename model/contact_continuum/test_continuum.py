from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest
from scipy.integrate import quad

from model.contact_continuum import continuum as c
from model.contact_scattering.onshell import threshold_enclosure, scattering_eigenvalue
from model.interacting_field.lattice import symbol, SIGMA


def test_exact_negative_green_sign_margin():
    cert = c.principal_value_certificate()
    assert cert['negative_slope_magnitude_lower'] == F(103, 1620)
    assert cert['negative_slope_magnitude_lower']-F(1, 16) == F(7, 6480)
    assert cert['all_nodes'] == 8


@pytest.mark.parametrize('x',[.01,.125,.5])
@pytest.mark.parametrize('degree',[0,1,2,5,10])
def test_principal_value_recurrence_against_independent_cauchy_quadrature(x, degree):
    value, error = quad(lambda t: -t**(2*degree)/(x+t),0.,1.,weight='cauchy',wvar=x,
                        epsabs=1e-11,epsrel=1e-11)
    assert abs(value-c.pv_monomial(degree,x)) < 2e-10
    assert value <= x**(2*degree)*math.atanh(x)/x+2e-10


@pytest.mark.parametrize('speed',[F(1,8),F(1),F(8)])
@pytest.mark.parametrize('ratio',[F(1,64),F(1,4),F(1,2)])
def test_exact_eight_node_density_obeys_coarse_bound(speed,ratio):
    energy = speed*ratio
    refined = threshold_enclosure(energy,speed=speed)
    coarse = c.scattering_bound(energy,speed=speed)
    assert refined['density_upper'] < coarse['density_upper']
    assert coarse['principal_green_upper'] < 0


@pytest.mark.parametrize('u',[0.,1e-6,1.,1e8])
def test_algebraic_scattering_bound_covers_entire_admitted_green_density_envelope(u):
    e, speed = F(1,32), F(1,8)
    bounds = c.scattering_bound(e,speed=speed,interaction=F(str(u)))
    for multiplier in [1.,2.,100.]:
        a = multiplier*float(bounds['principal_green_upper'])
        rho = float(bounds['density_upper'])
        # Using pi<22/7 along with the looser density coefficient would
        # lose the correlated pi cancellation. Use the actual bound with
        # pi retained, as in the analytic proof.
        rho = min(rho,float(e**2/speed**3)*4/(15*math.pi**2))
        s = scattering_eigenvalue(a,rho,interaction=u)
        assert abs(s-1) <= float(bounds['combined_bound'])+1e-14


def test_uniform_window_schedule_and_unbounded_repulsion_sequence():
    cert = c.spacing_certificate(F(1,16),F(1,1000))
    assert cert['error_upper'] <= F(1,1000)
    for j in range(5):
        a=cert['spacing_upper']/2**j
        result=c.scattering_bound(F(1,16),speed=F(1,8)/a,interaction=a**-10)
        assert result['combined_bound'] <= F(1,1000)/2**j
    assert c.scattering_bound(F(1,32),interaction=0)['combined_bound']==0


def test_all_eight_kinetic_nodes_against_actual_field_symbol():
    p=(F(1,2),F(-2,3),F(3,4));v=F(1,8)
    for a in [F(1,2),F(1,8),F(1,32)]:
        for bits in product((0,1),repeat=3):
            k=[math.pi*b+float(a*x) for b,x in zip(bits,p)]
            limit=-float(v)*sum(((-1)**b*float(x)*s for b,x,s in zip(bits,p,SIGMA)),
                               np.zeros((2,2),complex))
            err=np.linalg.norm(symbol(k,float(v/a))-limit,2)
            assert err <= float(c.kinetic_error(p,a,velocity=v))+1e-13


def test_broadened_integral_against_complete_three_dimensional_brillouin_sum():
    e,eta=.25,.5
    n=96;k=2*np.pi*(np.arange(n)+.5)/n;s=np.sin(k)**2
    z=complex(e,eta)
    direct=np.mean(z/(z*z-4*(s[:,None,None]+s[None,:,None]+s[None,None,:])))
    assert abs(c.broadened_green(e,eta,order=n)-direct)<1e-10
    assert abs(c.broadened_green(-e,eta,order=n)+direct.conjugate())<1e-10
    assert c.broadened_green(0,eta,order=n).imag < 0


def test_existing_fiber_bound_pair_scales_out_of_fixed_energy_window():
    from model.bound_pair.fiber import PairFiber
    ratios=[]
    for speed in [.125,.5,1.]:
        pair=PairFiber((0.,0.,0.),side=8,speed=speed,interaction=10*speed)
        pole=pair.pole()['energy']
        assert pole>2*math.sqrt(3)*speed
        ratios.append(pole/speed)
    assert max(ratios)-min(ratios)<1e-12


def test_green_scaling_uses_the_same_lattice_coefficient():
    unit=c.broadened_green(.25,.5)
    for speed in [.125,.5,8.]:
        assert abs(c.broadened_green(.25*speed,.5*speed,speed=speed)-unit/speed)<1e-13


@pytest.mark.parametrize('bad', [True,.1,float('nan'),-1])
def test_invalid_exact_inputs_rejected(bad):
    with pytest.raises(ValueError):c.green_upper(bad)


def test_domain_boundaries_are_not_silently_extrapolated():
    for e,s in [(0,1),(1,0),(F(3,4),1)]:
        with pytest.raises(ValueError):c.scattering_bound(e,speed=s)
    with pytest.raises(ValueError):c.scattering_bound(F(1,32),interaction=-1)
    with pytest.raises(ValueError):c.spacing_certificate(1,3)
    with pytest.raises(ValueError):c.kinetic_error([1,2],1)
    with pytest.raises(ValueError):c.broadened_green(1,0)
