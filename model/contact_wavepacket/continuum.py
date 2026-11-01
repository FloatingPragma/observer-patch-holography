"""All-energy fixed-wavepacket limit of the same K=0 contact action.

Exact rational certificates evaluate the proof in README.md. The independent
on-shell quadrature below is a floating diagnostic, never an interval proof.
"""
from __future__ import annotations
from fractions import Fraction as F
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import ellipkm1
from ..contact_continuum.continuum import rational
from ..contact_scattering.onshell import angular_density
from ..interacting_field.lattice import real


def slope_certificate():
    """a_1'(x)<=-1/16 on |x|<=1/4, retaining all eight nodes."""
    x=F(1,4); z=x*x/4
    correction=x*x*(3-z)/(18*(1-x*x)*(1-z)**2)
    exact_slope=F(1,12)-correction
    assert exact_slope>F(1,16)
    return {'scaled_energy_radius':x,'inverse_second_moment_lower':F(1,12),
            'positive_derivative_correction_upper':correction,
            'negative_slope_magnitude_lower':exact_slope,
            'rounded_slope_magnitude_lower':F(1,16),'all_nodes':8}


def window_certificate(energy_cap, *, speed, packet_density_cap=None, threshold=None):
    """Lebesgue integral of |S_c,U(E)-1|^2 on [-Emax,Emax].

    Uniform over all real finite U, including arbitrary U(c). A common
    energy-normalized bright-channel representation is explicit in README.
    This is not uniform operator-norm convergence over both energy signs.
    """
    e=rational(energy_cap,'energy cap',positive=True)
    c=rational(speed,'kinetic coefficient',positive=True)
    if e>c/4:raise ValueError('proved window requires Emax <= c/4')
    integral=min(8*e,F(1024,63)*e*e/c)
    result={'energy_cap':e,'speed':c,'integrated_squared_scattering_bound':integral,
            'all_real_contact_couplings':True}
    if packet_density_cap is not None:
        d=rational(packet_density_cap,'energy-normalized packet density cap')
        result['packet_squared_error_bound']=d*integral
    if threshold is not None:
        t=rational(threshold,'scattering threshold',positive=True)
        if t>2:raise ValueError('scattering threshold must be at most two')
        result['exceptional_energy_measure_upper']=min(2*e,integral/t**2)
    return result


def spacing_certificate(energy_cap, squared_error_tolerance, *, velocity,
                        packet_density_cap):
    e=rational(energy_cap,'energy cap',positive=True)
    tol=rational(squared_error_tolerance,'squared norm tolerance',positive=True)
    v=rational(velocity,'declared model velocity',positive=True)
    d=rational(packet_density_cap,'packet density cap',positive=True)
    spacing=min(v/(4*e),F(63,1024)*tol*v/(d*e*e))
    cert=window_certificate(e,speed=v/spacing,packet_density_cap=d)
    return {'spacing_upper':spacing,'velocity':v,'squared_error_tolerance':tol,
            'packet_squared_error_bound':cert['packet_squared_error_bound']}


def pv_correction_derivative_monomial(degree,x):
    """Diagnostic derivative of x^3 PV integral_0^1 t^(2n)/(x^2-t^2)."""
    if type(degree) is not int or not 0<=degree<=64:raise ValueError('degree 0..64 required')
    x=real(x,'positive scaled energy')
    if not 0<x<1:raise ValueError('0<x<1 required')
    positive=(2*degree+2)*x**(2*degree+1)*math.atanh(x)+x**(2*degree+2)/(1-x*x)
    negative=sum((2*k+3)*x**(2*k+2)/(2*(degree-1-k)+1) for k in range(degree))
    return positive-negative


def two_sine_density(s):
    """Diagnostic density of sin(k1)^2+sin(k2)^2 for uniform full torus.

    The arcsine convolution is 2 K(s(2-s))/pi^2. ellipkm1 evaluates its
    logarithmic singularity stably; independent tests integrate the original
    arcsine convolution. No such quadrature enters the exact certificate.
    """
    s=real(s,'sum of squared sines')
    if not 0<=s<=2:raise ValueError('s in [0,2] required')
    return float(2*ellipkm1((s-1)**2)/math.pi**2)


def principal_green_diagnostic(energy, *, speed=1.):
    """Actual on-shell real part, via analytic third-axis integration.

    For |E|/c<1/4, s0=(E/c)^2/4. Outside s>s0 the integrated third-axis
    resolvent is real, -x/[4 sqrt(s-s0)sqrt(s+1-s0)]. The substitution
    s=s0+u^2 removes its threshold endpoint. The interval s<s0 contributes
    only to the imaginary part. This numerical integral is not certified.
    """
    e=real(energy,'energy');c=real(speed,'speed')
    if c<=0 or abs(e)>c/4:raise ValueError('c>0 and |E|<=c/4 required')
    if e==0:return 0.
    x=e/c;s0=x*x/4;split=math.sqrt(1-s0)
    value,error=quad(lambda u:two_sine_density(s0+u*u)/math.sqrt(1+u*u),
                     0.,math.sqrt(2-s0),points=[split],epsabs=1e-10,epsrel=1e-10)
    if not math.isfinite(value) or error>1e-7:raise ArithmeticError('on-shell diagnostic unresolved')
    return -x*value/(2*c)


def signed_scattering(principal_green, density, interaction):
    a=real(principal_green,'real boundary value');rho=real(density,'density')
    u=real(interaction,'real contact coupling')
    if rho<=0:raise ValueError('strictly positive ac density required')
    if u==0:return 1.+0j
    # Divide by U to avoid overflow in large-coupling diagnostics.
    x=1/u-a;y=math.pi*rho
    if not math.isfinite(x):raise ValueError('inverse coupling outside floating domain')
    return complex((x-1j*y)/(x+1j*y))


def resonance_diagnostic(energy=-1/16, speeds=(1.,4.,16.,64.)):
    """Negative control: a tuned negative-energy resonance keeps S(E0)=-1.

    Its width collapses; the family still converges strongly on fixed packets.
    A positive-energy mirror uses attractive contact, outside the original
    repulsive source family but included in the quantified algebraic extension.
    """
    e=real(energy,'negative resonance energy')
    if e>=0:raise ValueError('negative target required for repulsive control')
    rows=[]
    for c in speeds:
        a=principal_green_diagnostic(e,speed=c)
        u=1/a
        rho=angular_density(e,speed=c)
        peak=signed_scattering(a,rho,u)
        offset=e/2
        off_a=principal_green_diagnostic(offset,speed=c)
        off=signed_scattering(off_a,angular_density(offset,speed=c),u)
        rows.append({'speed':c,'resonance_energy':e,'repulsive_U_diagnostic':u,
                     'U_over_c_squared_diagnostic':u/c**2,
                     'peak_scattering_distance':abs(peak-1),
                     'off_resonance_energy':offset,'off_resonance_distance':abs(off-1)})
    return rows
