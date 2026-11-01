"""Exact error budgets; separate, uncertified numerical Green diagnostics.

The comparison J_a=6/a**2-H_a reverses model time. Its two-body scattering
limit is nontrivial. This is not a many-body or relativistic QFT limit.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
import math
from scipy.integrate import quad
from scipy.optimize import brentq
from ..contact_continuum.continuum import rational
from ..contact_wavepacket.continuum import two_sine_density
from ..interacting_field.lattice import real


def green_certificate():
    """Evaluate rational majorants in README, using 3 < pi < 22/7."""
    outside = 16 + F(10, 9)
    inside = 32 + F(16, 9) + F(16, 27) + F(8, 21)
    scaled_inside = 18 + F(5, 12) + F(1, 24)
    # 9*(1+h/6)*(1+h/12), h<=1/8, bounds the exterior remainder.
    scaled_outside = 9 * F(49, 48) * F(97, 96) + F(5, 12) + F(1, 24)
    assert outside < 18 and inside < 36
    assert scaled_outside < 11 and scaled_inside < 19
    return {'ball_radius': F(1, 2), 'shell_radius_squared_max': F(1, 8),
            'angular_coefficients_upper': 1, 'W0_lower': F(2, 3), 'W0_upper': 5,
            'outside_W_remainder_coefficient': outside,
            'inside_W_remainder_coefficient': inside,
            'outside_scaled_green_coefficient': scaled_outside,
            'inside_scaled_green_coefficient': scaled_inside,
            'uniform_scaled_green_coefficient': 19, 'edge_valleys': 8}


def scattering_certificate(spacing, gamma, energy_cap):
    """Uniform scalar contact-channel scattering error on 0<=s<=energy_cap."""
    a = rational(spacing, 'spacing', positive=True)
    g = rational(gamma, 'renormalized inverse length', positive=True)
    cap = rational(energy_cap, 'energy cap', positive=True)
    if a*g > 1 or a*a*cap > F(1, 8):
        raise ValueError('requires a*gamma<=1 and a^2*energy_cap<=1/8')
    if 836*a*cap > 7*g:
        raise ValueError('requires 836*a*energy_cap<=7*gamma')
    return {'spacing': a, 'gamma': g, 'energy_cap': cap,
            'green_error_upper': 19*a*cap,
            'denominator_modulus_lower_after_a3': 7*g/44,
            'uniform_scattering_error_upper': min(F(2), F(1672, 7)*a*cap/g),
            'coupling_positive': True, 'all_edge_valleys': 8,
            'source_phase': '(gamma+i*sqrt(s))/(gamma-i*sqrt(s))',
            'reflected_phase': '(gamma-i*sqrt(s))/(gamma+i*sqrt(s))'}


def spacing_certificate(gamma, energy_cap, tolerance):
    g = rational(gamma, 'gamma', positive=True)
    cap = rational(energy_cap, 'energy cap', positive=True)
    eps = rational(tolerance, 'tolerance', positive=True)
    if eps > 2:
        raise ValueError('scattering tolerance must be at most two')
    a = min(1/g, 1/(8*(1+cap)), 7*eps*g/(1672*cap))
    return {'spacing_upper': a, 'tolerance': eps,
            **scattering_certificate(a, g, cap)}


def pole_certificate(spacing, gamma):
    a = rational(spacing, 'spacing', positive=True)
    g = rational(gamma, 'gamma', positive=True)
    if a*g > F(1, 256):
        raise ValueError('pole bracket requires a*gamma<=1/256')
    # Endpoint signs: leading magnitude is at least 7*gamma/22 at upper,
    # and 7*gamma/44 at lower. Error bounds are 76*a*gamma^2 and 19/4*a*gamma^2.
    assert F(76, 256) < F(7, 22)
    assert F(19, 1024) < F(7, 44)
    return {'spacing': a, 'gamma': g, 'unique_above_band_pole': True,
            'pole_offset_lower': g*g/4, 'pole_offset_upper': 4*g*g,
            'pole_offset_limit': g*g,
            'pole_offset_error_upper': F(5016, 7)*a*g**3,
            'reflected_bound_energy_limit': -g*g,
            'original_energy': '6/a^2 + pole_offset',
            'physical_rest_mass_inferred': False}


def dispersion_certificate(spacing, momentum):
    a = rational(spacing, 'spacing', positive=True)
    momentum = tuple(momentum)
    if any(type(x) not in (int, F) for x in momentum):
        raise ValueError('exact rational momentum components required')
    p = tuple(F(x) for x in momentum)
    if len(p) != 3:
        raise ValueError('three momentum components required')
    norm2 = sum(x*x for x in p)
    return {'limit_energy': norm2,
            'absolute_error_upper': a*a*(sum(x**4 for x in p)+norm2**2)/3,
            'all_edge_valleys': 8}


def _quad(fun, lo, hi, points=()):
    value, error = quad(fun, lo, hi, points=list(points), epsabs=2e-11,
                        epsrel=2e-11, limit=300)
    if not math.isfinite(value) or error > 1e-8:
        raise ArithmeticError('Green diagnostic quadrature unresolved')
    return value


@lru_cache(maxsize=1)
def W0_diagnostic():
    """Uncertified full-torus integral, not an imported empirical constant."""
    return _quad(lambda u: 2*two_sine_density(u*u)/math.sqrt(1+u*u),
                 0, math.sqrt(2), [1])


def W_diagnostic(eta):
    """Boundary value W(eta+i0); negative eta is limited to [-1/8,0)."""
    x = real(eta, 'resolvent offset')
    if x == 0:
        return complex(W0_diagnostic())
    if x > 0:
        value = _quad(lambda u: 2*u*two_sine_density(u*u)/
                      math.sqrt((u*u+x)*(u*u+1+x)), 0, math.sqrt(2), [1])
        return complex(value)
    if x < -1/8:
        raise ValueError('interior diagnostic requires eta>=-1/8')
    t = -x
    re = _quad(lambda u: 2*two_sine_density(t+u*u)/math.sqrt(1+u*u),
               0, math.sqrt(2-t), [math.sqrt(1-t)])
    im = -2*math.sqrt(t)*_quad(
        lambda theta: math.sin(theta)*two_sine_density(t*math.sin(theta)**2)/
        math.sqrt(1-t*math.cos(theta)**2), 0, math.pi/2)
    return complex(re, im)


def scaled_green_diagnostic(spacing, offset):
    """(g_a(B+offset+i0)-g_a(B))/a^3; independent of the coupling."""
    a = real(spacing, 'spacing'); z = real(offset, 'edge offset')
    if a <= 0 or a*a*abs(z) > 1/8:
        raise ValueError('requires a>0 and a^2*abs(offset)<=1/8')
    h = a*a*z
    return ((1+h/6)*W_diagnostic(h+h*h/12)-W0_diagnostic())/(2*a)


def scattering_diagnostic(spacing, gamma, energy):
    a = real(spacing, 'spacing'); g = real(gamma, 'gamma'); s = real(energy, 'energy')
    if a <= 0 or g <= 0 or s < 0 or a*g > 1:
        raise ValueError('requires a,gamma>0, s>=0 and a*gamma<=1')
    d = -g/math.pi - scaled_green_diagnostic(a, -s)
    source = d.conjugate()/d
    target = (g+1j*math.sqrt(s))/(g-1j*math.sqrt(s))
    return {'source_phase': source, 'reflected_phase': source.conjugate(),
            'source_target': target, 'reflected_target': target.conjugate(),
            'absolute_error': abs(source-target)}


def pole_diagnostic(spacing, gamma):
    a = real(spacing, 'spacing'); g = real(gamma, 'gamma')
    if a <= 0 or g <= 0 or a*g > 1/256:
        raise ValueError('requires 0<a*gamma<=1/256')
    root = brentq(lambda z: -g/math.pi-scaled_green_diagnostic(a, z).real,
                  g*g/4, 4*g*g, xtol=1e-10)
    return {'offset': root, 'target': g*g, 'absolute_error': abs(root-g*g)}
