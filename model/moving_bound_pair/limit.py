"""Exact moving-pole/propagator budgets and separate numerical diagnostics.

The single global comparison is J_a=6/a^2-H_a. The momentum-dependent free
edge is used in the proof, never subtracted as a different Hamiltonian.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
import math
import numpy as np
from scipy.optimize import brentq
from ..band_edge_contact import limit as edge
from ..contact_continuum.continuum import rational
from ..interacting_field.lattice import real


def green_comparison_certificate():
    # kappa=max|K_j|<=1/2; h=a^2 z<=1/8, z>=0.
    coefficient = F(5, 8) + F(49, 96)*(F(20, 7)+F(32, 315))
    assert coefficient < F(9, 4) < 3
    return {'maximum_lattice_momentum': F(1, 2),
            'maximum_scaled_exterior_offset': F(1, 8),
            'denominator_ratio_lower': F(7, 8),
            'squared_scaled_edge_lower': F(45, 16),
            'green_difference_coefficient': coefficient,
            'rounded_green_difference_coefficient': 3,
            'zero_momentum_exterior_coefficient': 11,
            'all_edge_valleys': 8}


def moving_pole_certificate(spacing, gamma, momentum_cap):
    a = rational(spacing, 'spacing', positive=True)
    g = rational(gamma, 'gamma', positive=True)
    p = rational(momentum_cap, 'componentwise rescaled model momentum cap')
    if 256*a*(g*g+p*p) > g:
        raise ValueError('requires 256*a*(gamma^2+momentum_cap^2)<=gamma')
    # AM-GM gives a*p<=1/512, and a*gamma<=1/256.
    assert a*p <= F(1, 512) and a*g <= F(1, 256)
    pole_error = F(22, 7)*a*g*(132*g*g+9*p*p)
    edge_error = a*a*p**4/4
    return {'spacing': a, 'gamma': g, 'momentum_cap': p,
            'unique_upper_pole_for_every_momentum': True,
            'pole_offset_lower': g*g/4, 'pole_offset_upper': 4*g*g,
            'pole_offset_limit': g*g, 'pole_offset_error_upper': pole_error,
            'center_of_mass_edge_error_upper': edge_error,
            'reflected_band_energy_error_upper': pole_error+edge_error,
            'reflected_band_energy_limit': '|P|^2/4-gamma^2',
            'same_global_energy_shift': '6/a^2',
            'physical_rest_mass_inferred': False}


def propagation_certificate(spacing, gamma, momentum_cap, duration):
    t = rational(duration, 'absolute duration')
    cert = moving_pole_certificate(spacing, gamma, momentum_cap)
    return {**cert, 'absolute_duration': t,
            'operator_error_upper': min(F(2), t*cert['reflected_band_energy_error_upper']),
            'comparison_map': 'explicit normalized resolvent bound-fiber isometry',
            'arbitrary_untouched_references': True}


def spacing_certificate(gamma, momentum_cap, duration, tolerance):
    g = rational(gamma, 'gamma', positive=True)
    p = rational(momentum_cap, 'momentum cap')
    t = rational(duration, 'absolute duration')
    eps = rational(tolerance, 'tolerance', positive=True)
    if eps > 2:
        raise ValueError('operator tolerance must be at most two')
    coefficient = F(22, 7)*g*(132*g*g+9*p*p)+p**4/4
    a = min(F(1), g/(256*(g*g+p*p)))
    if t:
        a = min(a, eps/(t*coefficient))
    return {'spacing_upper': a, 'tolerance': eps,
            **propagation_certificate(a, g, p, t)}


def exterior_transition_certificate(spacing, gamma, momentum_cap, offset_cap, root_gap):
    a = rational(spacing, 'spacing', positive=True)
    g = rational(gamma, 'gamma', positive=True)
    p = rational(momentum_cap, 'momentum cap')
    z = rational(offset_cap, 'exterior offset cap', positive=True)
    gap = rational(root_gap, 'distance of sqrt(z) from gamma', positive=True)
    if a*p > F(1, 2) or a*a*z > F(1, 8) or a*g > 1:
        raise ValueError('outside the proved momentum/energy/coupling domain')
    error = a*(11*z+3*p*p)
    if 44*error > 7*gap:
        raise ValueError('perturbation is too large for the supplied pole separation')
    return {'green_error_upper': error,
            'renormalized_transition_error_upper': F(968, 49)*error/gap**2,
            'comparison': 'a^3 tau_a(B_K+z) -> pi/(sqrt(z)-gamma)',
            'domain': '0<=z<=offset_cap and |sqrt(z)-gamma|>=root_gap',
            'on_shell_limit_claimed': False}


def _momentum(value):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked momentum forbidden')
    p = np.asarray(value)
    if p.shape != (3,) or p.dtype.kind not in 'fiu' or not np.isfinite(p).all():
        raise ValueError('three finite real momentum components required')
    return p.astype(float)


def free_edge(momentum, speed):
    k = _momentum(momentum); c = real(speed, 'speed')
    if c <= 0 or np.max(abs(k)) > math.pi/2:
        raise ValueError('requires c>0 and |K_j|<=pi/2')
    return 2*c*math.sqrt(float(np.sum(np.cos(k/2)**2)))


def edge_deficit(spacing, momentum):
    a = real(spacing, 'spacing'); p = _momentum(momentum)
    if a <= 0 or np.max(abs(a*p)) > .5:
        raise ValueError('requires a>0 and max|aP_j|<=1/2')
    sine2 = float(np.sum(np.sin(a*p/2)**2))
    alpha = math.sqrt(1-sine2/3)
    return 2*sine2/(a*a*(1+alpha))


def singlet_kernel(relative_momentum, total_momentum, energy, speed):
    """Exact rational Pauli resolvent formula evaluated in floating arithmetic."""
    q = np.asarray(relative_momentum, dtype=float); k = _momentum(total_momentum)
    e = real(energy, 'energy'); c = real(speed, 'speed')
    if q.shape[-1:] != (3,) or not np.isfinite(q).all() or c <= 0:
        raise ValueError('finite (...,3) relative momenta and positive speed required')
    if e <= free_edge(k, c):
        raise ValueError('exterior energy required')
    d = np.sin(q)*np.cos(k/2); v = np.cos(q)*np.sin(k/2)
    d2 = np.sum(d*d, axis=-1); v2 = np.sum(v*v, axis=-1)
    dv = np.sum(d*v, axis=-1)
    return e*(e*e-4*c*c*v2)/(e**4-4*c*c*(d2+v2)*e*e+16*c**4*dv*dv)


@lru_cache(maxsize=8)
def _pyramid_nodes(order):
    if type(order) is not int or not 8 <= order <= 96:
        raise ValueError('integer quadrature order 8..96 required')
    nodes, weights = np.polynomial.legendre.leggauss(order)
    radial = (nodes+1)*math.pi/4
    rr, uu, vv = np.meshgrid(radial, nodes, nodes, indexing='ij')
    weight = ((weights*math.pi/4)[:, None, None]*weights[None, :, None]*
              weights[None, None, :]*rr*rr)
    return rr, uu, vv, weight


def green_correction_diagnostic(spacing, momentum, offset, *, order=32):
    """(g_K(B_K+z)-g_0(B_0+z))/a^3, integrated over the full torus.

    Six cube pyramids remove the integrable origin singularity. Exact
    subtraction avoids cancellation of two large near-threshold Green
    functions. This is a diagnostic, not a certified quadrature interval.
    """
    a = real(spacing, 'spacing'); p = _momentum(momentum); z = real(offset, 'offset')
    if a <= 0 or z < 0 or a*a*z > 1/8 or np.max(abs(a*p)) > .5:
        raise ValueError('requires a>0,z>=0,a^2 z<=1/8,max|aP_j|<=1/2')
    k = a*p; s = np.sin(k/2); c = np.cos(k/2)
    alpha = math.sqrt(1-float(s@s)/3)
    da = -float(s@s)/(3*(1+alpha))
    h = a*a*z; tk2 = 3*alpha*alpha+alpha*h+h*h/12
    rr, uu, vv, weight = _pyramid_nodes(order)
    value = 0.
    for axis in range(3):
        vectors = [rr*uu, rr*vv]
        vectors.insert(axis, rr)
        r = np.stack(vectors, axis=-1)
        sr = np.sin(r); cr = np.cos(r)
        q = np.sum(sr*sr, axis=-1)
        v2 = np.sum((sr*s)**2, axis=-1)
        dv2 = np.sum(sr*cr*s*c, axis=-1)**2
        shape = -2*v2+dv2/tk2
        d0 = q+h+h*h/12; dk = d0+shape+da*h
        numerator = da*(q-h*h/12)-(1+h/6)*shape-(alpha+h/6)*v2/tk2*d0
        value += float(np.sum(weight*numerator/(dk*d0)))
    return value/(a*math.pi**3)  # 2 equal opposite pyramids and 1/(2a)


def scaled_green_diagnostic(spacing, momentum, offset, *, order=32):
    return (edge.scaled_green_diagnostic(spacing, offset).real+
            green_correction_diagnostic(spacing, momentum, offset, order=order))


def pole_diagnostic(spacing, gamma, momentum, *, order=32):
    a = real(spacing, 'spacing'); g = real(gamma, 'gamma'); p = _momentum(momentum)
    cap = float(np.max(abs(p)))
    if a <= 0 or g <= 0 or 256*a*(g*g+cap*cap) > g*(1+1e-14):
        raise ValueError('outside certified moving-pole bracket')
    delta = brentq(lambda z: -g/math.pi-scaled_green_diagnostic(a,p,z,order=order),
                    g*g/4,4*g*g,xtol=3e-10)
    energy = edge_deficit(a,p)-delta
    target = float(p@p)/4-g*g
    return {'offset':delta, 'reflected_band_energy':energy,
            'target':target, 'absolute_energy_error':abs(energy-target),
            'quadrature_order':order}
