"""Independent source-matrix and subtracted principal-value diagnostics.

Quadrature values are diagnostics, never the proof's error bounds.
"""
import math
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq

from ..band_edge_contact.limit import W0_diagnostic


def kernel(relative_momentum, total_momentum, energy, speed):
    q, k = np.asarray(relative_momentum, float), np.asarray(total_momentum, float)
    d, e = np.cos(k/2)*np.sin(q), np.sin(k/2)*np.cos(q)
    d2, e2, de = np.sum(d*d, axis=-1), np.sum(e*e, axis=-1), np.sum(d*e, axis=-1)
    return energy*(energy**2-4*speed**2*e2)/(energy**4-
        4*speed**2*(d2+e2)*energy**2+16*speed**4*de**2)


def angular_amplitude(x, omega, momentum, t):
    """A_direction and its y^2 derivative, without finite differences."""
    w, k = np.asarray(omega, float), np.asarray(momentum, float)
    q = 1-x*w*w
    b = np.cos(k/2)*np.sin(k/2)
    value = np.sum(b*w*np.sqrt(q))
    first = -.5*np.sum(b*w**3/np.sqrt(q))
    second = -.25*np.sum(b*w**5/q**1.5)
    f = np.sum(np.cos(k)*w*w)+value*value/t**2
    fp = 2*value*first/t**2
    fpp = 2*(first*first+value*second)/t**2
    D, Dp = f+x*fp, 2*fp+x*fpp
    J = float(np.prod(q**(-.5)))
    Jlog = .5*np.sum(w*w/q)
    e = np.sum(np.sin(k/2)**2*w*w)/t**2
    N, Np = 1-x*e, -e
    A = J*N/(math.sqrt(f)*D)
    derivative = A*(Jlog+Np/N-fp/(2*f)-Dp/D)/D
    return {'A': A, 'dA_dy2': derivative, 'f': f, 'df_dx': fp,
            'd2f_dx2': fpp, 'D': D, 'dD_dx': Dp}


@lru_cache(maxsize=8)
def quadrature(order):
    if type(order) is not int or not 12 <= order <= 96:
        raise ValueError('quadrature order 12..96 required')
    return np.polynomial.legendre.leggauss(order)


def radial(l, vector, k, t):
    r = l*np.asarray(vector)
    sn, cs = np.sin(r), np.cos(r)
    sj, cj = np.sin(k/2), np.cos(k/2)
    e2 = np.sum((sj*sn)**2)
    de = np.sum(sj*cj*sn*cs)
    shape = np.sum(np.cos(k)*sn*sn)+de*de/(t*t)
    derivative = (np.sum(2*np.cos(k)*sn*cs*vector)+
                  2*de*np.sum(sj*cj*(cs*cs-sn*sn)*vector)/(t*t))
    return shape, 1-e2/(t*t), derivative


def radial_pv(vector, momentum, t, lam, *, order=32):
    """One entire cube pyramid ray with its one interior pole subtracted."""
    k, vector = np.asarray(momentum), np.asarray(vector)
    root = brentq(lambda l: radial(l, vector, k, t)[0]-lam, 0., .25, xtol=1e-15)
    _, n0, derivative = radial(root, vector, k, t)
    b0 = 2*root**3*n0/derivative
    nodes, weights = quadrature(order)
    locations, weights = (nodes+1)*math.pi/4, weights*math.pi/4
    integrand = []
    for l in locations:
        f, n, _ = radial(l, vector, k, t)
        divisor = l*l-root*root
        b = l*l*n*divisor/(f-lam)
        integrand.append((b-b0)/divisor)
    real = float(weights@integrand)-b0*math.atanh(2*root/math.pi)/root
    imag = -math.pi*b0/(2*root)
    return complex(real, imag)


def scaled_green(spacing, momentum, energy, *, order=24):
    """Actual E=B_K-s+i0, via radial PV subtraction on six full pyramids.

    This never replaces the boundary value by a small imaginary broadening.
    """
    a, s = float(spacing), float(energy)
    p = np.asarray(momentum, float)
    if (not np.isfinite(a+s) or p.shape != (3,) or not np.isfinite(p).all() or
            a <= 0 or s <= 0 or a*a*s > 1/128 or np.max(abs(a*p)) > .25):
        raise ValueError('positive a,s and the proved momentum/energy domain required')
    k = a*p
    R = math.sqrt(float(np.sum(np.cos(k/2)**2)))
    t = R-a*a*s/(2*math.sqrt(3))
    lam = a*a*s*(R/math.sqrt(3)-a*a*s/12)
    nodes, weights = quadrature(order)
    total = 0j
    for axis in range(3):
        for u, wu in zip(nodes, weights, strict=True):
            for v, wv in zip(nodes, weights, strict=True):
                vector = [u, v]; vector.insert(axis, 1.)
                total += wu*wv*radial_pv(vector, k, t, lam, order=order)
    W = (2/math.pi**3)*total
    return (t/math.sqrt(3)*W-W0_diagnostic())/(2*a)


def scattering(spacing, gamma, momentum, energy, *, order=24):
    value = scaled_green(spacing, momentum, energy, order=order)
    denominator = -float(gamma)/math.pi-value
    reflected = denominator/denominator.conjugate()
    target = (gamma-1j*math.sqrt(energy))/(gamma+1j*math.sqrt(energy))
    return {'reflected_phase': reflected, 'target': target,
            'absolute_error': abs(reflected-target), 'PV_Green': value}
