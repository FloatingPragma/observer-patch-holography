"""Independent complex source contour diagnostics; these are not certificates."""
import cmath
import math

import numpy as np
from scipy.integrate import quad


def radial_integrand(radius, omega, momentum, spacing, energy):
    w, k = np.asarray(omega, float), np.asarray(momentum, complex)
    a, z, r = float(spacing), complex(energy), complex(radius)
    C, S = np.cos(k/2), np.sin(k/2)
    R = cmath.sqrt(np.sum(C*C))
    t = R-a*a*z/(2*math.sqrt(3))
    lam = R*R-t*t
    # The stable expression is essential at the actual fine lattice spacing.
    lam = (R/math.sqrt(3))*a*a*z-a**4*z*z/12
    u = r*w
    root = np.sqrt(1-u*u)
    F = np.sum(np.cos(k)*u*u)+(np.sum(C*S*u*root))**2/(t*t)
    N = 1-np.sum(S*S*u*u)/(t*t)
    J = 1/np.prod(root)
    return r*r*J*N/(F-lam)


def complex_integral(function, lower, upper, **kwargs):
    return (quad(lambda x: function(x).real, lower, upper, **kwargs)[0]+
            1j*quad(lambda x: function(x).imag, lower, upper, **kwargs)[0])


def deformed_inner_ray(omega, momentum, spacing, energy, *, orientation=1):
    phase = cmath.exp(orientation*1j*math.pi/4)
    ray = complex_integral(
        lambda x: phase*radial_integrand(phase*x, omega, momentum, spacing, energy),
        0, .5, epsabs=2e-12, epsrel=2e-12, points=[float(spacing)])
    arc = complex_integral(
        lambda theta: 1j*.5*cmath.exp(1j*theta)*radial_integrand(
            .5*cmath.exp(1j*theta), omega, momentum, spacing, energy),
        orientation*math.pi/4, 0, epsabs=2e-12, epsrel=2e-12)
    return ray+arc


def original_inner_ray(omega, momentum, spacing, energy):
    if complex(energy).imag >= 0:
        raise ValueError('positive-imaginary source energy requires Im(s)<0')
    a = float(spacing)
    return complex_integral(
        lambda r: radial_integrand(r, omega, momentum, spacing, energy),
        0, .5, epsabs=2e-12, epsrel=2e-12,
        points=[a/2, a, 2*a, 4*a], limit=150)
