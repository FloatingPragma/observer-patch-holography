"""Exact positive-energy scattering bounds for the unchanged Weyl/contact law.

The analytic proof is in README.md.  These rational routines evaluate its
constants; finite-broadening quadrature is an independent diagnostic only.
The limit has fixed model energy zero, K=0 and c=v/a.  It is not a theorem
about a filled vacuum, arbitrary counterterms, or the whole OPH class.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

import numpy as np

from ..interacting_field.lattice import real


def rational(value, name, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError(f'exact rational {name} required')
    value = F(value)
    if value < 0 or (positive and value == 0):
        raise ValueError(f'{name} must be {"positive" if positive else "nonnegative"}')
    return value


def principal_value_certificate():
    """A source-derived sign bound, with no numerical principal-value input.

    mu is the unit-mass contact spectral measure at c=1.  Its support is
    inside [-sqrt(12),sqrt(12)], so integral t^-2 dmu >= 1/12.
    The exact eight-node local density and pi>3 bound the only positive
    correction to a(x)/x by 8/405 for 0<x<=1/2.
    """
    ratio = F(1, 2)
    inverse_moment_lower = F(1, 12)
    correction_upper = 2*ratio**2/(36*(1-ratio**2/4)*(1-ratio**2))
    gap = inverse_moment_lower-correction_upper
    assert correction_upper == F(8, 405) and gap > F(1, 16)
    return {'ratio_upper': ratio, 'inverse_moment_lower': inverse_moment_lower,
            'positive_correction_upper': correction_upper,
            'negative_slope_magnitude_lower': gap,
            'rounded_negative_slope_magnitude_lower': F(1, 16),
            'all_nodes': 8, 'contact_measure_mass': F(1)}


def positive_energy_domain(energy, speed):
    energy = rational(energy, 'positive model energy', positive=True)
    speed = rational(speed, 'lattice kinetic coefficient', positive=True)
    if energy > speed/2:
        raise ValueError('proved domain requires 0 < E <= c/2')
    return energy, speed


def green_upper(energy, *, speed=F(1, 8)):
    """Certified upper bound for Re g_c(E+i0), strictly negative."""
    energy, speed = positive_energy_domain(energy, speed)
    return -energy/(16*speed**2)


def scattering_bound(energy, *, speed=F(1, 8), interaction=None):
    """Upper bound on ||S(E)-I||, uniformly in all repulsive U if omitted.

    The rank-one bright channel is the only nonidentity ac channel. Point
    spectrum is not part of S. No a-dependent choice of finite U>=0 escapes
    the universal bound. An optional specified U gives a second bound.
    """
    energy, speed = positive_energy_domain(energy, speed)
    universal = min(F(2), F(128, 45)*energy/speed)
    bounded_u = None
    if interaction is not None:
        interaction = rational(interaction, 'repulsive contact')
        bounded_u = F(8, 45)*interaction*energy**2/speed**3
    return {'universal_repulsive_bound': universal,
            'specified_repulsion_bound': bounded_u,
            'combined_bound': universal if bounded_u is None else min(universal, bounded_u),
            'principal_green_upper': green_upper(energy, speed=speed),
            'density_upper': F(4, 135)*energy**2/speed**3}


def spacing_certificate(energy_cap, tolerance, *, velocity=F(1, 8)):
    """A uniform positive-energy window certificate, c=v/a.

    Every 0<E<=Emax and 0<a<=spacing_upper obeys ||S-I||<=tolerance,
    independently of any nonnegative repulsion U(a). v has model units.
    """
    energy_cap = rational(energy_cap, 'energy cap', positive=True)
    tolerance = rational(tolerance, 'tolerance', positive=True)
    velocity = rational(velocity, 'declared long-wavelength velocity', positive=True)
    if tolerance > 2:
        raise ValueError('nontrivial scattering tolerance must be at most two')
    spacing = min(velocity/(2*energy_cap), F(45, 128)*tolerance*velocity/energy_cap)
    return {'spacing_upper': spacing, 'energy_cap': energy_cap,
            'velocity': velocity, 'all_repulsive_U': True,
            'error_upper': F(128, 45)*energy_cap*spacing/velocity}


def kinetic_error(momentum, spacing, *, velocity=F(1, 8)):
    """Each of all eight Weyl nodes: ||h(pi*b+a*p)+v sum (-1)^b sigma*p||.

    The common upper bound uses |sin x-x|<=|x|^3/6 and the triangle
    inequality. It is not a decoupling map that removes any node.
    """
    if type(momentum) not in (list, tuple) or len(momentum) != 3:
        raise ValueError('three exact signed momentum components required')
    if any(type(x) not in (int, F) for x in momentum):
        raise ValueError('exact signed momentum components required')
    spacing = rational(spacing, 'spacing', positive=True)
    velocity = rational(velocity, 'velocity', positive=True)
    return velocity*spacing**2*sum(abs(F(x))**3 for x in momentum)/6


def pv_monomial(degree, energy):
    """Floating diagnostic for PV integral_0^1 t^(2n)/(x^2-t^2) dt."""
    if type(degree) is not int or not 0 <= degree <= 64:
        raise ValueError('integer degree in 0..64 required')
    x = real(energy, 'scaled energy')
    if not 0 < x < 1:
        raise ValueError('scaled energy must lie strictly inside (0,1)')
    base = math.atanh(x)/x
    return x**(2*degree)*base-sum(x**(2*k)/(2*(degree-1-k)+1)
                                        for k in range(degree))


def broadened_green(energy, broadening, *, speed=1., order=96):
    """Independent two-dimensional quadrature after one exact sine integral.

    The returned value has eta>0; it neither evaluates nor certifies the
    boundary value used in the theorem. The complete Brillouin zone stays.
    """
    e, eta, c = real(energy, 'energy'), real(broadening, 'broadening'), real(speed, 'speed')
    if eta <= 0 or c <= 0 or type(order) is not int or not 8 <= order <= 512:
        raise ValueError('positive eta,c and integer quadrature order 8..512 required')
    z = complex(e, eta)
    k = 2*np.pi*(np.arange(order)+.5)/order
    sine = np.sin(k)**2
    a = z*z-4*c*c*(sine[:, None]+sine[None, :])
    # At E=0, +i0 is the branch continuous from E>0. Both negative-real
    # square roots are +i times their magnitude, making their product real
    # negative, as required by the original resolvent integral.
    denominator = np.sqrt(a+0j)*np.sqrt(a-4*c*c+0j)
    return complex(np.mean(z/denominator))
