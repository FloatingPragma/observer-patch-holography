"""Exact volume-dependent full-Fock stability and complete-history bounds."""
from fractions import Fraction as F
from math import prod


def _parameters(sides, speed, interaction):
    if type(sides) is not tuple or len(sides) != 3 or any(
            type(length) is not int or length < 1 for length in sides):
        raise ValueError('three positive integer side lengths required')
    if any(type(value) not in (int, F) or value < 0 for value in (speed, interaction)):
        raise ValueError('nonnegative exact rational coefficients required')
    volume = prod(sides)
    edges = sum((sides[axis]-1)*prod(sides[j] for j in range(3) if j != axis)
                for axis in range(3))
    return volume, edges, F(speed), F(interaction)


def stability_bounds(sides, speed=F(1, 8), interaction=F(1)):
    volume, edges, c, u = _parameters(sides, speed, interaction)
    return {'sites': volume, 'spatial_edges': edges,
            'energy_lower': -c*edges, 'energy_upper': c*edges+u*volume,
            'generator_norm_sum_upper': c*edges+u*volume}


def complete_history_bound(sides, segments, speed=F(1, 8), interaction=F(1)):
    """A segment is (signed model duration, positive integer split count)."""
    b = stability_bounds(sides, speed, interaction)['generator_norm_sum_upper']
    total = F()
    for duration, steps in segments:
        if type(duration) not in (int, F) or type(steps) is not int or steps < 1:
            raise ValueError('exact signed duration and positive integer split count required')
        total += b**3*abs(F(duration))**3/(3*steps**2)
    return min(F(1), total)


def exponential_upper(value, order=16):
    """Exact rational enclosure from a geometric Taylor tail and scaling."""
    if type(value) not in (int, F) or value < 0 or type(order) is not int or not 1 <= order <= 64:
        raise ValueError('nonnegative rational exponent and Taylor order in 1..64 required')
    reduced, squarings = F(value), 0
    # Keep this public arithmetic helper bounded independently of field size.
    if reduced > 64:
        raise ValueError('exponent exceeds the bounded exact-enclosure domain')
    while reduced > F(1, 2):
        reduced /= 2; squarings += 1
    term = total = F(1)
    for n in range(1, order+1):
        term *= reduced/n; total += term
    next_term = term*reduced/(order+1)
    upper = total+next_term/(1-reduced/(order+2))
    return upper**(1 << squarings)


def local_boundary_bound(radius, support_radius, support_sites, duration, speed=F(1, 8)):
    """Exact upper bound for ||tau_larger(A)-tau_radius(A)||, ||A||<=1.

    A is supported on support_sites sites within [-support_radius,support_radius]^3.
    The comparison is at fixed spacing on full CAR algebras and is uniform in
    larger-box size, occupation number, repulsive onsite U and reference size.
    It is a local-observable bound, not global-state or global-channel distance.
    """
    if (type(radius) is not int or type(support_radius) is not int
            or not 0 <= support_radius < radius <= 4096 or type(support_sites) is not int
            or not 1 <= support_sites <= (2*support_radius+1)**3
            or type(duration) not in (int, F) or type(speed) not in (int, F) or speed < 0):
        raise ValueError('nested integer cubes, support count and exact time/speed required')
    t, c = abs(F(duration)), F(speed)
    if not c or not t:
        return F(0)
    exponent = 48*c*t
    # The trivial automorphism bound avoids unbounded rational exponentiation.
    if exponent > 64:
        return F(2)
    boundary_edges = 6*(2*radius+1)**2
    bound = (c*support_sites*boundary_edges*t*F(1, 2**(radius-support_radius))
             *(exponential_upper(exponent)-1))
    return min(F(2), bound)
