"""Finite adjacent-rod survey on the transverse product, with integer readers."""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import ceil, isqrt, pi

from .geometry import Ball

# The closed seven-point Newton–Cotes rule on [0,1] is positive and exact
# through degree seven. Its nodes are material radial register addresses.
RADIAL = tuple((F(i, 6), F(w, 840)) for i, w in enumerate((41, 216, 27, 272, 27, 216, 41)))


def _atan_inverse(q, tolerance):
    power = F(1, q)
    acc = F(0)
    j = 0
    while abs(power)/(2*j+1) > tolerance:
        acc += power/(2*j+1)
        power *= F(-1, q*q)
        j += 1
    return acc


def angle_gain(ell):
    """Finite Machin sum for 2π, absolute error below ell^4/16."""
    tolerance = F(ell)**4/640
    return 32*_atan_inverse(5, tolerance) - 8*_atan_inverse(239, tolerance)


def _configuration(ell, u):
    ell = F(ell)
    if not 0 < ell <= F(1, 32):
        raise ValueError("the finite rod survey is enrolled at radii up to 1/32")
    # The registered exact frame components have denominators at most three.
    frame = tuple(F(x).limit_denominator(3) for x in u)
    if (len(frame) != 4 or any(abs(float(a)-float(b)) > 1e-14 for a,b in zip(frame, u, strict=True))
            or frame[0] <= 0 or frame[0]**2-sum(x*x for x in frame[1:]) != 1):
        raise ValueError("expected an enrolled unit timelike frame")
    b2 = 1+frame[1]**2+frame[2]**2
    if b2 > 9:
        raise ValueError("frame exceeds the enrolled transverse stretch")
    return ell, b2, ceil(1/ell)


def rod_ratio(z, n, precision):
    """Final separation / radial length. z=K s². Keep integer pointer ticks."""
    if not 0 <= z <= F(1, 2) or n <= 0 or precision < n**8:
        raise ValueError("rod recurrence outside its certified range")
    z = F(z)/(n*n)
    previous, current = 0, (2*precision+n)//(2*n)
    for _ in range(1, n):
        # Round z*current to the nearest integer without a float conversion.
        recoil = (2*z.numerator*current+z.denominator)//(2*z.denominator)
        previous, current = current, 2*current-previous-recoil
    return F(current, precision)


def _sqrt_down(x, ticks):
    return F(isqrt(x.numerator*ticks*ticks//x.denominator), ticks)


def _mean_density(K, r, b2, n):
    """Average six resolved rest-slice polar rays. Antipodes have equal counts."""
    if K == 0 or r == 0:
        return F(1)
    ticks = n**8
    # Principal directions of the projection have stretches b, 1, 0. Along
    # the third axis the transverse exponential is evaluated at its origin.
    along_b = rod_ratio(K*r*r*b2, n, ticks)
    along_one = rod_ratio(K*r*r, n, ticks)
    transverse = _sqrt_down(1+b2*(along_one*along_one-1), ticks)
    return (along_b+transverse+1)/3


def _volume_over_4pi(K, r, b2, n):
    return r**3*sum(w*t*t*_mean_density(K, r*t, b2, n) for t, w in RADIAL)


@dataclass(frozen=True)
class RodRead:
    radius: F
    stop_width: F
    area_over_4pi: F
    volume_over_4pi: F
    weighted_over_4pi: F


@lru_cache(maxsize=256)
def _stop(K, ell, b2, n):
    if not 0 <= K < 13:
        raise ValueError("the largest enrolled curvature is 4π")
    lo, hi = ell/2, 2*ell
    while hi-lo > ell**5:
        mid = (lo+hi)/2
        if _volume_over_4pi(K, mid, b2, n) < ell**3/3:
            lo = mid
        else:
            hi = mid
    r = (lo+hi)/2
    volume = _volume_over_4pi(K, r, b2, n)
    area = r*r*_mean_density(K, r, b2, n)
    weight = _weighted(K, r, b2, n)
    return RodRead(r, hi-lo, area, volume, weight)


def _weighted(K, r, b2, n):
    return r**4/2*sum(w*t*t*(1-t*t)*_mean_density(K, r*t, b2, n) for t, w in RADIAL)


def read_ball(K, ell, u):
    """All rational terminal readers with the common solid-angle factor removed."""
    ell, b2, n = _configuration(ell, u)
    return _stop(F(K), ell, b2, n)


def measured_ball(K, ell, u):
    """Floating presentation of the finite radius, volume and area registers."""
    value = read_ball(K, ell, u)
    return Ball(float(value.radius), 4*pi*float(value.volume_over_4pi),
                4*pi*float(value.area_over_4pi))


def weighted_volume(K, radius, u, *, ell):
    """Finite radial integration of the action's ball weight at a given rod stop."""
    ell, b2, n = _configuration(ell, u)
    r = F(radius)
    if not (0 < r <= 2*ell and 0 <= K < 13):
        raise ValueError("weighted rod reading outside the enrolled range")
    return 4*pi*float(_weighted(F(K), r, b2, n))
