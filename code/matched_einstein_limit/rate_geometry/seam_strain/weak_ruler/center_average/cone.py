"""Sharp coefficients for volume-matching before and after observer averaging."""

from math import inf, nextafter

from mpmath.ctx_iv import MPIntervalContext


_iv = MPIntervalContext()
_iv.dps = 50


def _bounds(x):
    return nextafter(float(x.a), -inf), nextafter(float(x.b), inf)


def _half_gap(alpha):
    alpha = _iv.mpf(alpha)
    if not (_iv.pi < alpha and alpha < 2 * _iv.pi):
        raise ValueError("the cone angle must lie strictly between π and 2π")
    return _iv.pi - alpha / 2


def _mean_volume(alpha):
    eta = _half_gap(alpha)
    return 2 * alpha / 15 * (eta / _iv.sin(eta)**2 - _iv.cot(eta))


def mean_volume_coefficient(alpha):
    """F such that the missing volume at a common radius is F r⁵ per unit hinge."""
    return float(_mean_volume(_iv.mpf(alpha)).mid)


def mean_volume_interval(alpha_lo, alpha_hi):
    """Outward bounds on F over the full cone-angle interval in radians."""
    return _bounds(_mean_volume(_iv.mpf([alpha_lo, alpha_hi])))


def _profile(eta, derivatives=False):
    """Volume and area fractions b(u), p(u), with u=distance-to-hinge/radius."""
    sn, cs = _iv.sin(eta), _iv.cos(eta)

    def fractions(u):
        if u == 0:
            values = (1 - eta / _iv.pi,) * 2
            return values + (1 / (2 * _iv.pi),) * 2 if derivatives else values
        if u == 1:
            values = (_iv.mpf(1),) * 2
            return values + (_iv.mpf(0),) * 2 if derivatives else values
        h = _iv.sqrt((1 - u) * (1 + u))
        y = u * sn
        lam = _iv.atan2(sn * h, cs)
        psi = _iv.atan2(h, u * cs)
        missing_area = lam - y * psi
        missing_volume = lam - y * (3 - y*y) * psi / 2 + u*y*cs*h / 2
        values = (1 - missing_volume / _iv.pi, 1 - missing_area / _iv.pi)
        if not derivatives:
            return values
        v = u * cs
        dp = (h - v * psi) / (2 * _iv.pi)
        db = (2*h**3 + 3*v**2*h - 3*v*(h**2 + v**2)*psi) / (4 * _iv.pi)
        return *values, db, dp

    return fractions


def individual_interval(alpha_lo, alpha_hi, bins=512):
    """Enclose the area deficit coefficient for individually matched ball volumes.

    The input interval is the full cone angle in radians. Inputs may be
    binary floats or exact decimal strings. Both returned floats round outward.
    """
    if _iv.mpf(alpha_lo) > _iv.mpf(alpha_hi) or not isinstance(bins, int) or bins < 1:
        raise ValueError("ordered angle endpoints and a positive bin count are required")
    profiles_lo = _profile(_half_gap(alpha_lo))
    profiles_hi = _profile(_half_gap(alpha_hi))
    two_thirds = _iv.mpf(2) / 3
    five_thirds = _iv.mpf(5) / 3
    lower = upper = _iv.mpf(0)
    b_lo, p_lo = profiles_lo(_iv.mpf(0))
    for i in range(bins):
        right = _iv.mpf(i + 1) / bins
        b_hi, p_hi = profiles_hi(right)
        weight = _iv.mpf(2*i + 1) / (2 * bins**2)
        loss_lo = max(0, (1 - p_hi * b_lo**(-two_thirds)).a)
        loss_hi = 1 - p_lo * b_hi**(-two_thirds)
        lower += weight * p_lo * b_hi**(-five_thirds) * loss_lo
        upper += weight * p_hi * b_lo**(-five_thirds) * loss_hi
        b_lo, p_lo = profiles_lo(right)
    return (_bounds(4 * _iv.pi * _iv.mpf(alpha_lo) * lower)[0],
            _bounds(4 * _iv.pi * _iv.mpf(alpha_hi) * upper)[1])


def individual_slope_interval(alpha_lo, alpha_hi, bins=4096):
    """Enclose dI/dα throughout the full cone-angle interval, in radians."""
    if _iv.mpf(alpha_lo) > _iv.mpf(alpha_hi) or not isinstance(bins, int) or bins < 1:
        raise ValueError("ordered angle endpoints and a positive bin count are required")
    alpha = _iv.mpf([alpha_lo, alpha_hi])
    eta_lo, eta_hi = _half_gap(alpha_hi), _half_gap(alpha_lo)
    low_profile, high_profile = _profile(eta_hi, True), _profile(eta_lo, True)
    total = _iv.mpf(0)
    b_lo, p_lo, db_hi, dp_hi = low_profile(_iv.mpf(0))
    for i in range(bins):
        right = _iv.mpf(i + 1) / bins
        b_hi, p_hi, db_lo, dp_lo = high_profile(right)
        p = _iv.mpf([p_lo.a, p_hi.b])
        b = _iv.mpf([b_lo.a, b_hi.b])
        # Both angle partials decrease as u or α increases.
        dp = _iv.mpf([dp_lo.a, dp_hi.b])
        db = _iv.mpf([db_lo.a, db_hi.b])
        gain = p**2 * b**(-_iv.mpf(7)/3)
        dgain = gain * (2*dp/p - 7*db/(3*b))
        weight = _iv.mpf(2*i + 1) / (2 * bins**2)
        total += weight * (gain + alpha*dgain)
        b_lo, p_lo, db_hi, dp_hi = low_profile(right)
    return _bounds(2 * _iv.pi - 4 * _iv.pi * total)


def individual_ball(alpha, offset, bins=64):
    """Enclose radius, area and volume when the target volume is 4π/3.

    The center's transverse distance is offset in units of the target flat
    radius. Here bins is the number of bisections of the normalized offset.
    """
    profile = _profile(_half_gap(alpha))
    x = _iv.mpf(offset)
    if x < 0 or not isinstance(bins, int) or bins < 1:
        raise ValueError("a nonnegative offset and a positive bisection count are required")
    lo, hi = _iv.mpf(0), _iv.mpf(1)
    third = _iv.mpf(1) / 3
    if x >= 1:
        lo = hi
    elif x == 0:
        hi = lo
    else:
        for _ in range(bins):
            mid = (lo + hi) / 2
            b, _ = profile(mid)
            error = mid * b**(-third) - x
            if error < 0:
                lo = mid
            elif error > 0:
                hi = mid
            else:
                # dx/du = p b^(-4/3) ≥ 1/2.
                lo = max(lo, (mid - 2*abs(error)).a)
                hi = min(hi, (mid + 2*abs(error)).b)
                break
    b_lo, p_lo = profile(lo)
    b_hi, p_hi = profile(hi)
    radius_lo, radius_hi = b_hi**(-third), b_lo**(-third)
    return {
        "radius": (_bounds(radius_lo)[0], _bounds(radius_hi)[1]),
        "area": (_bounds(4*_iv.pi * radius_lo**2 * p_lo)[0],
                 _bounds(4*_iv.pi * radius_hi**2 * p_hi)[1]),
        "volume": _bounds(4*_iv.pi / 3),
    }
