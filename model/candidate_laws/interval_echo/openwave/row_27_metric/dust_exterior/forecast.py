"""Enclose predictions on the static exterior without using the ray tracer."""

from math import inf, nextafter
from mpmath.ctx_iv import MPIntervalContext


def outward(x):
    return [nextafter(float(x.a), -inf), nextafter(float(x.b), inf)]


def calibration(ordinary, stock):
    """Invert an ordinary released-clock appointment including its meter error."""
    iv = MPIntervalContext()
    iv.dps = 45
    sample = ordinary["released_clocks"]["1"]["1/2"]
    tick = iv.mpf(2)**-44
    radius = sample["radius_ticks"] * tick
    err = sample["radius_error_upper_by_refinement"][-1] + tick / 2
    measured = radius + iv.mpf([-1, 1]) * err
    theta = iv.atan2(iv.sqrt(2 - measured), iv.sqrt(measured))
    mu = 16 * (theta + iv.sin(2 * theta)/2)**2
    gain = mu * stock.denominator / stock.numerator
    return {"ordinary_radius_ticks": sample["radius_ticks"], "radius_error_upper": outward(err)[1],
            "gain": outward(gain)}


def calibrated_mu(gain, mass):
    if not mass:
        return [0., 0.]
    iv = MPIntervalContext()
    iv.dps = 42
    return outward(iv.mpf(gain)*mass.numerator/mass.denominator)


def forecast(mu_bounds, sine, radius=4, panels=512):
    iv = MPIntervalContext()
    iv.dps = 42
    mu = iv.mpf(mu_bounds)
    s = iv.mpf(sine)
    f = 1 - 2*mu/radius
    b2 = (radius*s)**2/f
    # r/sqrt(f(r)) is increasing on the entire registered exterior r >= 1 > 3mu.
    def cubic(r):
        return r**3 - b2*(r-2*mu)
    lo, hi = 1., float(radius)
    if 3*mu.b >= lo or cubic(iv.mpf(lo)).b >= 0 or cubic(iv.mpf(hi)).a <= 0:
        raise ValueError("No certified exterior turning point")
    for _ in range(58):
        mid = (lo+hi)/2
        value = cubic(iv.mpf(mid))
        if value.b < 0:
            lo = mid
        elif value.a > 0:
            hi = mid
        else:
            break
    q = iv.mpf([lo, hi])
    k = 2*mu/q
    zmax = iv.sqrt(1-q/radius)
    width = zmax/panels

    def values(z, derivatives):
        z2 = z*z
        a = (2-3*k) + (3*k-1)*z2 - k*z2**2
        v = 1-z2
        d = 1-k+k*z2
        angle = 1/iv.sqrt(a)
        time = angle/(v**2*d)
        if not derivatives:
            return angle, time
        da = 2*(3*k-1)*z - 4*k*z**3
        dda = 2*(3*k-1) - 12*k*z2
        dl = -da/(2*a)
        ddl = -(dda/a-(da/a)**2)/2
        dt = dl + 4*z/v - 2*k*z/d
        ddt = ddl + 4*(1+z2)/v**2 - 2*k*(1-k-k*z2)/d**2
        return abs(angle*(dl**2+ddl)), abs(time*(dt**2+ddt))

    # The midpoint rule error per integrated length is <= width^2 sup|I''| / 24.
    sums = [iv.mpf(0), iv.mpf(0)]
    errors = [iv.mpf(0), iv.mpf(0)]
    groups = 32
    if panels % groups:
        raise ValueError("Panel count must be divisible by derivative groups")
    for i in range(groups):
        interval = zmax*iv.mpf([i, i+1])/groups
        for j, bound in enumerate(values(interval, True)):
            errors[j] += (zmax/groups)*width**2*bound/24
    for i in range(panels):
        for j, term in enumerate(values(width*(2*i+1)/2, False)):
            sums[j] += term
    angle, time = [value*width + iv.mpf([-1, 1])*error for value, error in zip(sums, errors)]
    sweep = 4*angle
    proper = 8*q*iv.sqrt((1-k)*f)*time
    local = iv.atan2(s, iv.sqrt(1-s**2))
    return {"impact_parameter": outward(iv.sqrt(b2)), "turning_radius": outward(q),
            "bending": outward(sweep+2*local-iv.pi),
            "round_trip": outward(proper), "delay": outward(proper-4*radius*iv.sin(sweep/2)),
            "clock_ratio_2_to_4": outward(iv.sqrt((1-mu)/f))}
