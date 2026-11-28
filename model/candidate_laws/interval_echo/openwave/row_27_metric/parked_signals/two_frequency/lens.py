"""Reconstruct the spherical clock/radar records and retrace each enrolled ray."""

from fractions import Fraction as F
from math import asin, atan2, cos, hypot, pi, sin, sqrt


def crossing(mu, angular, state, changed):
    def rhs(s):
        r, _, p = s
        v = -sqrt(2 * mu / r)
        c, dc = 1.0, 0.0
        if changed and r < 4.0:
            x = 1 - r / 4.0
            c = 1 + 2 * mu * x * x / r
            dc = -2 * mu * (x * x / r**2 + x / (2 * r))
        tangent = c * angular / r
        w = hypot(p, tangent)
        return (
            v + p / w,
            c * c * angular / (w * r * r),
            v * p / (2 * r) + tangent * tangent * (1 / r - dc / c) / w,
        )

    def step(s, h):
        a = rhs(s)
        b = rhs(tuple(x + h * y / 2 for x, y in zip(s, a)))
        c = rhs(tuple(x + h * y / 2 for x, y in zip(s, b)))
        d = rhs(tuple(x + h * y for x, y in zip(s, c)))
        return tuple(
            x + h * (u + 2 * v + 2 * w + y) / 6 for x, u, v, w, y in zip(s, a, b, c, d)
        )

    h = 1 / 128
    for i in range(32 * 128):
        new = step(state, h)
        if new[0] < 1:
            raise ArithmeticError("Enrolled ray left exterior")
        if state[0] < 4 <= new[0]:
            lo, hi = 0.0, h
            for _ in range(42):
                mid = (lo + hi) / 2
                if step(state, mid)[0] < 4:
                    lo = mid
                else:
                    hi = mid
            return i * h + hi, step(state, hi)
        state = new
    raise ArithmeticError("Missing mirror or missing return")


def pulse(mu, sine, changed):
    v = -sqrt(mu / 2)
    lapse = sqrt(1 - mu / 2)
    L = 4 * sine
    p = (-sqrt(1 - sine**2) - v) / lapse
    t, (_, phi, p) = crossing(mu, L, (4.0, 0.0, p), changed)
    w = hypot(p, sine)
    energy, outward = (w + v * p) / lapse, (p + v * w) / lapse
    tb, (_, home, _) = crossing(
        mu, -L, (4.0, phi, (-outward - v * energy) / lapse), changed
    )
    assert abs(home) < 1e-9
    return phi + asin(sine) + asin(sine / energy) - pi, (t + tb) * lapse, phi


def run(records):
    words = list(
        zip(
            records["published_source"]["moving_public"],
            records["published_source"]["home_public"],
        )
    )
    n = len(words)
    ordinary = F(sum(a for a, _ in words), n)
    together = F(sum(a * b for a, b in words), n)
    assert (n, ordinary, together) == (33, F(1, 3), F(8, 33))
    tick = records["angular_and_time_tick"]
    calib = records["ordinary_release_calibration"]
    measured = calib["ordinary_radius_ticks"] * tick

    def gain(radius):
        theta = atan2(sqrt(2 - radius), sqrt(radius))
        return 16 * (theta + sin(2 * theta) / 2) ** 2 / float(ordinary)

    low = gain(measured + calib["radius_error_upper"])
    high = gain(measured - calib["radius_error_upper"])
    assert abs(low - calib["gain"][0]) < 5e-13 and abs(high - calib["gain"][1]) < 5e-13
    inferred = gain(measured)
    result = {}
    for name, stock in (
        ("blank", F(0)),
        ("ordinary", ordinary),
        ("held", ordinary + together),
    ):
        arm = records["sources"][name]
        assert F(arm["stock"]) == stock
        mu = inferred * float(stock)
        a, b = arm["emitted_inner_ticks"]
        c, d = arm["received_outer_ticks"]
        freq = (b - a) / (d - c)
        assert abs(freq - sqrt((1 - mu) / (1 - mu / 2))) < 1e-10
        rows = []
        for s in arm["shots"]:
            values = {}
            for key, changed in (("reading", False), ("angular_light_control", True)):
                if changed and name != "held":
                    continue
                raw = s[key]
                bend, time, address = pulse(mu, s["local_sine"], changed)
                blank_bend, blank_time, blank_addr = pulse(0.0, cos(address / 2), False)
                assert abs(address - blank_addr) < 1e-9 and abs(blank_bend) < 1e-9
                delay = time - blank_time
                measured_angle = (
                    raw["mirror_address_ticks"]
                    + raw["emission_angle_ticks"]
                    + raw["reception_angle_ticks"]
                ) * tick - pi
                measured_delay = (
                    raw["echo_clock_ticks"] - raw["matched_blank_echo_clock_ticks"]
                ) * tick
                assert (
                    abs(measured_angle - bend) < 1e-9
                    and abs(measured_delay - delay) < 1e-9
                )
                for target, read in (("bending", bend), ("delay", delay)):
                    lo, hi = s["forecast"][target]
                    assert (lo <= read <= hi) != changed
                values[key] = [round(bend, 8), round(delay, 8)]
            rows.append(values)
        result[name] = {"frequency_ratio": round(freq, 10), "shots": rows}
    return result
