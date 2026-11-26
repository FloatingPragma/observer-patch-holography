"""Independent canonical light, mirror and freely released clock readouts."""

import argparse
from fractions import Fraction as F
from functools import cache
import json
from math import hypot, inf, isclose, nextafter, sqrt
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from model.candidate_laws.interval_echo_partition_repair.published_strings.archive import address_lotteries, posted_words
from .body import gravitational_mass

HERE = Path(__file__).resolve().parent
INNER, OUTER = 2, 4
TICK = 2.0**-44


def velocity(mu, r):
    return -sqrt(2 * mu / r)


def stationary(r, mu, mass=1):
    """A held clock has zero radial coordinate velocity in the moving carrier."""
    v = velocity(mu, r)
    if not -1 < v <= 0 or mass <= 0:
        raise ValueError("Stationary clocks require positive mass and an exterior position")
    p = -mass * v / sqrt((1 + v) * (1 - v))
    return p, mass / hypot(mass, p)


def rk4(field, start, length, steps, state):
    step = length / steps
    def advance(x, y, scale):
        return tuple(a + scale*b for a, b in zip(x, y, strict=True))
    for i in range(steps):
        x = start + i * step
        a = field(x, state)
        b = field(x + step/2, advance(state, a, step/2))
        c = field(x + step/2, advance(state, b, step/2))
        d = field(x + step, advance(state, c, step))
        state = tuple(v + step * (w + 2*y + 2*z + q) / 6
                      for v, w, y, z, q in zip(state, a, b, c, d, strict=True))
    return state


def light(mu, start, end, momentum, steps, launch_time=0.0):
    """Integrate the Hamilton clock and momentum to a receiver's measured radius."""
    def derivative(r, state):
        _, p = state
        v = velocity(mu, r)
        speed = v + (1 if p > 0 else -1)
        return 1 / speed, v * p / (2 * r * speed)
    if steps <= 0 or momentum == 0 or (end - start) * momentum <= 0:
        raise ValueError("A null pulse needs a directed path to a distinct receiver")
    if min(start, end) <= 2 * mu:
        raise ValueError("Both stationary null endpoints must lie outside the horizon")
    return rk4(derivative, start, end - start, steps, (launch_time, momentum))


def released(mu, mass, duration, steps):
    """Advance r,p,t by the massive clock's own elapsed proper time."""
    p0, _ = stationary(INNER, mu, mass)
    def derivative(_, state):
        r, p, _ = state
        v = velocity(mu, r)
        gamma = hypot(mass, p) / mass
        return v * gamma + p / mass, gamma * v * p / (2 * r), gamma
    return rk4(derivative, 0.0, duration, steps, (float(INNER), p0, 0.0))


def number(iv, q):
    q = F(q)
    return iv.mpf(q.numerator) / q.denominator


def outward_floats(lo, hi):
    return [nextafter(float(lo.a), -inf), nextafter(float(hi.b), inf)]


@cache
def ray_enclosures(mu):
    if mu == 0:
        return [[OUTER - INNER, OUTER - INNER]] * 2
    iv = MPIntervalContext()
    iv.dps = 50
    b = iv.sqrt(2 * number(iv, mu))
    a0, a1 = (iv.sqrt(r) for r in (INNER, OUTER))
    out = []
    for sign in (1, -1):
        x = OUTER - INNER + sign * 2*b*(a1-a0) + 2*b*b*iv.ln((a1-sign*b)/(a0-sign*b))
        out.append(outward_floats(x, x))
    return out


@cache
def release_enclosure(mu, duration):
    if mu == 0 or duration == 0:
        return [float(INNER)] * 2
    iv = MPIntervalContext()
    iv.dps = 55
    phase = number(iv, duration) * iv.sqrt(2 * number(iv, mu) / INNER**3)
    lo, hi = F(0), F(1, 2)
    if phase.b >= (number(iv, hi) + iv.sin(2*number(iv, hi))/2).a:
        raise ValueError("The registered probe appointment must lie on the enclosed arc")
    for _ in range(84):
        mid = (lo + hi) / 2
        m = number(iv, mid)
        f = m + iv.sin(2*m) / 2
        if f.b < phase.a:
            lo = mid
        elif f.a > phase.b:
            hi = mid
        else:
            break
    return outward_floats(INNER * iv.cos(number(iv, hi))**2, INNER * iv.cos(number(iv, lo))**2)


def distance_from_exact(value, bounds):
    # One outward binary64 rounding covers the final subtraction.
    largest = max(abs(value - edge) for edge in bounds)
    return nextafter(largest, inf) if largest else 0.0


def measured_mu_from_redshift(ratio):
    return (1 - ratio*ratio) / (2 * (1/INNER - ratio*ratio/OUTER))


def source_counts():
    words = posted_words()
    return {name: (sum(w*a for w, (a, _) in zip(p, words, strict=True)),
                   sum(w*a*b for w, (a, b) in zip(p, words, strict=True)))
            for name, p in address_lotteries(words).items()}


def survey(mu_light, mu_matter):
    _, inner_rate = stationary(INNER, mu_light)
    _, outer_rate = stationary(OUTER, mu_light)
    truth_out, truth_back = ray_enclosures(mu_light)
    ray_errors = []
    for steps in (4, 8, 16, 32, 64):
        flight, p1 = light(mu_light, INNER, OUTER, 1.0, steps)
        # A stationary mirror returns the same photon frequency in its rest frame.
        surface_v = velocity(mu_light, OUTER)
        pback = -p1 * (1 + surface_v) / (1 - surface_v)
        inward, p0 = light(mu_light, OUTER, INNER, pback, steps)
        q = -velocity(mu_light, INNER)
        h = (OUTER - INNER) / steps
        truncation = (OUTER - INNER) * h**4 * 24*q / (2880 * INNER**4 * (1-q)**5)
        ray_errors.append({
            "panels": steps, "outbound": flight, "inbound": inward, "exact_arithmetic_time_error_bound": truncation,
            "outbound_error_upper": distance_from_exact(flight, truth_out),
            "inbound_error_upper": distance_from_exact(inward, truth_back),
        })
    emitted = [0, round(1 / TICK)]
    second, _ = light(mu_light, INNER, OUTER, 3.0, steps, launch_time=1 / inner_rate)
    received = [round(outer_rate * t / TICK) for t in (flight, second)]
    radar = round((flight + inward) * inner_rate / TICK)
    ratio = (emitted[1] - emitted[0]) / (received[1] - received[0])
    # Frequency of a wave covector measured by the stationary observer.
    outer_energy = p1 * (1 + surface_v) / outer_rate
    inner_energy = (1 + velocity(mu_light, INNER)) / inner_rate
    returned_energy = -p0 * (1 - velocity(mu_light, INNER)) / inner_rate
    releases = {}
    for mass in (1, 3):
        appointment = {}
        for duration in (F(1, 2), F(1, 16), F(1, 32)):
            enclosed = release_enclosure(mu_matter, duration)
            levels = []
            for steps in (1, 2, 4, 8, 16):
                r, p, t = released(mu_matter, mass, float(duration), steps)
                levels.append(distance_from_exact(r, enclosed))
            observed_r = round(r / TICK)
            appointment[str(duration)] = {
                "radius_ticks": observed_r, "action_time_ticks": round(t / TICK),
                "probe_proper_ticks": round(float(duration) / TICK), "radius_error_upper_by_refinement": levels,
                "inferred_radial_acceleration": 2 * (observed_r * TICK - INNER) / float(duration)**2,
            }
        releases[str(mass)] = appointment
    raw = releases["1"]
    a = (4 * raw["1/32"]["inferred_radial_acceleration"] - raw["1/16"]["inferred_radial_acceleration"]) / 3
    return {"emitted_inner_ticks": emitted, "received_outer_ticks": received, "echo_inner_ticks": radar,
            "redshift_ratio": ratio, "spectroscopic_ratio": outer_energy / inner_energy,
            "mirror_recovered_frequency_relative_error": abs(returned_energy/inner_energy - 1),
            "signed_rays": ray_errors, "released_clocks": releases,
            "mu_from_clock_spacing": measured_mu_from_redshift(ratio),
            "mu_from_two_release_intervals": -INNER**2 * a}


def run():
    counts = source_counts()
    ordinary, _ = counts["ordinary"]
    matter, collar = counts["held"]
    total = matter + collar
    sources = {"blank": (F(0), F(0)), "ordinary": (ordinary, ordinary), "held": (total, total),
               "both_collar_blind": (matter, matter), "light_alone_collar_blind": (matter, total)}
    rows = {name: survey(gravitational_mass(ml), gravitational_mass(mm)) for name, (ml, mm) in sources.items()}
    fitted = {meter: rows["ordinary"][meter] / float(ordinary)
              for meter in ("mu_from_clock_spacing", "mu_from_two_release_intervals")}
    forecast = {}
    for meter, val in fitted.items():
        mu = val * float(total)
        forecast[meter] = {
            "mu": mu, "redshift_ratio": sqrt((1 - 2*mu/INNER) / (1 - 2*mu/OUTER)),
            "initial_radial_acceleration": -mu / INNER**2,
        }
    return {"schema": "oph.spherical.exterior_hamilton_reader.v1", "proper_tick": TICK,
            "areal_radii": [INNER, OUTER], "release_panels": [1, 2, 4, 8, 16],
            "ordinary_and_collar_counts": {k: [str(v) for v in pair] for k, pair in counts.items()},
            "readouts": rows, "ordinary_gain_by_meter": fitted,
            "held_forecasts_from_ordinary": forecast}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    path = HERE / "exterior_reader.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        check_archive(data, json.loads(path.read_text()))
    print("Verified canonical ray, mirror and released-clock records.")


def check_archive(actual, saved, path=()):
    """Verify every emitted field, allowing for binary64 library and tick rounding."""
    if isinstance(actual, dict) and isinstance(saved, dict) and actual.keys() == saved.keys():
        for key in actual:
            check_archive(actual[key], saved[key], path + (key,))
        return
    if isinstance(actual, list) and isinstance(saved, list) and len(actual) == len(saved):
        for new, old in zip(actual, saved, strict=True):
            check_archive(new, old, path)
        return
    if isinstance(actual, int) and isinstance(saved, int) and path[-1].endswith("_ticks"):
        if abs(actual - saved) <= 8:
            return
    elif isinstance(actual, float) and isinstance(saved, (float, int)):
        # The release inversion amplifies last-tick noise by O(1/duration²).
        release = "mu_from_two_release_intervals" in path or path[-1] == "inferred_radial_acceleration"
        clock = "mu_from_clock_spacing" in path or path[-1] == "redshift_ratio"
        tolerance = 1e-8 if release else 5e-12 if clock else 2e-14
        if isclose(actual, saved, rel_tol=1e-13, abs_tol=tolerance):
            return
    elif actual == saved:
        return
    raise ArithmeticError("Exterior record disagrees with archive at " + ".".join(path))


if __name__ == "__main__":
    main()
