"""Independent angular screens, affine clocks and dust-action stress on one outgoing ray."""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from math import prod
from pathlib import Path

from .body import Body, archive_bodies, release_calibration
from .survey import mass_tally


HERE = Path(__file__).resolve().parent
SCALES = (32, 64, 128)
END = F(1, 4)
CENTER = F(1, 8)
CHORDS = (F(14, 25), F(10, 13), F(16, 17))
CLOCK_STEP = F(1, 4)


def decimal(x):
    x = F(x)
    return Decimal(x.numerator) / Decimal(x.denominator)


def mass(x):
    return F(x).limit_denominator(66)


def acquire(body, den, include_collar=True):
    """Only the emitter uses the body dynamics to generate the finite detector tapes."""
    tick, phase_step = F(1, den**10), F(1, den**2)
    with localcontext() as ctx:
        ctx.prec = 65
        h = decimal(mass(body.gravitating_count) / 2).sqrt()
        quantum, step = decimal(tick), decimal(phase_step)

        def stamp(value):
            return int((value / quantum).to_integral_value())

        def proper(phase):
            return phase * (1 - h * phase / 2 + (h * phase)**2 / 12)

        clock = [stamp(proper(decimal(k * CLOCK_STEP))) for k in (-1, 0, 1, 2)]
        events = []
        for i in range(den // 4 + 1):
            target = decimal(F(i, den))
            phase = 2 * (1 - ((1 - 5 * h * target / 2).ln() / 5).exp()) / h if h else target
            phase_tick = stamp(phase)
            phase = phase_tick * quantum
            scale = (1 - h * phase / 2)**2
            screens = [stamp(decimal(chord) * scale * phase) for chord in CHORDS]
            shells = [stamp(scale / 4), stamp(scale / 2)]
            frequency = [stamp(proper(phase + k * step)) for k in (-2, -1, 1, 2)]
            events.append([phase_tick, *screens, *shells, *frequency])
    tickets = [mass_tally(mass(body.ordinary), den), mass_tally(mass(body.collar), den)]
    return {"clock_unit": str(tick), "metered_stocks": 2 if include_collar else 1,
            "shell_mass_tickets": [list(row) for row in tickets], "clock_calibration": clock, "events": events}


def rate_polynomial(clock, tick):
    """Read dτ/dη as a quadratic from four conformal-phase/clock appointments."""
    minus, zero, plus, twice = [k * tick for k in clock]
    s = CLOCK_STEP
    symmetric = (plus - minus) / (2 * s)
    quadratic = (plus + minus - 2 * zero) / (2 * s**2)
    cubic = (twice - zero - 2 * s * symmetric - 4 * s**2 * quadratic) / (6 * s**3)
    return symmetric - cubic * s**2, 2 * quadratic, 3 * cubic


def affine_clock(phase, coefficients):
    return sum(a * b * phase**(i + j + 1) / (i + j + 1)
               for i, a in enumerate(coefficients) for j, b in enumerate(coefficients)) / coefficients[0]


def gain_from_release(record, ordinary_mass):
    dt, tick = F(record["duration"]), F(record["tick"])
    initial = record["start"] * tick
    released = record["released"] * tick
    difference = initial - released
    factor = 2 / (ordinary_mass * dt**2)
    gain = factor * difference * initial**2
    # Twice integrated acceleration is a positive weighted mean of GM/r(t)².
    # The two observed endpoint radii enclose the whole monotone ordinary release.
    lo = factor * (difference - tick / 2) * (released - tick / 2)**2
    hi = factor * (difference + tick / 2) * initial**2
    error = max(gain - lo, hi - gain) + F(1, 10**25)
    return gain, error


def material_read(row, tick, den, initial_rate, matter, gain, gain_error):
    inner, outer = row[4] * tick, row[5] * tick
    a, b, c, d = [k * tick for k in row[6:10]]
    phase_step = F(1, den**2)
    local_rate = (a - 8 * b + 8 * c - d) / (12 * phase_step)
    local_error = 3 * tick / (4 * phase_step)
    initial_error = 10 * tick
    half = tick / 2

    def source(g, inside, outside, initial, local):
        return 6 * g * matter * (initial / local)**2 / (outside**3 - inside**3)

    estimate = source(gain, inner, outer, initial_rate, local_rate)
    lo = source(gain - gain_error, inner - half, outer + half, initial_rate - initial_error, local_rate + local_error)
    hi = source(gain + gain_error, inner + half, outer - half, initial_rate + initial_error, local_rate - local_error)
    # Convert the rounded phase setting to the nominal affine event.
    shift_error = 30 * tick
    return estimate, max(F(0), lo - shift_error), hi + shift_error, initial_rate / local_rate


def derivative_weight(power, n):
    # q=(1-5hλ/2)^(1/5) >= 0.9; each derivative lowers the exponent by 5.
    return abs(prod(power - 5 * j for j in range(n))) * F(9, 10)**(power - 5 * n)


def optical_stencil_bound(spacing, tick, screen):
    h = F(11, 20)
    # On [1/8-1/32,1/8+1/32], q>=.95 and D''''=−h³(624−504q)/(8q^18).
    q = F(19, 20)
    fourth = h**3 * (624 - 504 * q) / (8 * q**18)
    screen_error = 3 * tick
    return (2 * spacing**2 * fourth / 12 + 8 * screen_error / spacing**2 + 3 * screen_error) / screen


def integral_remainder(spacing):
    # D(λ) κTkk = (12 G M/h)(q^-8-q^-7); n>=1 removes the apparent 1/h.
    def bound(n):
        return 12 * F(13, 50) * F(3, 5) * F(11, 20)**(n - 1) / 2**n * (
            derivative_weight(-8, n) + derivative_weight(-7, n))

    return END * spacing**4 * (END * bound(4) + 4 * bound(3)) / 180


def decode(packet, den, release, ordinary_mass):
    """The decoder receives rod, clock and matter records plus the ordinary calibration."""
    ds = F(1, den)
    tick = F(packet["clock_unit"])
    rows = packet["events"]
    if len(rows) != END * den + 1 or any(len(row) != 10 or any(type(k) is not int for k in row) for row in rows):
        raise ValueError("invalid angular screen or clock tape")
    coefficients = rate_polynomial(packet["clock_calibration"], tick)
    tickets = packet["shell_mass_tickets"][:packet["metered_stocks"]]
    shell_mass = sum(F(outer - inner, n) for n, inner, outer in tickets)
    gain, gain_error = gain_from_release(release, ordinary_mass)
    if gain + gain_error > F(13, 50) or abs(coefficients[1]) + 64 * tick > F(11, 20):
        raise ValueError("measured ordinary coupling or clock slope exceeds the finite error domain")
    screens, work, affine_errors = [], [], []
    shear = F(0)
    clock_prediction_error = F(0)
    for i, row in enumerate(rows):
        phase = row[0] * tick
        if phase > F(1, 3) or phase < 0:
            raise ValueError("beam left the calibrated interior clock range")
        affine_errors.append(abs(affine_clock(phase, coefficients) - i * ds))
        measured = [row[j + 1] * tick / chord for j, chord in enumerate(CHORDS)]
        shear = max(shear, max(measured) - min(measured))
        screens.append(sum(measured) / len(measured))
        reading = material_read(row, tick, den, coefficients[0], shell_mass, gain, gain_error)
        work.append(reading)
        predicted_rate = sum(value * phase**j for j, value in enumerate(coefficients))
        clock_prediction_error = max(clock_prediction_error, abs(predicted_rate - coefficients[0] / reading[3]))
    if max(affine_errors) > 120 * tick or shear > 2 * tick or clock_prediction_error > 100 * tick + 3 * den**2 * tick / 4:
        raise ValueError("independent affine, frequency or transverse screen calibration failed")
    mid = int(CENTER * den)
    d = screens[mid]
    focusing = -2 * (screens[mid-1] - 2 * d + screens[mid+1]) / (d * ds**2)
    material, lo, hi, omega = work[mid]
    optical_error = optical_stencil_bound(ds, tick, d)
    local_bound = optical_error + max(material - lo, hi - material)

    integrals = [F(0), F(0), F(0)]
    for i, (screen, values) in enumerate(zip(screens, work, strict=True)):
        weight = 1 if i in (0, len(rows)-1) else 4 if i % 2 else 2
        factor = ds * weight * (END - i * ds) / 3
        integrals[0] += factor * screen * values[0]
        integrals[1] += factor * max(F(0), screen - 3 * tick) * values[1]
        integrals[2] += factor * (screen + 3 * tick) * values[2]
    endpoint = 2 * (END - screens[-1])
    ibound = integral_remainder(ds) + 6 * tick + max(integrals[0] - integrals[1], integrals[2] - integrals[0])
    return {key: float(value) for key, value in {
        "ordinary_release_gain": gain, "shell_rest_mass": shell_mass, "midpoint_frequency": omega,
        "largest_affine_clock_error": max(affine_errors), "largest_transverse_screen_difference": shear,
        "midpoint_focusing": focusing, "midpoint_work_forecast": material, "signed_point_residual": focusing - material,
        "point_error_bound": local_bound, "endpoint_angular_distance": screens[-1], "twice_affine_deficit": endpoint,
        "weighted_work_integral": integrals[0], "signed_integral_residual": endpoint - integrals[0],
        "integral_error_bound": ibound}.items()}


def measure():
    bodies = archive_bodies()
    ordinary_mass = mass(bodies["ordinary"].ordinary)
    held = bodies["held"]
    arms = {**bodies, "half_collar_force": Body(held.ordinary, held.collar, 0.5), "blind_collar_work": held}
    calibrations = {str(d): release_calibration(ordinary_mass, d) for d in SCALES}
    stocks = {}
    for name, body in arms.items():
        stocks[name] = {}
        for d in SCALES:
            raw = acquire(body, d, include_collar=name != "blind_collar_work")
            stocks[name][str(d)] = {"acquisition": raw, "readings": decode(raw, d, calibrations[str(d)], ordinary_mass)}
    signed = {}
    for d in SCALES:
        ordinary, combined = (stocks[name][str(d)]["readings"] for name in ("ordinary", "held"))
        signed[str(d)] = {k: combined[k] - ordinary[k] for k in (
            "midpoint_focusing", "midpoint_work_forecast", "twice_affine_deficit", "weighted_work_integral")}
    return {"schema": "oph.spherical_dust.null_beam.v1", "affine_midpoint": str(CENTER), "affine_endpoint": str(END),
            "screen_chords": [str(c) for c in CHORDS], "ordinary_release": calibrations, "stocks": stocks,
            "held_minus_ordinary": signed}


def verify():
    rebuilt = measure()
    if rebuilt != json.loads((HERE / "null_beam.json").read_text()):
        raise ValueError("spherical null-beam records differ from the independently recomputed optical/work audit")
    for name, arms in rebuilt["stocks"].items():
        for den, arm in arms.items():
            r = arm["readings"]
            inside = [abs(r["signed_point_residual"]) <= r["point_error_bound"],
                      abs(r["signed_integral_residual"]) <= r["integral_error_bound"]]
            if name in ("blank", "ordinary", "held") and not all(inside):
                raise ValueError(f"matched dust null response outside finite error bounds at {name}/{den}")
            if name in ("half_collar_force", "blind_collar_work") and int(den) == SCALES[-1] and any(inside):
                raise ValueError(f"null response failed to distinguish held control {name}")
    return rebuilt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    if parser.parse_args().write:
        (HERE / "null_beam.json").write_text(json.dumps(measure(), indent=2) + "\n")
    else:
        result = verify()
        print(json.dumps({k: v[str(SCALES[-1])]["readings"] for k, v in result["stocks"].items()}, indent=2))
