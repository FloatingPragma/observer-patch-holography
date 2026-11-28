"""Trace angular address, mirror echo and clock records on one spherical source."""

import argparse
from fractions import Fraction
import json
from math import asin, cos, hypot, pi, sin, sqrt
from pathlib import Path

from scipy.integrate import solve_ivp

from model.candidate_laws.interval_echo_published_records.spherical_dust import exterior_reader as meter
from model.candidate_laws.interval_echo_published_records.spherical_dust.body import gravitational_mass
from .forecast import calibrated_mu, calibration, forecast

RADIUS = meter.OUTER
TICK = meter.TICK
SINES = (0.5, 0.75)
ARCHIVE = Path(__file__).with_name("results.json")


def angular_factor(mu, r, change):
    if not change or r >= RADIUS:
        return 1., 0.
    x = 1 - r/RADIUS
    return 1 + 2*mu*x*x/r, -2*mu*(x*x/r**2 + 2*x/(RADIUS*r))


def arc(mu, momentum, angular, phi0, change):
    def rates(_, y):
        r, _, p = y
        c, dc = angular_factor(mu, r, change)
        v = meter.velocity(mu, r)
        tangential = c*angular/r
        w = hypot(p, tangential)
        return [v+p/w, c*c*angular/(w*r*r),
                v*p/(2*r) + tangential*tangential*(1/r-dc/c)/w]

    def receive(t, y):
        return y[0] - RADIUS

    def leave_exterior(t, y):
        return y[0] - 1

    def turn(t, y):
        return rates(t, y)[0]

    receive.direction, receive.terminal = 1, True
    leave_exterior.direction, leave_exterior.terminal = -1, True
    turn.direction = 1
    solution = solve_ivp(rates, (0, 32), [RADIUS, phi0, momentum], method="DOP853",
                         rtol=5e-13, atol=2e-14, max_step=0.25,
                         events=(receive, leave_exterior, turn))
    if not solution.success or len(solution.t_events[0]) != 1 or len(solution.t_events[2]) != 1:
        return None
    return solution.t_events[0][0], solution.y_events[0][0], solution.y_events[2][0][0]


def observe(mu, sine, change=False):
    v = meter.velocity(mu, RADIUS)
    _, clock = meter.stationary(RADIUS, mu)
    incoming_p = (-sqrt(1-sine*sine)-v)/clock
    angular = RADIUS*sine
    outward = arc(mu, incoming_p, angular, 0., change)
    if outward is None:
        return {"event": "first leg missed"}
    t1, (_, phi, p), closest = outward
    w = hypot(p, angular/RADIUS)
    # Boost the received covector to the static mirror tetrad. Reverse its local
    # spatial components and boost back, preserving the frequency it actually sees.
    frequency = (w+v*p)/clock
    radial = (p+v*w)/clock
    seen_sine = angular/(RADIUS*frequency)
    return_p = (-radial-v*frequency)/clock
    back = arc(mu, return_p, -angular, phi, change)
    if back is None:
        return {"event": "mirror leg missed"}
    t2, (_, home, _), _ = back
    address = round(phi/TICK)
    local_in = round(asin(sine)/TICK)
    local_out = round(asin(seen_sine)/TICK)
    echo = round((t1+t2)*clock/TICK)
    return {"event": "returned", "mirror_address_ticks": address, "echo_address_ticks": round(home/TICK),
            "emission_angle_ticks": local_in, "reception_angle_ticks": local_out,
            "echo_clock_ticks": echo, "bending": (address+local_in+local_out)*TICK-pi,
            "round_trip": echo*TICK,
            "minimum_radius_diagnostic": closest}


def paired_observation(mu, sine, change=False):
    reading = observe(mu, sine, change)
    if reading["event"] != "returned":
        return reading
    phi = reading["mirror_address_ticks"] * TICK
    # The blank apparatus aims at the same address, not the same launch direction.
    blank = observe(0., cos(phi/2))
    if blank["event"] != "returned":
        return {"event": "matched blank did not return"}
    offset = blank["mirror_address_ticks"] - reading["mirror_address_ticks"]
    if abs(offset*TICK) > 1e-10 or abs(blank["echo_address_ticks"]*TICK) > 1e-10:
        raise AssertionError("The blank missed the enrolled address aperture")
    if abs(blank["round_trip"] - 4*RADIUS*sin(phi/2)) > 1e-10:
        raise AssertionError("The flat radar disagrees with its independently known chord")
    reading["matched_blank_mirror_address_ticks"] = blank["mirror_address_ticks"]
    reading["matched_blank_echo_clock_ticks"] = blank["echo_clock_ticks"]
    reading["delay"] = (reading["echo_clock_ticks"] - blank["echo_clock_ticks"])*TICK
    return reading


def contains(bounds, value):
    return bounds[0] <= value <= bounds[1]


def run():
    baseline = meter.run()
    counts = meter.source_counts()
    ordinary = counts["ordinary"][0]
    cal = calibration(baseline["readouts"]["ordinary"], ordinary)
    rows = {}
    for name, mass in [("blank", Fraction(0)), ("ordinary", ordinary), ("held", sum(counts["held"]))]:
        mu = float(gravitational_mass(mass))
        bounds = calibrated_mu(cal["gain"], mass)
        radial = baseline["readouts"][name]
        shots = []
        for sine in SINES:
            predicted = forecast(bounds, sine)
            reading = paired_observation(mu, sine)
            if reading["event"] != "returned":
                raise AssertionError(f"{name}, {sine}: {reading}")
            for key in ("bending", "delay", "round_trip"):
                if not contains(predicted[key], reading[key]):
                    raise AssertionError((name, sine, key, predicted[key], reading[key]))
            # Clock-ratio quantization rounds two independent receiver dates.
            clock_window = predicted["clock_ratio_2_to_4"]
            if not clock_window[0]-2*TICK <= radial["redshift_ratio"] <= clock_window[1]+2*TICK:
                raise AssertionError((name, clock_window, radial["redshift_ratio"]))
            if abs(reading["echo_address_ticks"]*TICK) > 1e-9:
                raise AssertionError("Return missed its emitter")
            shot = {"local_sine": sine, "reading": reading, "forecast": predicted}
            if name == "held":
                altered = paired_observation(mu, sine, change=True)
                if altered["event"] != "returned" or any(contains(predicted[key], altered[key])
                        for key in ("bending", "delay")):
                    raise AssertionError("Altered light law passed the common-action forecast")
                shot["angular_light_control"] = altered
            shots.append(shot)
        rows[name] = {"stock": str(mass), "emitted_inner_ticks": radial["emitted_inner_ticks"],
                      "received_outer_ticks": radial["received_outer_ticks"],
                      "clock_ratio_2_to_4": radial["redshift_ratio"], "shots": shots}
    for i in range(len(SINES)):
        for key in ("bending", "delay"):
            values = [rows[n]["shots"][i]["reading"][key] for n in ("blank", "ordinary", "held")]
            if abs(values[0]) > 1e-10 or not values[0] < values[1] < values[2]:
                raise AssertionError((i, key, values))
    return {"schema": "oph.row27.spherical_exterior.v1", "angular_and_time_tick": TICK,
            "ordinary_release_calibration": cal, "sources": rows}


def check(actual, saved, path=""):
    if isinstance(actual, dict) and isinstance(saved, dict) and actual.keys() == saved.keys():
        for key in actual:
            check(actual[key], saved[key], f"{path}.{key}")
    elif isinstance(actual, list) and isinstance(saved, list) and len(actual) == len(saved):
        for a, b in zip(actual, saved):
            check(a, b, path)
    elif isinstance(actual, (float, int)) and isinstance(saved, (float, int)):
        tolerance = 128 if path.endswith("_ticks") else 2e-9
        if abs(actual-saved) > tolerance:
            raise AssertionError((path, actual, saved))
    elif actual != saved:
        raise AssertionError((path, actual, saved))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    output = run()
    if parser.parse_args().write:
        ARCHIVE.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    else:
        check(output, json.loads(ARCHIVE.read_text()))
    print("Verified six spherical lens, radar and clock comparisons and two angular-light controls.")
