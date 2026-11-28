"""Infer compact-source gravitational length from anonymous clock differences."""

import argparse
import json
from fractions import Fraction
from math import log, log1p, sqrt
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from model.candidate_laws.interval_echo_published_records.spherical_dust import (
    exterior_reader as lab,
)
from model.candidate_laws.interval_echo_published_records.spherical_dust.clock_pairs import (
    census,
)

HERE = Path(__file__).resolve().parent
TICK = 2.0**-36


def acquire(mu):
    inner = 2
    _, lapse = lab.stationary(inner, mu)
    links = []
    for j, outer in enumerate((3, 4)):
        _, remote = lab.stationary(outer, mu)
        t1, _ = lab.light(mu, inner, outer, 1.0, 512)
        t2, _ = lab.light(mu, inner, outer, 3.0, 512, launch_time=1 / lapse)
        received, _ = lab.light(mu, outer, inner, -1.0, 512, launch_time=t1)
        # Private clock origins never leave their respective clocks.
        send = [4321, 4321 + round(1 / TICK)]
        arrival = [(-7011 + 1200 * j) + round(t * remote / TICK) for t in (t1, t2)]
        back = 4321 + round(lapse * received / TICK)
        # Independent stationary exterior integration, used only to bound the
        # forward clock error before digitizing; these radii are absent at decode.
        exact = (
            2
            * lapse
            * (outer - inner + 2 * mu * log((outer - 2 * mu) / (inner - 2 * mu)))
        )
        assert abs(lapse * received - exact) < TICK
        assert abs(remote * (t2 - t1) - remote / lapse) < TICK
        links.append(dict(emission=send, reception=arrival, return_tick=back))
    return dict(tick=TICK, links=links)


def spans(tape):
    out = []
    for link in tape["links"]:
        a, b = link["emission"]
        c, d = link["reception"]
        local = b - a
        remote = d - c
        echo = link["return_tick"] - a
        if min(local, remote, echo) <= 0:
            raise ValueError("Nonpositive measured duration")
        out.append((local, remote, echo))
    return out


def integral(x, y):
    # Integral from 1 to y of du/[u(1-xu)^2]. The dimensionless lapse at
    # the emitter is sqrt(x), and y is squared reception/emission spacing.
    return (
        log(y)
        + log1p(x * (y - 1) / (1 - x * y))
        + x * (y - 1) / ((1 - x * y) * (1 - x))
    )


def solution(y1, y2, ratio):
    if not 1 < y1 < y2 or ratio <= log(y2) / log(y1):
        raise ValueError("No resolved positive ordered exterior")

    def f(x):
        return integral(x, y2) / integral(x, y1) - ratio

    return brentq(f, 0.0, (1 - 1e-12) / y2, xtol=1e-15)


def decode(tape):
    points = spans(tape)
    if len(points) != 2:
        raise ValueError("Two separately recorded mirror links are required")
    tick = tape["tick"]
    if not np.isfinite(tick) or tick <= 0:
        raise ValueError("Invalid proper-time tick")
    (a, b, t), (c, d, v) = points
    # Allow four ticks in every measured duration, covering both rounding and
    # the independent forward check. The enclosures below use floating roots.
    e = 4
    if min(a, b, c, d, t, v) <= e:
        raise ValueError("Insufficient timing resolution")
    y, Y = (b / a) ** 2, (d / c) ** 2
    ylo, yhi = ((b - e) / (a + e)) ** 2, ((b + e) / (a - e)) ** 2
    Ylo, Yhi = ((d - e) / (c + e)) ** 2, ((d + e) / (c - e)) ** 2
    if ylo <= 1 or yhi >= Ylo:
        return dict(
            resolved=False, reason="lapse gradient unresolved at this clock precision"
        )
    ratio = v / t
    x = solution(y, Y, ratio)
    xl = solution(ylo, Yhi, (v - e) / (t + e)) - 1e-13
    xh = solution(yhi, Ylo, (v + e) / (t - e)) + 1e-13
    if not 0 < xl <= x <= xh < 1 / Yhi:
        raise ValueError("Clock bounds cross the static exterior domain")

    def mass(tt, xx, yy):
        return tt * tick / (4 * sqrt(xx) * integral(xx, yy))

    mu = mass(t, x, y)
    low, high = mass(t - e, xh, yhi), mass(t + e, xl, ylo)
    radii = [2 * mu / (1 - x * z) for z in (1, y, Y)]
    ranges = [
        [2 * low / (1 - xl * z0), 2 * high / (1 - xh * z1)]
        for z0, z1 in ((1, 1), (ylo, yhi), (Ylo, Yhi))
    ]
    pairs = []
    for xx in ((3 * x) / 4 + x * x * y / 4, (3 * x + 1 / y) / 4):
        mm = mass(t, xx, y)
        pairs.append(dict(mu=mm, radii=[2 * mm / (1 - xx * z) for z in (1, y)]))
    return dict(
        resolved=True,
        mu=mu,
        mu_interval=[low, high],
        fitted_radii=radii,
        radius_intervals=ranges,
        two_station_alternatives=pairs,
    )


def probe_ratio(clock, release, tick, rstart):
    # Independent differences of area-radius tickets and probe proper time.
    data = []
    for tag in ("1/16", "1/32"):
        r = release[tag]
        displacement = (r["radius_ticks"] - rstart) * tick
        tau = r["probe_proper_ticks"] * tick
        data.append(2 * displacement / tau**2)
    acceleration = (4 * data[1] - data[0]) / 3
    return -acceleration * clock["fitted_radii"][0] ** 2 / clock["mu"]


def run():
    prior = json.loads((lab.HERE / "exterior_reader.json").read_text())
    counts = lab.source_counts()
    period, words = census.census()
    blank_counts = words["blank"]["word_counts"]
    out = {}
    for label, mat in (
        ("blank", "blank"),
        ("ordinary", "ordinary"),
        ("held", "held"),
        ("light_alone_collar_blind", "ordinary"),
    ):
        t = acquire(float(lab.gravitational_mass(sum(counts[mat]))))
        decoded = decode(t)
        actual = words["held" if label == "light_alone_collar_blind" else label]
        frequencies = [Fraction(z, period) for z in actual["word_counts"]]
        reference = [Fraction(z, period) for z in blank_counts]
        assert (
            frequencies[0] + frequencies[1]
            == frequencies[0] + frequencies[2]
            == Fraction(1, 2)
        )
        modular = sum(
            float(p - q) * (-log(float(q)))
            for p, q in zip(frequencies, reference, strict=True)
        )
        if decoded["resolved"]:
            release = prior["readouts"][label]["released_clocks"]
            tick = prior["proper_tick"]
            # Original radius instrument at the release event, used only for
            # its subsequent displacement. It never enters the clock inversion.
            r0 = round(prior["areal_radii"][0] / tick)
            decoded["relative_free_fall"] = {
                key: probe_ratio(decoded, row, tick, r0) for key, row in release.items()
            }
        out[label] = dict(
            clock_tape=t,
            inferred=decoded,
            word_counts=actual["word_counts"],
            modular_change_nats=modular,
            one_end_modular_change_nats=0.0,
        )
    ordinary = float(counts["ordinary"][0])
    gain = out["ordinary"]["inferred"]["mu"] / ordinary
    modular_gain = (
        out["ordinary"]["inferred"]["mu"] / out["ordinary"]["modular_change_nats"]
    )
    return dict(
        tick=TICK,
        ordinary_reference_count=ordinary,
        calibrated_gain=gain,
        mu_forecast_from_single_pair={
            k: modular_gain * v["modular_change_nats"] for k, v in out.items()
        },
        held_active_mass_in_reference_units=out["held"]["inferred"]["mu"] / gain,
        word_census_denominator=period,
        sources=out,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = run()
    if parser.parse_args().write:
        (HERE / "exterior.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
    for name, record in result["sources"].items():
        print(name, record["inferred"], record["modular_change_nats"])
