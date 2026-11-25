"""Held station clock identities and a five-to-eight-station geometry comparison."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from . import flight
from .. import live_sphere, live_survey

HERE = Path(__file__).parent


def identities(minus, plus):
    """Necessary cycle and tetrad identities, with a multiplicative error bound."""
    n = len(plus)
    cycles = {}
    for tag, directed in (("minus", minus), ("plus", plus)):
        z = np.log(directed + np.eye(n))
        skew = z - z.T
        # At least one held station in each tested loop.
        trials = [
            (abs(skew[i, j] + skew[j, k] + skew[k, i]), (i, j, k))
            for i, j, k in combinations(range(n), 3)
            if k >= 5
        ]
        gap, loop = max(trials)
        cycles[tag] = dict(
            log_cycle_error=float(gap),
            ports=[flight.STATIONS[i] for i in loop],
            per_reading_log_error_lower=float(gap / 6),
        )
    z = np.log((plus + plus.T + np.eye(n)) / (minus + minus.T + np.eye(n)))
    trials = []
    for i, j, k, last in combinations(range(n), 4):
        if last < 5:
            continue
        for a, b, c, d in ((i, j, k, last), (i, k, last, j), (i, last, j, k)):
            mismatch = abs(z[a, b] + z[c, d] - z[a, c] - z[b, d])
            trials.append((mismatch, (a, b, c, d)))
    gap, ports = max(trials)
    return dict(
        oriented=cycles,
        reciprocal_species=dict(
            log_tetrad_error=float(gap),
            ports=[flight.STATIONS[i] for i in ports],
            per_reading_log_error_lower=float(gap / 8),
        ),
    )


def analyze(record, divisor, continuous):
    key = "continuous_panels" if continuous else "panels"
    rows = {
        p["charge"]: p["states"]["reference"]
        for p in record[key]
        if p["divisor"] == divisor
    }
    clocks = {
        name: identities(*(np.asarray(rows[q][name]) for q in (-1, 1)))
        for name in ("capped_range", "first_capped_duration")
    }
    d = flight.pooled(record, divisor, continuous=continuous)["reference"]
    smooth = live_sphere.sphere_fit(d)
    area, volume, angle = live_survey.star(d[:5, :5])
    return dict(
        divisor=divisor,
        continuous=continuous,
        clock_tests=clocks,
        first_five=dict(
            area=float(area), volume=float(volume), apex_solid_angle=float(angle)
        ),
        fitted_curvature=smooth["curvature"],
        held_distance_error=smooth["maximum_held_distance_error"],
        six_station_obstruction=smooth["all_radii"],
    )


def run():
    record = json.loads((HERE / "flight.json").read_text())
    return [analyze(record, d, c) for c in (False, True) for d in flight.DIVISORS]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "comparison.json").write_text(data)
    else:
        print(data, end="")
