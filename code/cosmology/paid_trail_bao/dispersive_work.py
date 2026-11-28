"""Compare the below-band clock profile with four outer orbital work intervals."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np

from . import occupied_work as orbit
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    dispersive_clock as lab,
)

HERE = orbit.HERE


def pointwise(profile, observations):
    intervals = []
    for r, observed in zip(range(2, 6), observations, strict=True):
        low, high = lab.meter.difference_interval(profile[r], profile[r + 1])
        assert low > 0
        intervals.append(
            [
                float(np.nextafter(observed[0] / high, -np.inf)),
                float(np.nextafter(observed[1] / low, np.inf)),
            ]
        )
    return dict(
        amplitude_intervals=intervals,
        common_amplitude=max(a for a, _ in intervals) <= min(b for _, b in intervals),
    )


def continuum(before, marks, r):
    # With fixed positive source D the mixed diagonal is a sum of positive
    # products of resolvent entries, each increasing as omega moves up.
    spectral_floor = lab.proof.floor_band(before)
    endpoints = {}
    for k in range(-10, 6):
        gap = Fraction(2) ** k
        value, _ = lab.proof.mixed_columns(before, marks, spectral_floor - gap, gap)
        shells = lab.proof.shell_intervals(value, r)
        endpoints[gap] = shells[[2, 4, 6]]
    pieces = []
    gaps = list(endpoints)
    for lo, hi in zip(gaps[:-1], gaps[1:], strict=True):
        near, far = endpoints[lo], endpoints[hi]
        inner_floor = np.nextafter(far[0, 0] - near[1, 1], -np.inf)
        outer_cap = np.nextafter(near[1, 1] - far[2, 0], np.inf)
        assert inner_floor > 0 and outer_cap >= 0
        ratio = float(np.nextafter(outer_cap / inner_floor, np.inf))
        pieces.append(
            dict(
                detuning_below_floor=[str(lo), str(hi)],
                inner_floor=inner_floor,
                outer_upper=outer_cap,
                outer_per_inner_upper=ratio,
            )
        )
    return dict(
        ordinary_floor=str(spectral_floor),
        panels=pieces,
        maximum_ratio=max(p["outer_per_inner_upper"] for p in pieces),
    )


def run():
    clockfile = lab.flight.HERE / "dispersive_clock.json"
    info = json.loads(clockfile.read_bytes())
    board, _, edges, r, _ = lab.flight.chart()
    clamp = info["clamp"]
    radii, observations = orbit.empirical()
    # The atomic/stellar reductions, adopted distance and integration convention
    # enter these measured work windows. They remain fixed throughout the test.
    required = float(
        np.nextafter(
            sum(x[0] for x in observations[2:]) / sum(x[1] for x in observations[:2]),
            -np.inf,
        )
    )
    rows = {}
    for species, t in (("vacancy", 0), ("excess", 2)):
        before, marks = lab.flight.apparatus(board, edges, t, clamp, 1)
        original = {
            k: pointwise(np.array(v["mixed_shell"]), observations)
            for k, v in info["results"][species]["inactive"].items()
        }
        whole = continuum(before, marks, r[np.arange(board.n) != clamp])
        rows[species] = dict(original=original, continuous=whole)
    upper = max(v["continuous"]["maximum_ratio"] for v in rows.values())
    assert upper < required
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                clockfile,
                lab.flight.HERE / "DISPERSIVE_CLOCK.md",
                HERE / "OCCUPIED_WORK.md",
            )
        },
        galaxy="NGC2841",
        radii_kpc=radii,
        orbital_work_km2_s2=observations,
        required_outer_per_inner_lower=required,
        positive_mixture_upper=upper,
        results=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    output = run()
    path = HERE / "dispersive_work.json"
    if opt.write:
        path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    else:
        lab.flight.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), output
        )
    print(output["required_outer_per_inner_lower"], output["positive_mixture_upper"])
