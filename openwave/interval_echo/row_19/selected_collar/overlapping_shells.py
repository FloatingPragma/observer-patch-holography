"""Check whether the source reads at every center have a common positive rest-ball density."""

import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np

from model.experiments.interval_echo_einstein_station.collar_scaling.station_clock.observer_survey import (
    measure,
)
from . import resolved_shells as shells

HERE = shells.HERE


def profile(adj, home, reach):
    hist = [1]
    seen = {home}
    ring = {home}
    for _ in range(reach):
        ring = {v for u in ring for v in adj[u]} - seen
        seen.update(ring)
        hist.append(len(ring))
    return tuple(hist)


def kernel(hist, R, offset, lower):
    # Offset in sixteenths of a measured station hop. A nearest station
    # bounds all distances even for mass off the measured station graph.
    sq = (16 * R) ** 2
    return sum(
        n * max(0, sq - (16 * d + offset if lower else max(0, 16 * d - offset)) ** 2)
        for d, n in enumerate(hist)
    )


def acquire_centers(hist, total_sites, radius, word_chance, rng):
    # Draw the center class first, then two ordered distinct addresses.
    # The acquisition never uses the predicted average pair probability.
    items = list(hist.items())
    visits = rng.multinomial(
        shells.STARTS, np.array([c for _, c in items]) / sum(hist.values())
    )
    hits = 0
    for (h, _), shots in zip(items, visits, strict=True):
        size = sum(h[: radius + 1])
        first = rng.binomial(shots, size / total_sites)
        hits += int(rng.binomial(first, (size - 1) / (total_sites - 1)))
    positive = int(rng.binomial(hits, float(word_chance)))
    tape = [positive, hits - positive, shells.STARTS - hits]
    f = positive / shells.STARTS
    err = shells.reg.src.probability_error(np.asarray(f), shells.STARTS)
    return tape, shells.reg.band(f - err, f + err)


def run():
    shells.iv.dps = 35
    prior = shells.dist.read("resolved_shells")
    rng = np.random.default_rng(28115732)
    chance = prior["word_event_probabilities"]
    kappa = shells.iv.mpf(prior["conditional_entropy_slope"])
    rows = []
    for level, R in ((1, 1), (3, 2), (5, 4), (7, 8)):
        adj = measure.atlas(level)
        N = len(adj)
        # The twenty congruent original faces are related by rotations
        # preserving all station incidences. Keep every center in one face.
        hist = Counter(profile(adj, i, 3 * R + 1) for i in range(4**level))
        assert 20 * sum(hist.values()) == N
        if level <= 3:
            assert Counter(profile(adj, i, 3 * R + 1) for i in range(N)) == Counter(
                {h: 20 * c for h, c in hist.items()}
            )
        output = []
        for radius in (R, 2 * R):
            chance2 = F(
                sum(
                    c * sum(h[: radius + 1]) * (sum(h[: radius + 1]) - 1)
                    for h, c in hist.items()
                ),
                sum(hist.values()) * N * (N - 1),
            )
            fresh, lo = acquire_centers(hist, N, radius, F(chance["fresh"]), rng)
            shared, hi = acquire_centers(hist, N, radius, F(chance["shared"]), rng)
            source = (
                (hi - lo)
                * kappa
                / shells.interval(F(chance["shared"]) - F(chance["fresh"]))
            )
            output.append(
                {
                    "radius": radius,
                    "pair_probability": str(chance2),
                    "word_counts": {"fresh": fresh, "shared": shared},
                    "measured_entropy": shells.bounds(source),
                }
            )
        ratio = 2 * F(output[1]["pair_probability"]) / F(output[0]["pair_probability"])
        measured = (
            2
            * shells.iv.mpf(output[1]["measured_entropy"])
            / shells.iv.mpf(output[0]["measured_entropy"])
        )
        bands = {}
        for offset in range(16 * R):
            upper = max(
                F(kernel(h, 2 * R, offset, False), kernel(h, R, offset, True))
                for h in hist
            )
            bands[str(offset)] = {
                "max_larger_over_smaller": str(upper),
                "excluded_by_counts": bool(measured.a > shells.interval(upper).b),
            }
        rows.append(
            {
                "level": level,
                "cells": N,
                "inner_hops": R,
                "center_profile_multiplicities": [
                    [list(h), 20 * c] for h, c in sorted(hist.items())
                ],
                "observations": output,
                "source_moment_ratio": str(ratio),
                "measured_moment_ratio": shells.bounds(measured),
                "offsets_in_sixteenths": bands,
            }
        )
    return {
        "starts": shells.STARTS,
        "levels": rows,
        "seed": 28115732,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in ("resolved_shells.json", "OVERLAPPING_SHELLS.md")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = run()
    path = HERE / "overlapping_shells.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    for z in data["levels"]:
        print(
            z["level"],
            z["measured_moment_ratio"],
            float(F(z["offsets_in_sixteenths"]["8"]["max_larger_over_smaller"])),
            max(
                [
                    int(i)
                    for i, b in z["offsets_in_sixteenths"].items()
                    if b["excluded_by_counts"]
                ],
                default=-1,
            ),
        )
