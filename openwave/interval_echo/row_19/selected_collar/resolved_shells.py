"""Read a distributed modular increment on the protected refining screen and test radial positivity."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json

import numpy as np

from model.experiments.interval_echo_einstein_station.collar_scaling.station_clock import (
    measure as clock,
    word_entropy as words,
)
from . import distributed_source as dist

HERE, iv, reg = dist.HERE, dist.iv, dist.reg
bounds = dist.lab.six.bounds
SEED = 28112129
STARTS = 2**48


def interval(q):
    q = F(q)
    return iv.mpf(q.numerator) / q.denominator


def event_counts(p, N, rng):
    # Two-stage sampling keeps every rejected preparation as an explicit miss.
    pairs = int(rng.binomial(STARTS, float(N)))
    yes = int(rng.binomial(pairs, float(p)))
    counts = [yes, pairs - yes, STARTS - pairs]
    freq = yes / STARTS
    error = reg.src.probability_error(np.asarray(freq), STARTS)
    return counts, reg.band(freq - error, freq + error)


def run():
    iv.dps = 35
    table = json.loads((clock.HERE / "results.json").read_bytes())
    word = json.loads((clock.HERE / "word_entropy.json").read_bytes())
    kappa = iv.mpf(word["reversed_reference_conditional_derivative"])
    groups = words.probabilities()
    # Both arm probabilities refer to the same independently read word event.
    fresh = sum(F(p * c, 2592) for p, d, c in groups if d > 0)
    shared = sum(F((p + d) * c, 2592) for p, d, c in groups if d > 0)
    gain = kappa / interval(shared - fresh)
    rng = np.random.default_rng(SEED)
    vertices, edges = clock.icosahedron_from_coordinates()
    edges = set(edges)
    faces = [f for f in combinations(range(12), 3) if set(combinations(f, 2)) <= edges]
    panels = {row["level"]: row for row in table["midpoint_surveys"]}
    rows = []
    for level in range(max(panels) + 1):
        if level in panels:
            survey = panels[level]
            R = survey["radius"]
            observed = np.bincount(clock.dual_distances(faces))
            H = len(observed) - 1
            N = len(faces)
            assert H == survey["eccentricity"] and N == survey["cells"]
            if level < 3:
                assert list(observed) == table["hardware"][level]["distance_histogram"]
            sizes = np.cumsum(observed)
            probs = [F(int(n) * (int(n) - 1), N * (N - 1)) for n in sizes[: R + 2]]
            radii = [F(i, H) for i in range(R + 2)]
            assert all(clock.corner_ball(i) == sizes[i] for i in range(R + 1))
            # Factor out the same positive conditional source reading κ.
            moment = [r * p for r, p in zip(radii, probs, strict=True)]
            slopes = [
                (b - a) / (y**2 - x**2)
                for a, b, x, y in zip(
                    moment[:-1], moment[1:], radii[:-1], radii[1:], strict=True
                )
            ]
            weights = [s - z for s, z in zip(slopes, [F(0), *slopes[:-1]], strict=True)]
            # Weights are in units κ/π. Only points through the outer enrolled
            # radius enter the fit. Use the next hop only to bound enclosed mass.
            assert all(
                sum(
                    w * max(r * r - s * s, 0)
                    for s, w in zip(radii[:-1], weights, strict=True)
                )
                == y
                for r, y in zip(radii, moment, strict=True)
            )
            dense, pd = event_counts(fresh, probs[R], rng)
            marked, pm = event_counts(shared, probs[R], rng)
            q = kappa * interval(probs[R])
            ell = radii[R]
            # One step below the outer hop, the same sources are admitted
            # at all three sub-hop radii. F=q*r is strictly concave in r².
            probe = [(F(R - 1) + F(j, 4)) / H for j in (1, 2, 3)]
            a, b, c = probe
            subhop = probs[R - 1] * (1 / (b + c) - 1 / (a + b))
            rows.append(
                {
                    "level": level,
                    "cells": N,
                    "eccentricity": H,
                    "radius_hops": R,
                    "shell_populations": list(map(int, observed[: R + 2])),
                    "eligible_pair_probability": str(probs[R]),
                    "forecast_entropy_tangent": bounds(q),
                    "direct_entropy_tangent": bounds((pm - pd) * gain),
                    "word_counts": {"fresh": dense, "shared": marked},
                    "shell_weights_in_kappa_over_pi": [str(w) for w in weights[:R]],
                    "next_weight_in_kappa_over_pi": str(weights[R]),
                    "enclosed_per_radius_cubed": [
                        bounds(interval(slopes[k] / ell**3) * kappa / iv.pi)
                        for k in (R - 1, R)
                    ],
                    "subhop_radii": list(map(str, probe)),
                    "subhop_secant_change": bounds(interval(subhop) * kappa),
                    "normalized_pair_coefficient": str(probs[R] / ell**4),
                }
            )
        if level < max(panels):
            vertices, faces = clock._subdivide(vertices, faces)
    # n(R) ~ (5/4)R², H/m -> 6, N=20m² on the corner midpoint continuation.
    coefficient = F(81, 16)
    return {
        "preparation": "protected reversed word tag; shared minus fresh endpoints",
        "seed": SEED,
        "starts": STARTS,
        "word_event_probabilities": {"fresh": str(fresh), "shared": str(shared)},
        "conditional_entropy_slope": bounds(kappa),
        "levels": rows,
        "limit_entropy_per_fourth": bounds(interval(coefficient) * kappa),
        "equivalent_constant_stress": bounds(
            15 * interval(coefficient) * kappa / (8 * iv.pi**2)
        ),
        "limit_enclosed_per_cubed": bounds(
            5 * interval(coefficient) * kappa / (2 * iv.pi)
        ),
        "inputs": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "RESOLVED_SHELLS.md",
                clock.HERE / "results.json",
                clock.HERE / "word_entropy.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = run()
    path = HERE / "resolved_shells.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    for row in data["levels"]:
        print(
            row["level"],
            min(map(F, row["shell_weights_in_kappa_over_pi"])),
            row["forecast_entropy_tangent"],
            row["direct_entropy_tangent"],
            row["enclosed_per_radius_cubed"],
        )
