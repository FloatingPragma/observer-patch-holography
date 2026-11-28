"""Predict the spatial cut of sharp mobile doubles under native exchange."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, sqrt
from pathlib import Path

import numpy as np

from cosmology.screen_coherence import measure as screen
from cosmology.rar_deep_regime.lensing_profile import check
from . import measure

RECEIPT = measure.HERE / "stationary.json"
COPIES, STEPS, SEED = 4096, 11, 14031127


def equilibrium(n, k):
    """First two moments of boundary length for uniform k-subsets, on a simple cubic graph."""
    edges = 3 * n // 2
    shared = 3 * n
    p = Fraction(2 * k * (n - k), n * (n - 1))
    both = Fraction(
        4 * k * (k - 1) * (n - k) * (n - k - 1), n * (n - 1) * (n - 2) * (n - 3)
    )
    variance = (
        edges * p * (1 - p)
        + 2 * shared * (p / 2 - p**2)
        + 2 * (comb(edges, 2) - shared) * (both - p**2)
    )
    return {
        "mean": float(edges * p),
        "variance": float(variance),
        "exact_mean": str(edges * p),
        "exact_variance": str(variance),
    }


def exchange_moments(first, second, x, y):
    # Different matched pairs receive independent fair exchanges. Two ends of
    # the same matched pair share one exchange and have invariant joint product.
    within = second[x, y].copy()
    pair_mean = (first[x] + first[y]) / 2
    first[x] = first[y] = pair_mean
    second[x, :] = second[y, :] = (second[x, :] + second[y, :]) / 2
    second[:, x] = second[:, y] = (second[:, x] + second[:, y]) / 2
    second[x, y] = second[y, x] = within
    np.fill_diagonal(second, first)


def forecast(board):
    groups = [tuple(m[:2]) for m in screen.observations.schedules(board)[0]]
    near, _ = screen.near_half(board.cells, groups)
    first = np.isin(np.arange(board.cells), near).astype(float)
    second = np.outer(first, first)
    for t in range(STEPS):
        exchange_moments(first, second, *groups[t % len(groups)])
    np.testing.assert_allclose(
        second.sum(axis=1), len(near) * first, atol=1e-10, rtol=0
    )
    x, y = (np.asarray(board.external) // 12).T
    marginal_sum = (first[x] + first[y]).sum()
    exact = marginal_sum - 2 * second[x, y].sum()
    independent = marginal_sum - 2 * (first[x] * first[y]).sum()
    return groups, near, x, y, float(exact), float(independent)


def sample_cuts(board, level, groups, near, x, y):
    rng = np.random.default_rng(np.random.SeedSequence([SEED, level]))
    copies = np.zeros((COPIES, board.cells), dtype=bool)
    copies[:, near] = True
    for t in range(STEPS):
        a, b = groups[t % len(groups)]
        exchange = rng.integers(0, 2, (COPIES, len(a)), dtype=np.uint8).astype(bool)
        screen.paired_moves(copies, a, b, exchange)
    np.testing.assert_array_equal(copies.sum(axis=1), len(near))
    values = np.concatenate(
        [(c[:, x] != c[:, y]).sum(axis=1) for c in np.split(copies, 64)]
    )
    return (
        np.bincount(values).tolist(),
        float(values.mean()),
        float(values.std(ddof=1) / sqrt(COPIES)),
    )


def archive_source(path=None):
    published = json.loads((screen.HERE / "measurements.json").read_text())["archive"][
        "2"
    ]
    row = next(r for r in published["rows"] if r["service"] == published["services"])
    n, k = published["cells"], published["initial_doubles"]
    # The pre-existing covariance is centered at 0.16, not at k/n.
    p11 = (
        row["equal_neighbor_cov"]
        + 2 * screen.DENSITY * row["matter_mean"]
        - screen.DENSITY**2
    )
    cut = round((3 * n) * (k / n - p11))
    if path:
        if sha256(path.read_bytes()).hexdigest() != published["npz_sha256"]:
            raise ValueError("Unexpected mobile-double archive")
        with np.load(path, allow_pickle=False) as raw:
            x, y = (raw["external"] // 12).T
            bits = np.unpackbits(raw["double"][-1])
            assert bits.sum() == k and np.count_nonzero(bits[x] != bits[y]) == cut
    return {
        "cells": n,
        "doubles": k,
        "service": row["service"],
        "direct_cut": cut,
        "stationary": equilibrium(n, k),
        "archive_sha256": published["npz_sha256"],
        "summary_sha256": sha256(
            (screen.HERE / "measurements.json").read_bytes()
        ).hexdigest(),
    }


def run(path=None):
    instrument = screen.observations.instruments()[1]
    np.testing.assert_allclose(instrument[[0, 2], [2, 0]], 0.5, rtol=0, atol=1e-13)
    returned = measure.admission.material_panel(1, 1, 8)["records"][-1]["d"][-1]
    per_seam = float(measure.charge(measure.UNIFORM_PLACEMENT, returned))
    data = []
    for level in (1, 2, 3):
        board = (
            measure.load_refinement(measure.HERE / "FEDERATION_LEVEL3.npz", level=3)
            if level == 3
            else measure.federation(level)
        )
        groups, near, x, y, exact, independent = forecast(board)
        hist, mean_cut, se = sample_cuts(board, level, groups, near, x, y)
        stationary = equilibrium(board.cells, len(near))
        data.append(
            {
                "cells": board.cells,
                "doubles_each_trial": len(near),
                "stationary": stationary,
                "stationary_charge": stationary["mean"] * per_seam,
                "forecast_cut": exact,
                "independent_marginals_cut": independent,
                "forecast_charge": exact * per_seam,
                "cut_histogram": hist,
                "sample_mean_cut": mean_cut,
                "sample_se_cut": se,
                "sample_conditional_charge": mean_cut * per_seam,
            }
        )
    old = archive_source(path)
    old["conditional_charge"] = per_seam * old["direct_cut"]
    return {
        "protocol_sha256": sha256(
            (measure.HERE / "STATIONARY.md").read_bytes()
        ).hexdigest(),
        "samples_per_board": COPIES,
        "held_service": STEPS,
        "seed": SEED,
        "gain_per_mismatched_seam": per_seam,
        "panels": data,
        "published_archive": old,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archive", type=Path, help="optional original level-2 universe.npz"
    )
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run(args.archive)
    if args.write:
        if not args.archive:
            parser.error("The original archive is required to create the receipt")
        with RECEIPT.open("x") as out:
            json.dump(result, out, allow_nan=False, indent=2)
            out.write("\n")
    else:
        check(json.loads(RECEIPT.read_text()), result)
    for row in result["panels"]:
        print(
            row["cells"],
            row["forecast_cut"],
            row["sample_mean_cut"],
            row["sample_se_cut"],
            row["independent_marginals_cut"],
            row["stationary"]["mean"],
        )
    print(result["published_archive"])
