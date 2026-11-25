"""Certified range and Ramsey records for a single mediator between written cores."""

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import asin, atan2, cos, log, pi, sqrt
from pathlib import Path

import numpy as np
import sympy as sp
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from scipy.stats import binom

from model.candidate_laws.interval_echo.candidate import KAPPA
from model.executable_model import federation
from ..core_band_bound.measure import native_columns
from ..mediator import ground

HERE = Path(__file__).resolve().parent
ORDER = 24
SCALE = 10**9
SHOTS = 2**20
SEED = 311125320


def require(value, reason):
    if not value:
        raise ArithmeticError(reason)


def supports(edges):
    return sorted({p for edge in edges for p in edge})


def shape_survey(board, two_seams):
    shapes = set()
    for a in range(board.n):
        local = sorted(b for b in board.neighbours[a] if a // 12 == b // 12)
        groups = combinations(local, 2) if two_seams else ((b,) for b in local)
        for ends in groups:
            if two_seams and ends[1] not in board.neighbours[ends[0]]:
                continue
            shapes.add(tuple(sorted(tuple(sorted((a, b))) for b in ends)))
    first = min(shapes)

    def signature(edges):
        # Keep the centre and the outer port degrees distinct.
        flat = [p for e in edges for p in e]
        return sorted((flat.count(p), int(board.degree[p])) for p in set(flat))

    survey = shortest_path(
        board.adjacency(), directed=False, unweighted=True, indices=supports(first)
    ).min(axis=0)
    held = {}
    for edges in sorted(shapes):
        distance = int(survey[supports(edges)].min())
        if distance in (1, 2, 6, 9) and signature(edges) == signature(first):
            held.setdefault(distance, edges)
    require(len(held) == 4, "missing enrolled shape or distance")
    return first, held


def defect_twice(edges, strength):
    pts = supports(edges)
    index = {p: i for i, p in enumerate(pts)}
    w = np.zeros((len(pts), len(pts)), dtype=int)
    for a, b in edges:
        i, j = index[a], index[b]
        w[i, i] += 2 * strength + 1
        w[j, j] += 2 * strength + 1
        w[i, j] -= 1
        w[j, i] -= 1
    return w


def local_powers(lap, pts):
    # J = 6I - L has nonnegative integer entries and row sum six.
    n = lap.shape[0]
    walk = 6 * sparse.eye(n, dtype=np.int64, format="csr") - lap
    require(walk.min() >= 0 and (walk.sum(axis=1) == 6).all(), "invalid walk")
    require(6**ORDER <= np.iinfo(np.int64).max, "integer overflow bound")
    columns = np.zeros((n, len(pts)), dtype=np.int64)
    columns[pts, np.arange(len(pts))] = 1
    moments = [columns[pts].astype(object)]
    for _ in range(ORDER):
        columns = walk @ columns
        moments.append(columns[pts].astype(object))
    return moments


def certify_side(t, w2, moments, upper):
    """Test a strict rational sign of W - W R(t) W on the supported ports."""
    arg = 2 * t + 6
    a, b = arg.numerator, arg.denominator
    m = len(moments) - 1
    denominator = a ** (m + 1)
    truncated = sum(2 * b * b**k * a ** (m - k) * v for k, v in enumerate(moments))
    tail = (Fraction(6, 1) / arg) ** (m + 1) / t
    identity = np.eye(len(w2), dtype=int).astype(object)
    # For the upper binding endpoint require a strictly positive matrix even
    # at R_m + eps I. For the lower endpoint find a negative minor at R_m - eps I.
    offset = (1 if upper else -1) * tail.numerator * denominator
    bound = truncated * tail.denominator + offset * identity
    matrix = sp.Matrix(2 * denominator * tail.denominator * w2 - w2 @ bound @ w2)
    minors = [int(sp.sign(matrix[:i, :i].det())) for i in range(1, len(w2) + 1)]
    require(
        all(x > 0 for x in minors) if upper else -1 in minors,
        f"inconclusive ground-energy certificate at {t}",
    )
    return minors


def binding_interval(lap, edges, strength):
    # The sparse eigensolver proposes an interval. Exact resolvent signs
    # certify it includes the bottom of the full material operator.
    energy, residual = ground(lap / 2, edges, strength)
    require(residual < 1e-11, "compiled stationary source residual")
    center = round(-energy * SCALE)
    endpoints = [center - 1, center + 1]
    powers = local_powers(lap, supports(edges))
    w2 = defect_twice(edges, strength).astype(object)
    signs = [
        certify_side(Fraction(x, SCALE), w2, powers, upper)
        for x, upper in zip(endpoints, (False, True), strict=True)
    ]
    return {"binding_nanounits": endpoints, "leading_minor_signs": signs}


def envelope(d, first_degree, second_degree, amount):
    k = Fraction(KAPPA * amount**2)
    gap = k - Fraction(5, 2)
    nu = (k + 1) * sqrt(first_degree * second_degree)
    return float(nu * (1 + 6 / gap) * (3 / (gap + 3)) ** d)


def max_touch(edges):
    ports = [p for e in edges for p in e]
    return max(ports.count(p) for p in set(ports))


def meter(gain, wait, phase, rng):
    """All counts compatible with one binomial quantile and certified phase."""
    x = wait * (gain[0] + gain[1]) / 2
    probability = (1 + cos(x - phase)) / 2
    # |dP/dphase| <= 1/2 globally. Both ideal arm ground energies
    # lie in the certified intervals. A small floating allowance is included.
    # The two compiled arms each have energy residual below 1e-11.
    error = wait * ((gain[1] - gain[0]) / 4 + 2e-11) + 2e-14
    lo, hi = max(0, probability - error), min(1, probability + error)
    quantile = rng.random()
    low, high = (int(binom.ppf(quantile, SHOTS, p)) for p in (lo, hi))
    require(low == high, "certified clock interval admits two sample counts")
    return {"plus": low, "minus": SHOTS - low}


def calibration(gain, rng):
    banks, predictions = [], []
    for wait in (1 / 8, 8):
        banks.append([meter(gain, wait, phase, rng) for phase in (0, pi / 2)])
        real, imag = (2 * r["plus"] / SHOTS - 1 for r in banks[-1])
        angle = atan2(imag, real)
        if predictions:
            angle += 2 * pi * round((wait * predictions[0] - angle) / (2 * pi))
        predictions.append(angle / wait)
    # One 99% family over eight calibration quadratures across both shapes.
    phase_error = asin(sqrt(2) * sqrt(2 * log(16 / 0.01) / SHOTS))
    require(phase_error * (8 / (1 / 8) + 1) < pi, "ambiguous phase unwrap")
    return {
        "short_long_quadratures": banks,
        "gain": predictions[-1],
        "allowance": phase_error / 8,
    }


def run():
    native_columns()
    board = federation(2)
    lap = board.laplacian().astype(np.int64)
    rng = np.random.default_rng(SEED)
    control_rng = np.random.default_rng(SEED + 1)
    shapes = {}
    for name in ("seam", "two_seams"):
        first, held = shape_survey(board, name == "two_seams")
        first_binding = binding_interval(lap, first, KAPPA)
        rows = {}
        for d, second in sorted(held.items()):
            other = binding_interval(lap, second, KAPPA)
            combined = binding_interval(lap, first + second, KAPPA)
            a, b, ab = (
                u["binding_nanounits"] for u in (first_binding, other, combined)
            )
            gain = [
                max(0, ab[0] - max(a[1], b[1])) / SCALE,
                (ab[1] - max(a[0], b[0])) / SCALE,
            ]
            limit = envelope(d, max_touch(first), max_touch(second), 1)
            require(0 <= gain[0] <= gain[1] <= limit, "gain escaped envelope")
            rows[d] = {
                "second_edges": [list(e) for e in second],
                "second_binding": other,
                "joint_binding": combined,
                "gain_interval": gain,
                "upper_envelope": limit,
            }
        cal = calibration(rows[2]["gain_interval"], rng)
        ref = cal["gain"]
        require(ref > cal["allowance"], "three-link calibration is inconclusive")
        for d in (1, 6, 9):
            row = rows[d]
            wait = 1.0 if d == 1 else pi * (d + 1) / (3 * ref)
            row["wait"] = wait
            # pi recombination reads the dark outcome of the clock.
            row["dark_bright"] = meter(row["gain_interval"], wait, pi, rng)
            if d != 1:
                row["coulomb_dark_floor"] = cos(pi * cal["allowance"] / (2 * ref)) ** 2
        # In the changed law each nonwritten edge contributes only its graph
        # Laplacian. A constant vector has exactly zero energy on every panel.
        cut_energies = {}
        for edges in [first, *held.values(), *(first + e for e in held.values())]:
            pts = supports(edges)
            residual = lap @ np.ones(board.n, dtype=int)
            residual[pts] -= defect_twice(edges, 0) @ np.ones(len(pts), dtype=int)
            require((residual == 0).all(), "changed-law constant mode was lost")
            cut_energies[edges] = float(residual.sum()) / (2 * board.n)
        for d in (1, 6, 9):
            second = held[d]
            changed_gain = (
                min(cut_energies[first], cut_energies[second])
                - cut_energies[first + second]
            )
            row = rows[d]
            row["zero_contrast_dark_bright"] = meter(
                [changed_gain] * 2, row["wait"], pi, control_rng
            )
            row["receiver_pi_zero_contrast"] = meter(
                [changed_gain] * 2, row["wait"], 0, control_rng
            )
        shapes[name] = {
            "first_edges": [list(e) for e in first],
            "first_binding": first_binding,
            "calibration": cal,
            "distance_panels": {str(k): v for k, v in rows.items()},
        }
    return {
        "cells": int(board.cells),
        "shots_per_bank": SHOTS,
        "seed": SEED,
        "integer_moment_order": ORDER,
        "shapes": shapes,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as out:
            out.write(result)
    else:
        require(
            json.loads(result) == json.loads((HERE / "results.json").read_text()),
            "receipt differs",
        )
        print(
            "Certified the mediator range and reproduced the 320-cell receiver banks."
        )


if __name__ == "__main__":
    main()
