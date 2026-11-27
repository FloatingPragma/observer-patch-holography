"""Test a shared distribution of elapsed times against paired material clock records."""

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull
from scipy.stats import hypergeom

from model.experiments.count_intervals import clopper_pearson
from .measure import SHOTS, TIMER_SHOTS, TIMES

HERE = Path(__file__).resolve().parent
SETTINGS = ["blank"] + [
    f"{mode}{d}" for mode in ("free", "park") for d in (0, 1, 2, 3, 14)
]


def hands(times):
    return np.column_stack((np.cos(32 * times), np.cos(33 * times)))


@lru_cache(maxsize=4)
def time_polygon(limit):
    # An inner polygon suffices: each vertex is an actual admissible phase.
    times = np.linspace(0, limit, 4097)
    points = hands(times)
    vertices = ConvexHull(points).vertices
    return times[vertices], points[vertices]


def cut(polygon, axis, bound, sign):
    inside = []
    for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
        a, b = sign * (start[axis] - bound), sign * (end[axis] - bound)
        if a >= 0:
            inside.append(start)
        if a * b < 0:
            inside.append(start + a / (a - b) * (end - start))
    return np.array(inside)


def witness(low, high, ceiling):
    times, vertices = time_polygon(min(np.pi, ceiling))
    region = vertices
    for axis in (0, 1):
        for bound, sign in ((low[axis], 1), (high[axis], -1)):
            region = cut(region, axis, bound, sign)
            if not len(region):
                return None
    target = np.mean(region, axis=0)
    # A triangle fan represents any point of the convex polygon by three times.
    for j in range(1, len(vertices) - 1):
        triangle = np.array([0, j, j + 1])
        v, w = vertices[triangle[1:]] - vertices[0]
        a, b = np.linalg.solve(np.column_stack((v, w)), target - vertices[0])
        weights = np.array([1 - a - b, a, b])
        if weights.min() < -1e-10:
            continue
        weights = np.maximum(weights, 0)
        weights /= weights.sum()
        phases = times[triangle]
        recovered = weights @ hands(phases)
        if not np.all((low - 1e-10 <= recovered) & (recovered <= high + 1e-10)):
            raise ArithmeticError(
                "The constructed phase mixture misses a clock interval"
            )
        return {
            "times": phases.tolist(),
            "weights": weights.tolist(),
            "flip_probabilities": ((1 - recovered) / 2).tolist(),
        }
    raise ArithmeticError("Could not triangulate feasible intersection")


def compare(archive):
    appointments = len(TIMES)
    trials = appointments * len(SETTINGS)
    # Four binomial readings per panel and appointment, with two tails each.
    tail = 0.01 / (8 * trials)
    summary, source_tests = {}, {}
    example = None
    for setting in SETTINGS:
        bands, reference, arrivals = [], [], []
        for kind in ("even", "odd"):
            key = f"k8_{kind}_{setting}"
            counts = archive[key]
            if (counts.shape != (appointments, 4) or not np.issubdtype(counts.dtype, np.integer)
                    or np.any((counts < 0) | (counts > SHOTS)) or np.any(counts.sum(axis=1) != SHOTS)):
                raise ValueError(f"Incomplete coin and endpoint record: {key}")
            timer_counts = archive[key + "_timer"]
            if timer_counts.shape != (appointments,) or not np.issubdtype(timer_counts.dtype, np.integer):
                raise ValueError(f"Incomplete timing record: {key}")
            bands.append(clopper_pearson(counts[:, 2:].sum(axis=1), SHOTS, tail))
            timer = clopper_pearson(timer_counts, TIMER_SHOTS, tail)
            reference.append(16 * np.arcsin(np.sqrt(timer)))
            arrivals.append(counts[:, [1, 3]].sum(axis=1))
        low = np.column_stack([1 - 2 * upper for _, upper in bands])
        high = np.column_stack([1 - 2 * lower for lower, _ in bands])
        lower_timer = np.maximum(reference[0][0], reference[1][0])
        upper_timer = np.minimum(reference[0][1], reference[1][1])
        if np.any(lower_timer > upper_timer):
            raise ValueError("The timing references disagree")
        found = 0
        for i, (lo, hi, limit) in enumerate(zip(low, high, upper_timer, strict=True)):
            sample = witness(lo, hi, limit)
            if sample is not None:
                found += 1
            if setting == "free1" and i == 108:
                example = sample
        summary[setting] = found
        a, b = arrivals
        # Equality of the endpoint probability under a change of clock species.
        # Conditional on a+b, equal shot counts make both Fisher tails symmetric.
        p = np.minimum(1, 2 * hypergeom.cdf(np.minimum(a, b), 2 * SHOTS, a + b, SHOTS))
        i = int(np.argmin(p))
        source_tests[setting] = {
            "appointment_index": i,
            "even_at_clock": int(a[i]),
            "odd_at_clock": int(b[i]),
            "bonferroni_p": min(1.0, float(p[i]) * trials),
        }
    return {
        "appointments_per_setting": appointments,
        "settings": len(SETTINGS),
        "joint_count_band_coverage": 0.99,
        "count_bands_with_witness": summary,
        "example_free1_index108": example,
        "source_species_tests": source_tests,
    }


if __name__ == "__main__":
    with np.load(HERE / "counts.npz") as counts:
        print(json.dumps(compare(counts), indent=2))
