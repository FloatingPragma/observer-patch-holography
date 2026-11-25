"""Piecewise Euclidean three-spheres for the same six calibrated radar stations."""

from itertools import combinations
from pathlib import Path
import json
from math import ceil, floor

from mpmath import iv
import numpy as np

from .. import curved

HERE = Path(__file__).resolve().parent
ALL = tuple(range(6))
ANCHORS = tuple(range(4))
REPORT = HERE / "results.json"


def raw_bounds():
    iv.dps = 35
    original = json.loads((HERE.parent / "curved.json").read_text())

    def interval(s):
        return iv.mpf(s.strip("[]").split(","))

    return np.array(
        [[interval(x) for x in row] for row in original["directed_mean_range_bounds"]],
        dtype=object,
    )


def observed_squared(means):
    """Freeze the four original gains. Each outside clock enrolls against the four."""
    logs = np.array(
        [[iv.ln(means[i, j] / means[j, i]) if i != j else 0 for j in ALL] for i in ALL],
        dtype=object,
    )
    gains = [sum(logs[i, j] for j in ANCHORS) / 4 for i in ANCHORS]
    gains += [sum(logs[i, j] + gains[j] for j in ANCHORS) / 4 for i in (4, 5)]
    return np.array(
        [
            [
                means[i, j] * means[j, i] * iv.exp(-gains[i] - gains[j])
                if i != j
                else iv.mpf(0)
                for j in ALL
            ]
            for i in ALL
        ],
        dtype=object,
    )


def dot(square, center, j, k):
    return (square[center, j] + square[center, k] - square[j, k]) / 2


def face_area(square, face):
    a, b, c = face
    return iv.sqrt(square[a, b] * square[a, c] - dot(square, a, b, c) ** 2) / 2


def solid_angle(square, center, others):
    lengths = [iv.sqrt(square[center, j]) for j in others]
    gram = np.array(
        [[dot(square, center, j, k) for j in others] for k in others], dtype=object
    )
    determinant = curved.det_interval(gram)
    if determinant.a <= 0:
        raise ValueError("degenerate incident tetrahedron")
    denominator = lengths[0] * lengths[1] * lengths[2] + sum(
        lengths[i] * gram[(i + 1) % 3, (i + 2) % 3] for i in range(3)
    )
    if denominator.a <= 0:
        raise ValueError("trihedral angle outside acute denominator domain")
    return 2 * iv.atan2(iv.sqrt(determinant), denominator), iv.sqrt(determinant) / 6


def interval_numbers(value):
    # Round decimal millionths outwards; padding also dominates binary conversion.
    return [
        (floor(float(value.a) * 1e6) - 1) / 1e6,
        (ceil(float(value.b) * 1e6) + 1) / 1e6,
    ]


def fillings(square):
    """Bound each metric for a pair of missing triangles containing the original survey tetrahedron."""
    tetra = {}
    for sites in combinations(ALL, 4):
        angles, altitudes = {}, []
        for i in sites:
            rest = tuple(j for j in sites if j != i)
            angles[i], volume = solid_angle(square, i, rest)
            altitudes.append(3 * volume / face_area(square, rest))
        tetra[sites] = (angles, volume, min(float(h.a) for h in altitudes))
    rows = []
    for inside in combinations(ALL, 3):
        outside = tuple(i for i in ALL if i not in inside)
        if 0 not in inside or ((4 in inside) == (5 in inside)):
            continue
        chambers = [sites for sites in tetra if len(set(sites) & set(inside)) == 2]
        assert len(chambers) == 9 and ANCHORS in chambers
        omega = [sum(tetra[t][0][i] for t in chambers if i in t) for i in ALL]
        assert min(tetra[t][2] for t in chambers) > 1
        rows.append(
            dict(
                triangles=[list(inside), list(outside)],
                solid_angles=[interval_numbers(o) for o in omega],
                total_volume=interval_numbers(sum(tetra[t][1] for t in chambers)),
                angle_ratio_0_to_1=interval_numbers(omega[0] / omega[1]),
            )
        )
    return rows, interval_numbers(tetra[ANCHORS][1])


def certify(means):
    iv.dps = 35
    sq = observed_squared(means)
    # On the centered subspace of dimension five, the squared equilateral
    # length matrix has Gram 9I/2. Bound the symmetric perturbation by its
    # maximum absolute row sum (both centering projections have norm one).
    budget = max(sum(abs(sq[i, j] - 9).b for j in ALL if i != j) for i in ALL)
    gram_lower = (9 - budget) / 2
    if gram_lower.a <= 0:
        raise ValueError("six distances lack a certified Euclidean simplex realization")
    rows, volume = fillings(sq)
    return dict(
        centered_gram_lower=interval_numbers(gram_lower)[0],
        original_four_volume=volume,
        fillings=rows,
    )


if __name__ == "__main__":
    print(json.dumps(certify(raw_bounds()), indent=2))
