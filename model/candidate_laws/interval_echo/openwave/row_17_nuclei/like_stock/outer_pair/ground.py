"""Enclose true ground shell weights with deposited positive trial vectors."""

import json
from fractions import Fraction

import numpy as np

from .measure import Assembly, HERE, distances


def negative_inertia(matrix):
    """Signs of successive nonzero leading minors via fraction-free elimination."""
    a = [r.copy() for r in matrix]
    previous = 1
    count = 0
    for j in range(len(a)):
        pivot = a[j][j]
        if pivot == 0:
            raise ValueError("zero leading minor")
        count += (pivot < 0) != (previous < 0)
        for i in range(j + 1, len(a)):
            for k in range(i, len(a)):
                value = pivot * a[i][k] - a[i][j] * a[j][k]
                q, r = divmod(value, previous)
                if r:
                    raise ArithmeticError("inexact Bareiss division")
                a[i][k] = a[k][i] = q
        previous = pivot
    return count


def trial(n, row, threshold):
    assembly = Assembly(n, 0, 8)
    h = 2 * assembly.h.toarray()
    assert np.array_equal(h, np.rint(h))
    twice = h.astype(int).tolist()
    inertia_matrix = [
        [
            v * threshold.denominator - (2 * threshold.numerator if i == j else 0)
            for j, v in enumerate(coefficients)
        ]
        for i, coefficients in enumerate(twice)
    ]
    if negative_inertia(inertia_matrix) != 1:
        raise ValueError("threshold does not isolate the lowest eigenvector")
    one = n == 1
    x = row["positive_one" if one else "positive_two"]
    assert len(x) == len(h) and min(x) > 0
    hx = [
        sum(v * w for v, w in zip(coefficients, x, strict=True))
        for coefficients in twice
    ]
    ratios = [Fraction(t, 2 * s) for t, s in zip(hx, x, strict=True)]
    lo, hi = min(ratios), max(ratios)
    assert list(map(str, [lo, hi])) == row["one" if one else "two"] and hi < threshold
    # For normalized x, ||(H-lo)x|| <= hi-lo componentwise in magnitude.
    # The orthogonal complement has energy >= threshold.
    eps = (hi - lo) / (threshold - lo)
    norm = sum(v * v for v in x)
    p = [Fraction(v * v, norm) for v in x]
    return assembly.w, p, eps


def certificate():
    row = json.loads((HERE.parent / "two_counts.json").read_text())["rows"][0]
    assert row["k"] == 8
    dist = distances()
    words, prob, eps = trial(2, row, Fraction("84.80"))
    joined = [
        sum(
            p
            for w, p in zip(words, prob, strict=True)
            if dist[tuple(np.flatnonzero(w == 1))] == d
        )
        for d in (1, 2, 3)
    ]
    words, prob, eps1 = trial(1, row, Fraction("42.46"))
    addr = [int(np.argmax(w)) for w in words]
    copies = [
        sum(
            p * q
            for i, p in zip(addr, prob, strict=True)
            for j, q in zip(addr, prob, strict=True)
            if dist[i, j] == d
        )
        for d in (0, 1, 2, 3)
    ]

    # Independent tensor products have trace distance at most 2*eps1.
    # Intersect with distinct-address support before taking the conditional ratio.
    def conditional(selected, total, error):
        return [
            float((selected - error) / (total + error)),
            float((selected + error) / (total - error)),
        ]

    far = sum(joined[1:])
    return dict(
        joined_far=[float(far - eps), float(far + eps)],
        independent_far=conditional(sum(copies[2:]), sum(copies[1:]), 2 * eps1),
        ground_distances=[float(eps), float(eps1)],
    )


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2))
