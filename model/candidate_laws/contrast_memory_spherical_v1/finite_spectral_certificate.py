"""Integer certificates for the sharp-pair bound on the three enrolled federations."""

from fractions import Fraction

import numpy as np
from scipy.sparse.csgraph import connected_components

from . import theorem

ADJACENT_CEILING = Fraction(5701, 1000)
COUPLING_CEILING = Fraction(347, 20)  # 4 ||QKP||²


def positive_leading_minors(matrix):
    """Check Sylvester's criterion by fraction-free symmetric elimination."""
    if not np.array_equal(matrix, matrix.T):
        return False
    remaining = np.asarray(matrix, dtype=object)
    previous = 1
    while remaining.size:
        pivot = remaining[0, 0]
        if pivot <= 0:
            return False
        numerator = pivot * remaining[1:, 1:] - np.outer(remaining[1:, 0], remaining[0, 1:])
        if np.any(numerator % previous):
            raise ArithmeticError("Nonintegral Bareiss division")
        remaining = numerator // previous
        previous = pivot
    return True


def integer_twice(block):
    """The exchange blocks are integer half-counts of moves and incident seams."""
    scaled = 2 * block
    result = scaled.astype(np.int64)
    if not np.array_equal(scaled.data, result.data):
        raise ValueError("Exchange entries must be half-integers")
    return result


def certify(car):
    adjacent, coupling = theorem.blocks(car)
    twice_a = integer_twice(adjacent)
    twice_b = integer_twice(coupling)
    components, labels = connected_components(twice_a, directed=False)
    ordered = np.argsort(labels, kind="stable")
    pieces = np.split(ordered, np.flatnonzero(np.diff(labels[ordered])) + 1)
    seen = set()
    for rows in pieces:
        block = twice_a[rows][:, rows].toarray()
        key = block.tobytes()
        if key in seen:
            continue
        seen.add(key)
        # 1000(5.701 I - A) = 5701 I - 500(2A).
        if not positive_leading_minors(5701 * np.eye(len(rows), dtype=np.int64) - 500 * block):
            raise ArithmeticError("The adjacent block exceeds the proposed ceiling")

    gram = (twice_b.T @ twice_b).tocsr()
    if np.any(gram.data < 0):
        raise ArithmeticError("The comparison needs an entrywise nonnegative Gram")
    # Deterministic positive test vector. Any successful iterate is sufficient:
    # Mx < (347/20)x bounds the spectral radius by a weighted row-sum norm.
    x = np.full(gram.shape[0], 10000, dtype=np.int64)
    largest_row_sum = int(np.max(gram.sum(axis=1)))
    max_entry = np.iinfo(np.int64).max // max(347, 20 * largest_row_sum)
    for updates in range(5001):
        if int(x.max()) > max_entry:
            raise ArithmeticError("Integer witness exceeds int64 product range")
        mx = gram @ x
        product = 20 * mx
        slack = 347 * x - product
        if int(slack.min()) > 0:
            # A positive difference puts the Rayleigh quotient strictly above
            # the old 17.3 estimate. Use unbounded integers for these dot products.
            legacy_excess = sum(int(v) * (10 * int(w) - 173 * int(v))
                                for v, w in zip(x, mx, strict=True))
            return {"cells": car.cells, "seams": gram.shape[0], "adjacent_components": components,
                    "distinct_adjacent_blocks": len(seen), "witness_updates": updates,
                    "minimum_integer_slack": int(slack.min()), "old_rayleigh_excess": legacy_excess}
        x = 10000 + product // 347
    raise ArithmeticError("No positive integer witness for the proposed Gram ceiling")


def sharp_pair_bound(kappa):
    """Probability after free evolution, for the three boards checked by certify()."""
    gap = 2 * kappa - ADJACENT_CEILING
    return min(Fraction(1), COUPLING_CEILING / gap**2) if gap > 0 else Fraction(1)


if __name__ == "__main__":
    from model.executable_model import federation

    for level in range(3):
        print(certify(federation(level)))
