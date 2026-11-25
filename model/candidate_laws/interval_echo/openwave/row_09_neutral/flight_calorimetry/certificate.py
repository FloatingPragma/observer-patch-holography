"""Exact zero-phase root brackets and charge-to-pair thresholds on this finite ring."""

from collections import Counter
from itertools import combinations
import numpy as np
import sympy as sp

from .apparatus import RING, sector


def exact_support():
    words, h, costs, shift = sector("pair", 8)
    # The fourteen translations put each ordered pair in one of 13 orbits.
    # Column sums from one representative give the quotient in integer 2H units.
    relative = np.array(
        [(RING.index(b) - RING.index(a)) % len(RING) - 1 for a, b in words]
    )
    twice = np.rint(2 * h).astype(int)
    np.testing.assert_array_equal(2 * h, twice)
    reduced = np.zeros((13, 13), int)
    for d in range(13):
        columns = np.flatnonzero(relative == d)
        for e in range(13):
            sums = twice[np.ix_(relative == e, columns)].sum(axis=0)
            assert np.all(sums == sums[0])
            reduced[e, d] = sums[0]
    q = sp.Matrix(reduced).charpoly().as_poly()
    brackets = [(49_066_660, 49_066_670), (50_062_495, 50_062_505)]
    counts = [
        int(q.count_roots(sp.Rational(2 * a, 10**6), sp.Rational(2 * b, 10**6)))
        for a, b in brackets
    ]
    assert counts == [1, 1]
    assert q.count_roots(2 * 36, sp.oo) == 2
    assert set(costs.tolist()) == {32, 48}

    # An arbitrary charged final state containing also an opposite-load pair
    # has two identical signs and one opposite sign. Count *every* placement.
    size = len(RING)
    three = Counter()
    for a, b in combinations(range(size), 2):
        for c in range(size):
            if c in (a, b):
                continue
            word = np.ones(size, dtype=int)
            word[[a, b]] = 0
            word[c] = 2
            three[int(8 * np.sum((word - np.roll(word, 1)) ** 2))] += 1
    assert min(three) == 32
    for kind in ("vacancy", "excess"):
        _, single, diagonal, _ = sector(kind, 8)
        upper = np.max(
            np.diag(single) + abs(single - np.diag(np.diag(single))).sum(axis=1)
        )
        assert upper == 18 and set(diagonal.tolist()) == {16}

    return dict(
        zero_phase_quotient_twice=reduced.tolist(),
        energy_brackets_micro_units=brackets,
        roots_per_bracket=counts,
        other_zero_phase_roots_above_36=0,
        all_pair_contrast_histogram=dict(Counter(costs.tolist())),
        charged_single_full_spectrum_interval=[16, 18],
        charged_daughter_plus_pair_contrast_histogram=dict(sorted(three.items())),
    )
