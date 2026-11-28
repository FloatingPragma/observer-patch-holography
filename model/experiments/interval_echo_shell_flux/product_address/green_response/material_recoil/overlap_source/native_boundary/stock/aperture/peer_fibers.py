"""Enumerate finite native return populations at fixed ordinary observer totals."""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import math

import numpy as np
from sympy import Matrix, ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from . import collective_source as source


@lru_cache(maxsize=1)
def apparatus():
    energy = source.tidal.lab.energy
    _, w, c, _, codes, h, ids, groups = energy.pair.setup()
    ordinary = np.vstack((w[:, :4].T, h.diagonal()[codes], energy.returned(w)[-1])).T
    averages = [
        tuple(Q(int(sum(row)), int(sum(ids == i))) for row in ordinary[ids == i].T)
        for i in range(groups)
    ]
    keys = [
        (int(cell), averages[g], int(head == 0))
        for cell, g, head in zip(c, ids, w[:, -1], strict=True)
    ]
    partial = np.isin(w[:, -1], source.PARTIAL)
    counts = Counter(
        key for key, exclude in zip(keys, partial, strict=True) if not exclude
    )
    types = sorted(counts)
    lookup = {key: i for i, key in enumerate(types)}
    cat = np.array(
        [-1 if exc else lookup[key] for key, exc in zip(keys, partial, strict=True)]
    )
    assert all(len(set(cat[ids == i])) == 1 for i in range(groups))
    n = np.array([counts[key] for key in types])
    protected = np.bincount(c[partial], minlength=4)
    unit = math.gcd(*n, *protected)
    # Include both head marginals. Each observable in the ordinary chart is
    # a component average measured on exchange-equilibrated native copies.
    rows = Matrix([[1, i // 2, i % 2, *features] for i, features, _ in types]).T
    den = math.lcm(*(Q(x).denominator for x in rows))
    A = DomainMatrix.from_Matrix(rows * den).convert_to(ZZ)
    D, U, V = smith_normal_decomp(A)
    assert U * A * V == D and abs(V.det()) == 1
    basis = V.to_Matrix()[:, A.rank() :].T.lll()
    assert rows * basis.T == Matrix.zeros(rows.rows, basis.rows)
    B = np.array(basis.tolist(), np.int64)
    return w, c, ids, cat, types, n, protected, unit, B


def projections(B, reference):
    # Project x=reference+B.T*t>=0 successively eliminating the shortest
    # integer basis directions. Keep only the strongest bound per parallel
    # normal. Each projected inequality is a nonnegative combination.
    rows = [
        (tuple(map(int, row)), int(base))
        for row, base in zip(B.T, reference, strict=True)
    ]
    stages = []
    for k in range(B.shape[0]):
        stages.append(rows)
        if k == B.shape[0] - 1:
            break
        plus = [r for r in rows if r[0][0] > 0]
        minus = [r for r in rows if r[0][0] < 0]
        reduced = {}

        def add(row, offset):
            g = math.gcd(*row)
            if not g:
                assert offset >= 0
                return
            normal = tuple(v // g for v in row)
            # Variables are integers, so floor division of the offset is exact.
            value = offset // g
            reduced[normal] = min(value, reduced.get(normal, value))

        for row, base in rows:
            if row[0] == 0:
                add(row[1:], base)
        for (p, pc), (m, mc) in product(plus, minus):
            add(
                tuple(p[0] * v - m[0] * u for u, v in zip(p[1:], m[1:], strict=True)),
                p[0] * mc - m[0] * pc,
            )
        rows = list(reduced.items())
    return stages


def populations(B, reference, skip=0):
    stages = projections(B, reference)

    def walk(i, tail):
        lo, hi = None, None
        for row, base in stages[i]:
            rem = base + sum(a * b for a, b in zip(row[1:], tail, strict=True))
            if row[0] > 0:
                x = -(rem // row[0])
                lo = x if lo is None else max(lo, x)
            elif row[0] < 0:
                x = rem // -row[0]
                hi = x if hi is None else min(hi, x)
            elif rem < 0:
                return
        assert lo is not None and hi is not None
        for x in range(lo, hi + 1):
            t = (x, *tail)
            if i > skip:
                yield from walk(i - 1, t)
            else:
                yield reference + np.array(t) @ B[skip:]

    yield from walk(len(stages) - 1, ())


BLANK_PAIRS = ((0, 1), (2, 3), (11, 12))


def batches(B, base, weights, types, collapse):
    """Yield native sequence counts and first moments on each projected fiber."""
    n = int(base.sum())
    fac = [math.factorial(i) for i in range(n + 1)]
    powers = [[int(d) ** v for v in range(n + 1)] for d in weights]
    if not collapse:
        for v in populations(B, base):
            assert min(v) >= 0
            number = fac[-1] // math.prod(fac[x] for x in v)
            number *= math.prod(p[x] for p, x in zip(powers, v, strict=True))
            yield v, number, number * v.astype(object), 1
        return
    # Two unimodular basis directions redistribute a fixed number of
    # blanks among three families without changing their totals or cells.
    pairs = BLANK_PAIRS
    rest = [i for i in range(len(base)) if all(i not in pair for pair in pairs)]
    for i, j in pairs:
        assert types[i][0] == types[j][0] and (types[i][-1], types[j][-1]) == (0, 1)
        assert weights[i] == 2 * weights[j] and np.array_equal(B[:2, i], -B[:2, j])
    slopes = B[:2, [j for _, j in pairs]]
    assert (
        np.all(slopes.sum(1) == 0)
        and abs(int(Matrix(slopes[:, :2].tolist()).det())) == 1
    )
    assert not B[:2, rest].any()
    for v in populations(B, base, skip=2):
        totals = [int(v[i] + v[j]) for i, j in pairs]
        blanks = sum(int(v[j]) for _, j in pairs)
        total = sum(totals)
        if min(totals) < 0 or min(v[rest]) < 0 or not 0 <= blanks <= total:
            continue
        # The ratio of the native multiplicities of blank and closed is
        # identical in each family. Vandermonde sums the binomial factors.
        number = fac[-1] // math.prod(fac[x] for x in (*totals, *v[rest]))
        labels = math.prod(
            powers[i][x] for (i, _), x in zip(pairs, totals, strict=True)
        )
        number *= (
            (labels >> blanks)
            * math.comb(total, blanks)
            * math.prod(powers[i][v[i]] for i in rest)
        )
        sums = number * v.astype(object)
        for (i, j), x in zip(pairs, totals, strict=True):
            b, rem = divmod(number * blanks * x, max(total, 1))
            assert rem == 0
            sums[i], sums[j] = number * x - b, b
        # Number of bounded triples of blank counts (inclusion-exclusion).
        count = 0
        for mask in range(8):
            rem = blanks - sum(x + 1 for i, x in enumerate(totals) if mask >> i & 1)
            if rem >= 0:
                count += (-1) ** mask.bit_count() * math.comb(rem + 2, 2)
        assert count > 0
        yield v, number, sums, count
