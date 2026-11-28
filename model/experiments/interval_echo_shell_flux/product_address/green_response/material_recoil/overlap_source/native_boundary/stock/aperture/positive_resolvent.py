"""Integer residual enclosures of the four native source columns below the band."""

from fractions import Fraction

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu

DEN = 1024
BITS = 320
S = 1 << BITS


def integer_shift(twice, omega):
    z = omega * DEN
    assert z.denominator == 1
    out = (DEN // 2) * twice.copy()
    out.flat[:: len(out) + 1] -= z.numerator
    return out


def matvec(m, x):
    rows, cols = np.nonzero(m)
    y = np.zeros(x.shape, dtype=object)
    np.add.at(y, rows, m[rows, cols].astype(object)[:, None] * x[cols])
    return y


def floor_band(twice):
    # The floating edge selects a rational test frequency. The positive
    # supersolution certifies the spectral lower bound using only integers.
    e = np.linalg.eigvalsh(twice / 2)[0]
    candidate = Fraction(int(np.floor(e * DEN)), DEN)
    for _ in range(3):
        m = integer_shift(twice, candidate)
        y = splu(csc_matrix(m, dtype=float)).solve(np.ones(len(m)))
        v = np.array([round(x * 2**40) for x in y], dtype=object)[:, None]
        if v.min() > 0 and matvec(m, v).min() > 0:
            return candidate
        candidate -= Fraction(1, DEN)
    raise ArithmeticError("no positive certificate for the selected edge")


def columns(twice, omega, gap, sites):
    # The smallest eigenvalue of H-omega is at least gap. Use a floating
    # factorization for corrections and exact integer residuals for acceptance.
    m = integer_shift(twice, omega)
    factor = splu(csc_matrix(m, dtype=float))
    rhs = np.zeros((len(m), len(sites)), dtype=object)
    rhs[sites, np.arange(len(sites))] = S * DEN
    x = np.zeros(rhs.shape, dtype=object)
    for _ in range(12):
        residual = rhs - matvec(m, x)
        bound = Fraction(int(max(np.abs(residual).sum(0))), DEN) / gap
        ticks = -(-bound.numerator // bound.denominator)
        if ticks < 2**30:
            return np.maximum(x - ticks, 0), x + ticks, ticks
        dx = factor.solve(np.asarray(residual / S, dtype=float))
        x += np.array([[round(float(v) * S) for v in row] for row in dx], dtype=object)
    raise ArithmeticError("integer residual did not reach the requested accuracy")


def diagonal_product(first, middle, last, rows, cols, src):
    return ((first[:, rows] @ middle[np.ix_(src, cols)]) * last[:, cols]).sum(1)


def mixed_columns(twice, marks, omega, gap):
    off = twice.copy()
    np.fill_diagonal(off, 0)
    assert off.max() == 0
    assert connected_components(csc_matrix(off != 0), directed=False)[0] == 1
    sites = np.flatnonzero(marks.any(0))
    A, B = [np.flatnonzero(row[sites]) for row in marks]
    bounds = {}
    maximum = 0
    for a in (0, 32, 33):
        for b in (0, 32, 33):
            h = twice + np.diag(2 * (a * marks[0] + b * marks[1]))
            lo, hi, error = columns(h, omega, gap, sites)
            bounds[a, b] = (lo, hi)
            maximum = max(maximum, error)
    ends = []
    for index in (0, 1):

        def g(a, b):
            return bounds[a, b][index]

        mixed = np.zeros(len(twice), dtype=object)
        for a in (32, 33):
            for b in (32, 33):
                # R00 B R0b A Ra0 + R0b A Ra0 B Rab is the full
                # rectangular difference. All entries are nonnegative.
                left = diagonal_product(g(0, 0), g(0, b), g(a, 0), B, A, sites[B])
                right = diagonal_product(g(0, b), g(a, 0), g(a, b), A, B, sites[A])
                mixed += a * b * (left + right)
        ends.append(mixed)
    assert np.all(ends[0] > 0) and np.all(ends[1] >= ends[0])
    return ends, maximum


def shell_intervals(values, radii):
    # The conditional magnitude-two head puts weight 1/4 on each positive
    # removal frequency. The result is the mixed increment PER source D.
    denom = 16 * S**3
    out = []
    for d in range(max(radii) + 1):
        take = radii == d
        q = [Fraction(sum(x[take]), denom * int(take.sum())) for x in values]
        out.append(
            [
                float(np.nextafter(float(q[0]), -np.inf)),
                float(np.nextafter(float(q[1]), np.inf)),
            ]
        )
    return np.array(out)
