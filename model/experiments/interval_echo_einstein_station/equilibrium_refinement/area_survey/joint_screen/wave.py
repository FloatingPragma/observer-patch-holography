"""Bound the half-time native propagator in the enrolled 276-word sector."""

from fractions import Fraction
from math import factorial

import numpy as np

UNIT = 2**80
BASE = 2**21
DEGREE = 29
SQUARES = 7


def split(x):
    # Four small limbs including a signed high limb. Binary64 exactly represents
    # all of the integer multiply-adds in each BLAS dot product.
    parts = []
    for _ in range(3):
        parts.append(np.asarray(x % BASE, dtype=float, order="C"))
        x = x // BASE
    assert np.max(abs(x)) < BASE
    parts.append(np.asarray(x, dtype=float, order="C"))
    return parts


def product(x, y):
    n = len(x[0])
    assert n * BASE**2 < 2**53 and 4 * n * BASE**2 < 2**62
    result = np.zeros((n, n), dtype=object)
    for degree in range(6, -1, -1):
        limb = np.zeros((n, n), dtype=np.int64)
        for i in range(max(0, degree - 3), min(3, degree) + 1):
            limb += (x[i] @ y[degree - i]).astype(np.int64)
        result = result * BASE + limb.astype(object)
    return result


def kernel(h):
    n = len(h)
    eye = np.eye(n, dtype=int)
    doubled = np.rint(2 * h).astype(np.int64)
    assert np.array_equal(doubled, 2 * h)
    doubled -= int(np.rint(np.trace(doubled) / n)) * eye
    rows = [
        [(int(doubled[i, j]), j) for j in np.flatnonzero(doubled[i])] for i in range(n)
    ]
    scale = 4 * 2**SQUARES
    radius = int(np.abs(doubled).sum(axis=1).max())
    assert radius <= scale
    # exp(-i*2H/(4*2^s)) before squaring. An irrelevant scalar phase was removed.
    den = factorial(DEGREE) * scale**DEGREE
    powers = eye.astype(object)
    a, b = den * powers, np.zeros_like(powers)
    for k in range(1, DEGREE + 1):
        powers = np.array(
            [sum(coef * powers[j] for coef, j in row) for row in rows], dtype=object
        )
        weight = (den // (factorial(k) * scale**k)) * (1 if k % 4 in (0, 3) else -1)
        if k % 2:
            b += weight * powers
        else:
            a += weight * powers
    a = (a * UNIT + den // 2) // den
    b = (b * UNIT + den // 2) // den
    # Taylor's integral remainder on the real spectrum and complex entry rounding.
    error = Fraction(radius, scale) ** (DEGREE + 1) / factorial(DEGREE + 1) + Fraction(
        n, UNIT
    )
    # Work with a simple dyadic majorant after each squaring to avoid enormous fractions.
    for _ in range(SQUARES):
        x, y = split(a), split(b)
        ab = product(x, y)
        a = (product(x, x) - product(y, y) + UNIT // 2) // UNIT
        b = (ab + ab.T + UNIT // 2) // UNIT
        assert np.array_equal(a, a.T) and np.array_equal(b, b.T)
        f = UNIT * (2 * error + error**2) + n
        error = Fraction((f.numerator + f.denominator - 1) // f.denominator, UNIT)
    row_error = 2 * error + error**2
    assert row_error < Fraction(2, 10**19)
    return a * a + b * b, row_error
