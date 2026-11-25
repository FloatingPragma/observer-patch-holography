"""Isolate the source-selected pair bright eigenvalue with rational arithmetic."""

from fractions import Fraction
from math import prod

import numpy as np
import sympy as sp

from model.candidate_laws.interval_echo import currents
from model.executable_model import cell


def evaluate(coefficients, interval):
    low = high = Fraction(0)
    for coefficient in reversed(coefficients):
        products = [prod((a, b)) for a in (low, high) for b in interval]
        low = min(products) + int(coefficient)
        high = max(products) + int(coefficient)
    return low, high


def pair_bright():
    """Apply F(K)/(K-t) to the native odd star; F vanishes on its cyclic space."""
    car = cell()
    sector = currents.PairSector(car)
    p = currents.SOURCE_PORT
    opposite = currents._geometry(car)[0][p]
    sharp = np.zeros(sector.P, dtype=np.int64)
    for port, sign in ((p, 1), (opposite, -1)):
        for neighbour in car.neighbours[port]:
            sharp[sector.index(port, neighbour)] = sign
    operator = 2 * sector.hamiltonian(currents.KAPPA).toarray()
    k = operator.astype(np.int64)
    np.testing.assert_array_equal(k, operator)
    poly = sp.Poly.from_list([1, 34, -128, -6904, 41936], sp.Symbol("t"))
    coeff = list(map(int, poly.all_coeffs()))
    vectors = [sharp.astype(object)]
    for _ in range(4):
        vectors.append(k @ vectors[-1])
    if np.any(sum(a * v for a, v in zip(coeff, reversed(vectors), strict=True))):
        raise ValueError("native source has spectrum outside the enrolled quartic")
    interval, mult = poly.intervals(eps=sp.Rational(1, 10**25))[0]
    lo, hi = map(Fraction, interval)
    if not -28 < lo < hi < -27 or mult != 1:
        raise ValueError("lowest enrolled pair line is not isolated")
    projected = np.array(
        [sum(coeff[j] * vectors[3 - i - j] for j in range(4 - i)) for i in range(4)],
        dtype=object,
    ).T
    number = np.isin(sector.aa, [p, opposite]).astype(int) + np.isin(
        sector.bb, [p, opposite]
    )
    squares = np.array([np.convolve(row, row) for row in projected], dtype=object)
    den = evaluate(squares.sum(axis=0), (lo, hi))
    num = evaluate(number @ squares, (lo, hi))
    if min(*den, *num) <= 0:
        raise ValueError("spectral projection has zero or unresolved norm")
    lower, upper = num[0] / den[1], num[1] / den[0]
    if not Fraction(719927, 1000000) < lower <= upper < Fraction(719929, 1000000):
        raise ValueError("bright eigenvalue outside the reported interval")
    return lower, upper


if __name__ == "__main__":
    print("Pair bright eigenvalue:", *map(float, pair_bright()))
