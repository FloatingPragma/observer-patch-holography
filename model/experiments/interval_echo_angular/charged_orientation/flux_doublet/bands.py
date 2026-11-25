"""Certified band multiplicities using dyadic enclosures and integer inertia."""

from math import isqrt

import numpy as np

GRID = 2**24


def integer_roots():
    # sin(pi/10)=(sqrt(5)-1)/4, cos(pi/10)^2=(5+sqrt(5))/8.
    s5 = isqrt(5 * GRID**2)
    c = isqrt(((5 * GRID + s5) * GRID) // 8)
    s = (s5 - GRID) // 4
    x, y = GRID, 0
    roots = []
    for _ in range(20):
        roots.append((x, y))
        x, y = (x * c - y * s + GRID // 2) // GRID, (x * s + y * c + GRID // 2) // GRID
    return roots


def inertia(matrix):
    """Sign changes of the exact leading principal minors (no zero pivots)."""
    work = matrix.copy()
    previous = 1
    signs = [1]
    for i in range(len(matrix)):
        pivot = work[i, i]
        assert pivot != 0
        signs.append(1 if pivot > 0 else -1)
        for j in range(i + 1, len(matrix)):
            for k in range(i + 1, len(matrix)):
                numerator = pivot * work[j, k] - work[j, i] * work[i, k]
                assert numerator % previous == 0
                work[j, k] = numerator // previous
        previous = pivot
    return sum(a != b for a, b in zip(signs, signs[1:]))


def band_counts(link, seams):
    n = len(link)
    real = np.eye(n, dtype=object) * (85 * GRID)
    imag = np.zeros((n, n), dtype=object)
    roots = integer_roots()
    for a, b in seams:
        c, s = roots[int(link[a, b]) % 20]
        real[b, a] = real[a, b] = -c
        imag[b, a], imag[a, b] = -s, s
    doubled = np.block([[real, -imag], [imag, real]])
    rows = []
    # Relative to the full H approximation error < 250/GRID < 1/1024.
    for cut in (41, 43):
        counts = []
        for side in (-1, 1):
            shifted = doubled - np.eye(2 * n, dtype=object) * (
                2 * GRID * cut + side * (2 * GRID // 1024)
            )
            neg = inertia(shifted)
            assert neg % 2 == 0
            counts.append(neg // 2)
        assert counts[0] == counts[1]
        rows.append(dict(cut=cut, counts_below=counts[0]))
    return dict(
        grid=GRID,
        operator_error_bound=250 / GRID,
        bracket_radius=1 / 1024,
        inertia=rows,
    )
