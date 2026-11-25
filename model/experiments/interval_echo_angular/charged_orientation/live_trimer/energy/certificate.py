"""Integer cyclicity and dyadic enclosures for a finite real symmetric generator."""

from fractions import Fraction
from math import isqrt

import numpy as np


def cyclic(matrix, start, prime=65521):
    assert all(prime % k for k in range(2, isqrt(prime) + 1))
    n = len(matrix)
    a = matrix.astype(np.int64)
    v = np.eye(n, dtype=np.int64)[:, start]
    cols = []
    for _ in range(n):
        cols.append(v)
        v = (a @ v) % prime
    w = np.column_stack([*cols, v])
    det = 1
    for j in range(n):
        active = np.flatnonzero(w[j:, j])
        if not len(active):
            raise ArithmeticError("Source is not cyclic modulo the chosen prime")
        pivot = j + active[0]
        if pivot != j:
            w[[j, pivot]] = w[[pivot, j]]
            det = -det
        b = int(w[j, j])
        det = det * b % prime
        w[j] = w[j] * pow(b, -1, prime) % prime
        w[j + 1 :] = (w[j + 1 :] - w[j + 1 :, j, None] * w[j]) % prime
    assert det != 0
    relation = np.zeros(n, dtype=np.int64)
    for j in reversed(range(n)):
        relation[j] = (w[j, -1] - w[j, j + 1 : n] @ relation[j + 1 :]) % prime
    np.testing.assert_array_equal(np.column_stack(cols) @ relation % prime, v)
    return dict(
        prime=prime,
        source_index_in_component=start,
        dimension=n,
        det_mod_prime=int(det),
        monic_charpoly_ascending=[*((-relation) % prime).tolist(), 1],
    )


def ordered_bound(h, values, frame):
    """Weyl bound from exact residual and Gram error of a rounded diagonalizer."""
    scale = 2**46
    ticks = 2**43
    x = np.rint(scale * frame).astype(np.int64)
    d = np.rint(ticks * values).astype(np.int64)
    assert int(np.max(np.sum(abs(h), axis=1))) * int(abs(x).max()) < 2**62
    hx = (h @ x).astype(object)
    x = x.astype(object)
    g = x.T @ x - scale**2 * np.eye(len(h), dtype=object)
    gram = Fraction(int(max(np.sum(abs(g), axis=1))), scale**2)
    residual = hx * ticks - 8 * x * d
    r = max(int(max(np.sum(abs(residual), axis=k))) for k in (0, 1))
    err = Fraction(r, 8 * scale * ticks)
    nh = Fraction(int(np.max(np.sum(abs(h), axis=1))), 8)
    nd = Fraction(int(max(abs(d))), ticks)
    assert gram < Fraction(1, 2)
    bound = err + (nh + nd) * gram
    return (
        dict(
            frame_scale=scale,
            energy_ticks=ticks,
            gram_bound=str(gram),
            residual_bound=str(err),
            ordered_error_bound=str(bound),
        ),
        d,
        bound,
    )
