"""A single-dwell rank witness from exact leading minors of a matrix power series."""

from math import factorial
import json

import numpy as np

from . import fixed

P = fixed.PRIME
MAX_DEGREE = 17


def pivot_minor(matrix):
    a = matrix.copy()
    row_indices = np.arange(a.shape[0])
    pivots = []
    rank = 0
    for j in range(a.shape[1]):
        good = np.flatnonzero(a[rank:, j])
        if not len(good):
            continue
        pivot = rank + int(good[0])
        a[[rank, pivot]] = a[[pivot, rank]]
        row_indices[[rank, pivot]] = row_indices[[pivot, rank]]
        a[rank] = a[rank] * pow(int(a[rank, j]), -1, P) % P
        a[rank + 1 :] = (a[rank + 1 :] - a[rank + 1 :, j, None] * a[rank]) % P
        pivots.append(j)
        rank += 1
    return row_indices[:rank], np.asarray(pivots, dtype=int)


def inverse(matrix):
    n = len(matrix)
    work = np.hstack([matrix, np.eye(n, dtype="int64")])
    for i in range(n):
        good = i + int(np.flatnonzero(work[i:, i])[0])
        work[[i, good]] = work[[good, i]]
        pivot = work[i] * pow(int(work[i, i]), -1, P) % P
        work = (work - work[:, i, None] * pivot) % P
        work[i] = pivot
    return work[:, n:]


def multiply(x, y):
    # All entries are residues in [0,65520], and inner dimension <=400.
    # Every binary64 multiply and sum is therefore an exact integer (<2^53).
    assert x.shape[1] * P * P < 2**53
    return (x.astype(float) @ y.astype(float)).astype("int64") % P


def stages(series):
    original_rows = np.arange(len(series[0]))
    picked = []
    counts = []
    for degree in range(6, MAX_DEGREE + 1):
        r, c = pivot_minor(series[0])
        counts.append([degree, len(c)])
        picked.extend(original_rows[r].tolist())
        if len(c) == series[0].shape[1]:
            return counts, picked
        if not len(c):
            series = series[1:]
            continue
        rr = np.setdiff1d(np.arange(series[0].shape[0]), r)
        cc = np.setdiff1d(np.arange(series[0].shape[1]), c)
        original_rows = original_rows[rr]
        a = [m[np.ix_(r, c)] for m in series]
        d = [m[np.ix_(rr, cc)] for m in series]
        b = [m[np.ix_(r, cc)] for m in series]
        bottom = [m[np.ix_(rr, c)] for m in series]
        inv = inverse(a[0])
        quot = []
        reduced = []
        np.testing.assert_array_equal(
            multiply(inv, a[0]), np.eye(len(c), dtype="int64")
        )
        for n in range(len(series)):
            target = b[n].copy()
            for k in range(1, n + 1):
                target -= multiply(a[k], quot[n - k])
            quot.append(multiply(inv, target % P))
            schur = d[n].copy()
            for k in range(n + 1):
                schur -= multiply(bottom[k], quot[n - k])
            reduced.append(schur % P)
        assert np.all(reduced[0] == 0)
        series = reduced[1:]
    raise ValueError("coefficient range exhausted before full column rank")


def witness():
    b = fixed.camera.front.optics.pairing.formation.federation(1)
    route = fixed.camera.front.optics.pairing.formation.region(b)
    rings = fixed.camera.screens(b, route)
    rows = {}
    for odd in (False, True):
        powers = {}
        for q in (-1, 1):
            h = fixed.camera.front.native.generator(b, route, q, odd)
            powers[q] = fixed.powers(h, rings, MAX_DEGREE)
        # Full charge-resolved analyzer freedom is the strongest case.
        series = []
        for n in range(6, MAX_DEGREE + 1):
            m = np.column_stack(
                [fixed.constraint_rows(powers[q], n)[1] for q in (-1, 1)]
            )
            series.append(m * pow(factorial(n), -1, P) % P)
        counts, picked = stages(series)
        assert len(set(picked)) == 320
        val = sum(n * r for n, r in counts)
        assert val < P
        rows[str(odd)] = dict(
            pivot_multiplicities=counts, leading_power=val, pivot_rows=picked
        )
    return dict(prime=P, maximum_jet=MAX_DEGREE, sectors=rows)


if __name__ == "__main__":
    print(json.dumps(witness(), indent=2))
