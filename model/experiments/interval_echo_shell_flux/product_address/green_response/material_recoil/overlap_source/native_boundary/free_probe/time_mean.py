"""Retain degenerate energy coherences and compute exact small-time radial moments."""

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial

import numpy as np
from scipy.linalg import eigh

from . import measure as m


def coefficients(h, radius, indices):
    g = 4 * h
    integer = g.astype(np.int64)
    assert (integer != g).nnz == 0
    integer.setdiag(integer.diagonal() - int(np.median(integer.diagonal())))
    integer.eliminate_zeros()
    columns = [
        list(
            zip(
                integer.indices[integer.indptr[i] : integer.indptr[i + 1]].tolist(),
                integer.data[integer.indptr[i] : integer.indptr[i + 1]].tolist(),
                strict=True,
            )
        )
        for i in range(h.shape[0])
    ]
    radius = radius.tolist()
    results = []
    for src in indices:
        powers = [{src: 1}]
        for _ in range(12):
            new = defaultdict(int)
            for i, c in powers[-1].items():
                for j, k in columns[i]:
                    new[j] += c * k
            powers.append({i: v for i, v in new.items() if v})
        row = []
        for p in range(0, 13, 2):
            terms = 0
            for j in range(p + 1):
                inner = sum(
                    v * powers[p - j].get(i, 0) * radius[i]
                    for i, v in powers[j].items()
                )
                terms += (-1) ** (p // 2 + j) * comb(p, j) * inner
            row.append(Fraction(terms, 4**p * factorial(p)))
        results.append(row)
    return np.array(results, dtype=object)


def sector(board, seam, dist, n, k, written, starts):
    h, words, bits = m.generator(board, seam, n, k, written)
    dim = 4 if written else 1
    idx = [words.index(tuple(sorted(w))) * dim + mem for w, mem in starts]
    radius = np.repeat(bits @ dist, dim)
    e, v = eigh(h.toarray())
    previous = None
    for tol in (1e-8, 1e-10):
        splits = np.flatnonzero(np.diff(e) > tol) + 1
        groups = np.split(np.arange(len(e)), splits)
        mean = np.zeros(len(idx))
        for group in groups:
            states = v[:, group] @ v[idx, :][:, group].T
            mean += radius @ (states**2)
        if previous is not None:
            np.testing.assert_allclose(previous, mean, atol=1e-11, rtol=0)
        previous = mean
    return (
        mean,
        coefficients(h, radius, idx),
        dict(
            groups=len(groups),
            smallest_separation=float(np.diff(e)[splits - 1].min()),
            largest_within_group=float(max(np.ptp(e[g]) for g in groups)),
        ),
    )


def run():
    board, seam, dist, origins = m.make_board()
    result = {}
    for name, k in (("native", 8), ("no_contrast", 0), ("ordinary", 8)):
        written = name != "ordinary"
        _, branches = m.restricted(k)
        source = (
            [((seam[p],), c) for c, _ in branches for _, p in branches]
            if written
            else [((p,), 0) for p in seam]
        )
        unit = source if written else [*source, *[((p,), 0) for p in origins]]
        joint = [((*s, p), c) for p in origins for s, c in source]
        single, c1, g1 = sector(board, seam, dist, 1, k, written, unit)
        pair, c2, g2 = sector(board, seam, dist, 2, k, written, joint)
        rows = []
        for i, port in enumerate(origins):
            if written:
                z = [-1, 1, 1, -1]
                val = float((pair[4 * i : 4 * i + 4] - single) @ z / 4)
                coeff = (z @ (c2[4 * i : 4 * i + 4] - c1)) / 4
            else:
                val = float(
                    pair[2 * i : 2 * i + 2].mean() - single[:2].mean() - single[2 + i]
                )
                coeff = (c2[2 * i] + c2[2 * i + 1] - c1[0] - c1[1]) / 2 - c1[2 + i]
            rows.append(
                dict(
                    port=port,
                    time_mean=val,
                    taylor_even_coefficients=[str(x) for x in coeff],
                )
            )
        result[name] = dict(rows=rows, spectral_groups=[g1, g2])
    return dict(
        protocol_sha256=sha256((m.HERE / "TIME_MEAN.md").read_bytes()).hexdigest(),
        panels=result,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    path = m.HERE / "time_mean.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        m.check(json.loads(path.read_text()), data)
    for tag, v in data["panels"].items():
        print(tag, v["rows"][0])
