"""A deterministic outward bound for both native connected-radius signs at time four."""

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import factorial, isqrt

import numpy as np
from scipy import sparse

from . import measure as m


GRID = 2**47
ORDER = 13
BILLION = 10**9


def integer_generator(board, seam, n):
    law = m.source.LAW
    local = m.source.formation.emission.four_h(8)
    memories = [1 + i for i, (_, d) in enumerate(law.LABELS) if abs(d) == 1]
    tags = {
        z: (j, law.SRC[k])
        for j, memory in enumerate(memories)
        for k in range(9)
        if max(law.SRC[k]) <= 1
        for z in [9 * memory + k]
    }
    words = list(combinations(range(board.n), n))
    pos = {w: i for i, w in enumerate(words)}
    rows = []
    cols = []
    values = []

    def emit(a, b, v):
        if v:
            rows.append(a)
            cols.append(b)
            values.append(int(v))

    for i, w in enumerate(words):
        occupied = set(w)
        private = 32 * sum(6 - int(board.degree[p]) for p in w)
        ordinary = 0
        for a, b in board.seams:
            if (a, b) == seam or ((a in occupied) == (b in occupied)):
                continue
            ordinary += 34
            j = pos[tuple(sorted(occupied ^ {a, b}))]
            for mem in range(4):
                emit(4 * j + mem, 4 * i + mem, -2)
        pair = tuple(int(p in occupied) for p in seam)
        outside = occupied - set(seam)
        for mem, label in enumerate(memories):
            src = 4 * i + mem
            emit(src, src, private + ordinary)
            column = local[:, 9 * label + law.SRC.index(pair)]
            for j in np.flatnonzero(column):
                other, word = tags[int(j)]
                assert sum(word) == sum(pair)
                out = tuple(
                    sorted(outside | {p for p, b in zip(seam, word, strict=True) if b})
                )
                emit(4 * pos[out] + other, src, column[j])
    g = sparse.csr_matrix(
        (np.array(values, dtype=np.int64), (rows, cols)),
        shape=(4 * len(words), 4 * len(words)),
    )
    g.eliminate_zeros()
    assert (g != g.T).nnz == 0
    simulator, w, _ = m.fp.generator(board, seam, n, 8, True)
    assert w == words and (g != 4 * simulator).nnz == 0
    return g, words


def distance_from_seam(board, seam):
    distances = {i: 0 for i in seam}
    frontier = list(seam)
    while frontier:
        nxt = []
        for p in frontier:
            for q in board.neighbours[p]:
                if q not in distances:
                    distances[q] = distances[p] + 1
                    nxt.append(q)
        frontier = nxt
    return [distances[i] for i in range(board.n)]


def mirror(board, seam, origins, words, g, radius):
    # Antipodal addresses in each of the original 12-port cells and swap owners.
    port = [12 * (1 - (p // 12)) + ((p % 12) ^ 3) for p in range(board.n)]
    assert [port[p] for p in seam] == list(reversed(seam))
    assert [port[p] for p in origins] == list(reversed(origins))
    assert {tuple(sorted((port[a], port[b]))) for a, b in board.seams} == set(
        board.seams
    )
    pos = {w: i for i, w in enumerate(words)}
    perm = [
        4 * pos[tuple(sorted(port[p] for p in w))] + (mem ^ 2)
        for w in words
        for mem in range(4)
    ]
    assert (g != g[perm, :][:, perm]).nnz == 0
    np.testing.assert_array_equal(radius, radius[perm])
    return perm


def bound_columns(g, indices):
    diag = g.diagonal()
    scalar = int((diag.min() + diag.max()) // 2)
    centered = g - scalar * sparse.eye(g.shape[0], dtype=np.int64, format="csr")
    norm = int(abs(centered).sum(axis=1).max())
    steps = 1 << (2 * norm - 1).bit_length()
    x = Fraction(norm, steps)
    assert x <= Fraction(1, 2)
    # At total time 4, exp(-4iH)=exp(-i(4H)). Dropping a common phase
    # preserves every diagonal and interference expectation.
    rounding = Fraction(isqrt(2 * g.shape[0]) + 1, 2 * GRID)
    tail = x ** (ORDER + 1) / factorial(ORDER + 1)
    error = steps * (2 * tail + 2 * ORDER * rounding)
    assert error < 1
    real = np.zeros((g.shape[0], len(indices)), dtype=np.int64)
    imag = np.zeros_like(real)
    real[indices, np.arange(len(indices))] = GRID
    for _ in range(steps):
        term_real = real.copy()
        term_imag = imag.copy()
        for k in range(1, ORDER + 1):
            den = steps * k
            assert (
                norm * max(int(abs(term_real).max()), int(abs(term_imag).max())) + den
                < 2**63
            )
            term_real, term_imag = (
                (centered @ term_imag + den // 2) // den,
                (-(centered @ term_real) + den // 2) // den,
            )
            assert (
                max(int(abs(real).max()), int(abs(imag).max()))
                + max(int(abs(term_real).max()), int(abs(term_imag).max()))
                < 2**63
            )
            real += term_real
            imag += term_imag
    return (
        real,
        imag,
        error,
        dict(
            scalar_removed_from_4H=scalar,
            absolute_row_bound=norm,
            steps=steps,
            per_column_error=str(error),
            per_column_error_ceiling_billionths=ceil(error * BILLION),
        ),
    )


def ceil(x):
    return -((-x.numerator) // x.denominator)


def bounds(center, error):
    lower = BILLION * (center - error)
    upper = BILLION * (center + error)
    return [lower.numerator // lower.denominator, ceil(upper)]


def run():
    board, seam, dist, origins = m.fp.make_board()
    exact_dist = distance_from_seam(board, seam)
    np.testing.assert_array_equal(dist, exact_dist)
    densities = {
        "target": np.diag([1, -1, -1, 1]),
        "phase_feedback": np.zeros((4, 4), dtype=int),
    }
    densities["phase_feedback"][1, 2] = densities["phase_feedback"][2, 1] = -1
    for name, delta4 in densities.items():
        lo, hi = m.stocks()[name]
        np.testing.assert_array_equal(delta4, 4 * (hi[1] - lo[1]))
        np.testing.assert_array_equal(delta4, delta4[::-1, ::-1])
    values = {name: [] for name in densities}
    sectors = []
    for total in (1, 2):
        g, words = integer_generator(board, seam, total)
        extra = () if total == 1 else (origins[0],)
        other = () if total == 1 else (origins[1],)

        def sources(ports):
            return [
                4 * words.index(tuple(sorted((p, *ports)))) + mem
                for mem in (0, 2)
                for p in seam
            ]

        ids = sources(extra)
        radius = np.repeat([sum(exact_dist[p] for p in w) for w in words], 4)
        perm = mirror(board, seam, origins, words, g, radius)
        assert [perm[i] for i in ids] == sources(other)[::-1]
        re, im, error, metadata = bound_columns(g, ids)
        # Arbitrary-precision contractions; only Re(U* R U) enters real states.
        rr = re.astype(object)
        ii = im.astype(object)
        r = radius[:, None].astype(object)
        gram = rr.T @ (r * rr) + ii.T @ (r * ii)
        for name, delta4 in densities.items():
            center = Fraction(int((delta4 * gram.T).sum()), 4 * GRID**2)
            bound = (
                Fraction(int(abs(delta4).sum()) * int(radius.max()), 4)
                * error
                * (2 + error)
            )
            values[name].append((center, bound))
        sectors.append(
            dict(
                total=total,
                dimension=g.shape[0],
                source_indices=ids,
                observable_norm=int(radius.max()),
                **metadata,
            )
        )
    result = {}
    for name, ((small, err1), (large, err2)) in values.items():
        center = large - small
        error = err1 + err2
        result[name] = dict(
            fixed_point_center=str(center),
            absolute_error_bound=str(error),
            radius_billionths=bounds(center, error),
        )
    assert (
        result["target"]["radius_billionths"][0] > 0
        and result["phase_feedback"]["radius_billionths"][1] < 0
    )
    return dict(
        protocol_sha256=sha256(
            (m.HERE / "TIME_FOUR_CERTIFICATE.md").read_bytes()
        ).hexdigest(),
        time=4,
        launch=origins[0],
        exact_reflected_launch=origins[1],
        scale=GRID,
        order=ORDER,
        sectors=sectors,
        contrasts=result,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = run()
    path = m.HERE / "time_four.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        assert json.loads(path.read_text()) == data
    print(data["sectors"])
    print({k: v["radius_billionths"] for k, v in data["contrasts"].items()})
