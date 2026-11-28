"""The original radius read on the clock-equipped board at its registered time."""

import argparse
from functools import cache
from hashlib import sha256
from itertools import combinations, product
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from . import signal as s


@cache
def parts(n, stiffness=8):
    board, seam, dist, origins = s.m.fp.make_board()
    marked = (seam, (0, 1), (2, 3))
    local = s.local_four_h(stiffness)
    words = list(combinations(range(board.n), n))
    ix = {w: i for i, w in enumerate(words)}
    pairs = ((0, 0), (1, 0), (0, 1), (1, 1))
    rows = [[] for _ in range(4)]
    cols = [[] for _ in range(4)]
    vals = [[] for _ in range(4)]

    def add(k, i, j, v):
        if v:
            rows[k].append(i)
            cols[k].append(j)
            vals[k].append(int(v))

    for i, w in enumerate(words):
        occupied = set(w)
        private = 4 * stiffness * sum(6 - int(board.degree[p]) for p in w)
        diag = private
        for e in board.seams:
            if e in marked or len(occupied & set(e)) != 1:
                continue
            diag += 4 * stiffness + 2
            j = ix[tuple(sorted(occupied ^ set(e)))]
            for cc in range(8):
                add(0, 8 * j + cc, 8 * i + cc, -2)
        for cc in range(8):
            add(0, 8 * i + cc, 8 * i + cc, diag)
        for k, e in enumerate(marked):
            mask = 1 << k
            outside = occupied - set(e)
            pw = pairs.index(tuple(int(p in occupied) for p in e))
            for cc in range(8):
                for where, op in ((0, local[0]), (k + 1, local[1] - local[0])):
                    col = op[:, 2 * pw + ((cc >> k) & 1)]
                    for j in np.flatnonzero(col):
                        tw, tc = divmod(int(j), 2)
                        y = ix[
                            tuple(
                                sorted(
                                    outside
                                    | {
                                        p
                                        for p, b in zip(e, pairs[tw], strict=True)
                                        if b
                                    }
                                )
                            )
                        ]
                        add(where, 8 * y + (cc & ~mask) + mask * tc, 8 * i + cc, col[j])
    components = [
        sparse.csr_matrix(
            (np.array(v, dtype=int), (r, c)), shape=(8 * len(words), 8 * len(words))
        )
        for v, r, c in zip(vals, rows, cols, strict=True)
    ]
    for op in components:
        assert (op != op.T).nnz == 0
    starts = {}
    for extra in (None, *origins):
        count = n - int(extra is not None)
        if count not in (0, 1):
            continue
        sites = [()] if not count else [(p,) for p in seam]
        starts[extra] = [
            8 * ix[tuple(sorted((*w, *(() if extra is None else (extra,)))))]
            + c
            + 2 * nc
            + 4 * fc
            for c in (0, 1)
            for w in sites
            for nc, fc in product(range(2), repeat=2)
        ]
    return components, starts, np.repeat([sum(dist[p] for p in w) for w in words], 8)


def run():
    stocks, _ = s.get_arms()
    radius = {}
    residual = 0.0
    for n in (0, 1, 2):
        comp, starts, r = parts(n)
        all_ix = sorted({i for ix in starts.values() for i in ix})
        initial = np.zeros((len(r), len(all_ix)))
        initial[all_ix, np.arange(len(all_ix))] = 1
        for par in range(8):
            h = comp[0].copy()
            for k in range(3):
                if par & (1 << k):
                    h += comp[k + 1]
            if n <= 1:
                np.testing.assert_array_equal(h.toarray(), s.model(n, par)[0])
            center = int(np.mean(h.diagonal()))
            reduced = -1j * (h - center * sparse.eye(len(r), format="csr"))
            # h = 4H so this evolves to time four; global phases drop out
            # separately within each parity block for the matter-only reader.
            psi = expm_multiply(reduced, initial, traceA=reduced.diagonal().sum())
            residual = max(
                residual, float(abs(psi.conj().T @ psi - np.eye(len(all_ix))).max())
            )
            for extra, indices in starts.items():
                count = n - int(extra is not None)
                block = psi[:, [all_ix.index(j) for j in indices]]
                gram = (block.conj().T @ (r[:, None] * block)).real
                key = (count, extra)
                radius[key] = radius.get(key, 0) + gram / 8
    assert residual < 1e-10
    _, _, _, origins = s.m.fp.make_board()

    def read(name, extra):
        return sum(
            np.sum(np.kron(v, np.eye(4) / 4) * radius[n, extra])
            for n, v in enumerate(stocks[name])
        )

    out = {}
    for key, (lo, hi) in s.ARMS.items():
        rows = [
            float((read(hi, p) - read(hi, None)) - (read(lo, p) - read(lo, None)))
            for p in origins
        ]
        out[key] = dict(mirror_launches=rows, mean=float(np.mean(rows)))
    return dict(
        protocol_sha256=sha256((s.HERE / "SAME_BOARD.md").read_bytes()).hexdigest(),
        time=4,
        source_near_far=[(1, 14), (0, 1), (2, 3)],
        ports=origins,
        full_dimensions=[8 * len(parts(n)[2]) for n in (0, 1, 2)],
        radius=out,
        max_gram_residual=residual,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = run()
    path = s.HERE / "flight.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2) + "\n")
    else:
        s.m.fp.check(json.loads(path.read_text()), json.loads(json.dumps(out)))
    print(out)
