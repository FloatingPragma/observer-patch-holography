"""Source-only prediction of the initial current through every exterior edge."""

import json
from hashlib import sha256
import argparse

import numpy as np
from scipy import sparse

from . import mediated_release as lab

HERE = lab.HERE


def edge_certificate(board, marks, i, j, h):
    n = board.n
    rows = np.arange(len(i))
    wave = sparse.csr_matrix(
        (np.ones(len(i), dtype=np.int64), (rows, j)), shape=(len(i), n)
    )
    twice = (2 * h).astype(np.int64)
    assert (twice != 2 * h).nnz == 0
    hw = (twice @ wave).tocsr()
    # The coordinate action of 2H on every conditional release W_z has
    # one free bath Laplacian plus twice the native contrast for that z.
    _, lap, m = lab.cold.motion(board, marks, 1, lab.K)
    adj = (sparse.diags(np.asarray(lap.diagonal()), format="csr") - lap).astype(
        np.int64
    )
    pair = np.asarray(adj[i, j]).ravel()
    save = lab.K * (4 * m[i] + m[j] + 4 * pair)
    expected = lap.astype(np.int64)[j, :] - wave.multiply((2 * save)[:, None])
    assert (hw != expected).nnz == 0

    # Sum all relative initial currents for each ordinary address z on
    # each physical directed seam. Only the two original bath tags can
    # contribute in each configuration pair, so sparse products suffice.
    def index(a, b):
        return a * (n - 1) + b - (b > a)

    carrier = np.arange(n)
    selected = [(u, v) for u, v in board.seams if (u, v) not in marks]
    all_indices = []
    all_targets = []
    for u, v in selected:
        tags = carrier[(carrier != u) & (carrier != v)]
        all_indices.append(np.r_[index(u, tags), index(u, v)])
        all_targets.append(np.r_[index(v, tags), index(v, u)])
    left = np.asarray(all_indices).ravel()
    right = np.asarray(all_targets).ravel()
    raw = wave[left].multiply(hw[right]) - wave[right].multiply(hw[left])
    add = sparse.csr_matrix(
        (
            np.ones(len(left), dtype=np.int64),
            (np.repeat(np.arange(len(selected)), n - 1), np.arange(len(left))),
        )
    )
    observed = (add @ raw).toarray()
    a, b = np.array(selected).T
    ident = np.eye(n, dtype=np.int64)
    derived = 8 * lab.K * (m[a] - m[b])[:, None] * (1 - ident[a] - ident[b])
    derived += 8 * lab.K * (adj[a] - adj[b]).toarray() + (8 * lab.K + 1) * (
        ident[a] - ident[b]
    )
    np.testing.assert_array_equal(observed, derived)
    return dict(
        live_edges=len(selected), cohorts=n, integer_edge_equalities=derived.size
    )


def predict(record):
    board = lab.cold.federation(0)
    marks, _ = lab.cold.read_map(board, record["spectral_tape"])
    shells = lab.cold.distances(board, [x for edge in marks for x in edge])
    gain = np.array(record["decoded"]["ordinary_neighbor_gain"])
    shape = np.zeros(board.n)
    term_banks = []
    for name, mask in (("blank", []), ("written", marks)):
        c = record["source_counts"][name]
        p = np.array(c) / sum(c)
        adj = board.adjacency().tolil()
        for a, b in mask:
            adj[a, b] = 0
            adj[b, a] = 0
        adj = adj.tocsr()
        g = gain.mean()
        local_field = (g - 0.5) * (adj @ p) + g * p
        shape += (1 if name == "written" else -1) * local_field
        # On any exterior cut both source sectors have the same incident
        # neighbor lists. Each surface contraction is a single histogram
        # mean. Bound it directly, avoiding a separate error per edge.
        term_banks.append((p, adj))
    intervals = []
    for r in range(1, shells.max() + 1):
        divergence = np.zeros(board.n)
        for u, v in board.seams:
            if (shells[u] <= r) != (shells[v] <= r):
                sign = 1 if shells[u] <= r else -1
                divergence[u] += sign
                divergence[v] -= sign
        candidates = []
        for g in gain:
            ci = []
            for p, adj in term_banks:
                kernel = ((g - 0.5) * (adj @ divergence) + g * divergence) / (
                    board.n - 1
                )
                center = kernel @ p
                half = np.ptp(kernel) * np.sqrt(
                    lab.RISK / (2 * sum(record["source_counts"]["blank"]))
                )
                ci.append(center + np.array([-1, 1]) * half)
            candidates.append(lab.delta(ci[1], ci[0]))
        intervals.append([min(v[0] for v in candidates), max(v[1] for v in candidates)])
    return dict(
        exterior_cumulative_change=intervals,
        shell_strength=np.bincount(shells, weights=shape).tolist(),
    )


def run():
    source = HERE / "mediated_release.json"
    record = json.loads(source.read_bytes())
    board = lab.cold.federation(0)
    marks, _ = lab.cold.read_map(board, record["spectral_tape"])
    audits = {}
    for name, mask in (("blank", []), ("written", marks)):
        h, i, j = lab.ordinary.two_loads(board, lab.K, mask)
        audits[name] = edge_certificate(board, mask, i, j, h)
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (source, HERE / "CLOCK_TO_FLIGHT.md")
        },
        native=audits,
        from_counts=predict(record),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "medium_field.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        lab.cold.apparatus.check(json.loads(path.read_bytes()), out)
    print(out["native"], out["from_counts"]["exterior_cumulative_change"][:2])
