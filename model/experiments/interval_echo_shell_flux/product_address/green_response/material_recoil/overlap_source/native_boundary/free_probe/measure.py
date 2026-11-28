"""Response to adding one ordinary load to a freely evolving written source."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from .. import free
from ..source import LAW
from ...closed_region.measure import check, layout

HERE = Path(__file__).resolve().parent
TIMES = (0, 1, 4, 16, 64)
PAIRED = np.array([0.5, 0, 0, 0.5])
ERASED = np.full(4, 0.25)


def make_board():
    _, edges, _, _ = layout(2)
    b = free.Carrier("closed two-owner assembly", 24, edges)
    seam = min(b.external)
    dist = (
        shortest_path(
            b.adjacency(), directed=False, unweighted=True, indices=list(seam)
        )
        .min(axis=0)
        .astype(int)
    )
    probes = [
        p
        for owner in range(b.cells)
        for p in range(12 * owner, 12 * owner + 12)
        if dist[p] == max(dist[12 * owner : 12 * owner + 12])
    ]
    return b, seam, dist, probes


def restricted(k):
    full = free.source.formation.emission.four_h(k)
    memories = [1 + i for i, (_, d) in enumerate(LAW.LABELS) if abs(d) == 1]
    basis = [
        9 * m + LAW.SRC.index((a, b)) for a in (0, 1) for b in (0, 1) for m in memories
    ]
    assert not np.any(full[np.ix_(np.setdiff1d(np.arange(99), basis), basis)])
    local = full[np.ix_(basis, basis)] / 4
    local = local.reshape(4, 4, 4, 4)
    for i, j in np.ndindex(4, 4):
        if i.bit_count() != j.bit_count():
            assert not np.any(local[i, :, j, :])
    branches = []
    for co in np.flatnonzero(abs(LAW.echo()[:, LAW.SRC.index((0, 1))]) > 1e-8):
        branches.append((memories.index(int(co // 9)), LAW.SRC[co % 9].index(1)))
    return local, branches


def generator(board, seam, n, k, written):
    local, _ = restricted(k)
    words = list(combinations(range(board.n), n))
    idx = {w: i for i, w in enumerate(words)}
    bits = np.zeros((len(words), board.n), dtype=int)
    for i, w in enumerate(words):
        bits[i, list(w)] = 1
    edges = np.array([s for s in board.seams if not written or s != seam])
    x, y = edges.T
    unlike = bits[:, x] != bits[:, y]
    potentials = (k + 0.5) * unlike.sum(axis=1) + k * (bits @ (6 - board.degree))
    row = []
    col = []
    for i, edge_id in zip(*unlike.nonzero(), strict=True):
        x, y = edges[edge_id]
        dest = tuple(sorted(set(words[i]) ^ {x, y}))
        row.append(i)
        col.append(idx[dest])
    base = sparse.diags(potentials, format="csr") + sparse.csr_matrix(
        (np.full(len(row), -0.5), (row, col)), shape=(len(words), len(words))
    )
    if not written:
        return base, words, bits
    mat = sparse.kron(base, sparse.eye(4), format="csr")
    a, b = seam
    row = []
    col = []
    val = []
    for j, w in enumerate(words):
        outside = set(w) - set(seam)
        src = 2 * bits[j, a] + bits[j, b]
        for dst in range(4):
            block = local[dst, :, src, :]
            if not np.any(block):
                continue
            target = idx[
                tuple(
                    sorted(outside | {p for i, p in enumerate(seam) if dst & (2 >> i)})
                )
            ]
            u, v = block.nonzero()
            row.extend(4 * target + u)
            col.extend(4 * j + v)
            val.extend(block[u, v])
    mat += sparse.csr_matrix((val, (row, col)), shape=mat.shape)
    assert (mat != mat.T).nnz == 0
    return mat, words, bits


def propagation(h, indices):
    start = np.zeros((h.shape[0], len(indices)))
    start[indices, np.arange(len(indices))] = 1
    e, v = eigh(h.toarray())
    coeff = v[indices, :].T
    initial_energy = np.sum(e[:, None] * coeff**2, axis=0)
    residual = 0.0
    waves = []
    state = start.astype(complex)
    center = h.diagonal().mean()
    gen = -1j * (h - center * sparse.eye(h.shape[0], format="csr"))
    prev = 0
    for t in TIMES:
        if t > prev:
            state = expm_multiply((t - prev) * gen, state, traceA=0.0)
        check_state = v @ (np.exp(-1j * (e - center) * t)[:, None] * coeff)
        residual = max(
            residual, float(np.max(np.linalg.norm(check_state - state, axis=0)))
        )
        np.testing.assert_allclose((abs(state) ** 2).sum(axis=0), 1, atol=2e-11)
        np.testing.assert_allclose(
            np.sum(state.conj() * (h @ state), axis=0).real,
            initial_energy,
            rtol=0,
            atol=5e-10,
        )
        waves.append(state.copy())
        prev = t
    assert residual < 1e-10
    return waves, residual


def occupations(h, words, bits, starts, dim):
    ids = [words.index(tuple(sorted(w))) * dim + mem for w, mem in starts]
    waves, err = propagation(h, ids)
    # Actual number at every port, including both identical loads.
    reads = np.stack(
        [bits.T @ (abs(w.reshape(len(words), dim, -1)) ** 2).sum(axis=1) for w in waves]
    )
    np.testing.assert_allclose(reads.sum(axis=1), len(words[0]), atol=1e-11)
    return reads, h.diagonal()[ids], err


def readout(d, dist):
    shell = np.array([d[:, dist == i].sum(axis=1) for i in range(max(dist) + 1)]).T
    np.testing.assert_allclose(shell.sum(axis=1), 0, atol=2e-11)
    return dict(
        port_contrast=d.tolist(),
        shell_contrast=shell.tolist(),
        radial_moment=(d @ dist).tolist(),
        inner_two_shells=shell[:, :2].sum(axis=1).tolist(),
    )


def material(k):
    board, seam, dist, origins = make_board()
    _, branches = restricted(k)
    starts = [((seam[load],), m) for m, _ in branches for _, load in branches]
    h1, w1, b1 = generator(board, seam, 1, k, True)
    if k == 8:
        original = free.operators(board)
        assert (h1 != original[0]).nnz == 0
        np.testing.assert_array_equal(
            original[6].argmax(axis=0), [w1.index(w) * 4 + m for w, m in starts]
        )
    single, e1, err1 = occupations(h1, w1, b1, starts, 4)
    h2, w2, b2 = generator(board, seam, 2, k, True)
    launched = [((*w, p), m) for p in origins for w, m in starts]
    two, e2, err2 = occupations(h2, w2, b2, launched, 4)
    e2 = e2.reshape(len(origins), 4)
    results = []
    for i, port in enumerate(origins):
        panel = two[:, :, 4 * i : 4 * i + 4]
        paired = panel @ PAIRED
        erased = panel @ ERASED
        connected = (erased - single @ ERASED) - (paired - single @ PAIRED)
        out = readout(connected, dist)
        source_delta = float(e1 @ (ERASED - PAIRED))
        assert abs(e2[i] @ (ERASED - PAIRED) - source_delta) < 1e-13
        assert np.ptp(e2[i] - e1) < 1e-13
        out.update(
            origin=port,
            with_extra_paired=paired.tolist(),
            with_extra_erased=erased.tolist(),
            written_source_energy_increment=source_delta,
            inserted_energy_per_unit=float((e2[i] - e1)[0]),
        )
        results.append(out)
    return dict(
        states=[h1.shape[0], h2.shape[0]],
        seam=seam,
        dist=dist.tolist(),
        no_extra_paired=(single @ PAIRED).tolist(),
        no_extra_erased=(single @ ERASED).tolist(),
        max_independent_wave_discrepancy=max(err1, err2),
        probes=results,
    )


def ordinary():
    board, seam, dist, origins = make_board()
    h1, w1, b1 = generator(board, seam, 1, 8, False)
    h2, w2, b2 = generator(board, seam, 2, 8, False)
    one, _, err1 = occupations(h1, w1, b1, [((p,), 0) for p in (*seam, *origins)], 1)
    both, _, err2 = occupations(
        h2, w2, b2, [((s, p), 0) for p in origins for s in seam], 1
    )
    src = one[:, :, :2].mean(axis=2)
    out = []
    for i, p in enumerate(origins):
        pair = both[:, :, 2 * i : 2 * i + 2].mean(axis=2)
        probe = one[:, :, 2 + i]
        row = readout(pair - src - probe, dist)
        row.update(origin=p, source_and_extra=pair.tolist(), extra_alone=probe.tolist())
        out.append(row)
    return dict(
        states=[h1.shape[0], h2.shape[0]],
        source_alone=src.tolist(),
        probes=out,
        max_independent_wave_discrepancy=max(err1, err2),
    )


def run():
    return dict(
        protocol_sha256=sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        times=TIMES,
        written={str(k): material(k) for k in (8, 0)},
        ordinary=ordinary(),
        independent_board_control=0,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = json.loads(json.dumps(run()))
    for tag, rows in [*data["written"].items(), ("ordinary", data["ordinary"])]:
        print(
            tag, [(x["origin"], x["radial_moment"]) for x in rows["probes"]], flush=True
        )
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n")
    else:
        check(json.loads(path.read_text()), data)
