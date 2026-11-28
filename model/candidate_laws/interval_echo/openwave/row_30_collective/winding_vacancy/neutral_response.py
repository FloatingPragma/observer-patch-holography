"""Finite-momentum neutral density spectroscopy in a homogeneous winding sector."""

import argparse
import json
from hashlib import sha256

import networkx as nx
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh

from . import geometry, momentum, reservoir

HERE = reservoir.HERE
LOCK = "8303071a19ce2ab1814bcef2ac3b851ca4254fbbbcb59041f42ca7c9fbecf19c"


def uniform(n):
    L = n // 2
    records = []
    for heads in range(1 << L):
        bits = [(heads >> i) & 1 for i in range(L)]
        occ = np.array([heads << n], dtype=np.int64)
        for i in range(L):
            N = 1 + bits[i - 1] - bits[i]
            if N == 2:
                occ |= 3 << (2 * i)
            elif N == 1:
                occ = np.r_[occ | (1 << (2 * i)), occ | (2 << (2 * i))]
        records.append(occ)
    all_words = np.sort(np.concatenate(records))
    assert len(all_words) == 3**L + 1
    assert np.all(momentum.charges(all_words, n) == 1)
    return all_words


def translations(words, n):
    L = n // 2
    r = momentum.rotate(words & ((1 << n) - 1), 2, n) + (
        momentum.rotate(words >> n, 1, L) << n
    )
    shifts = np.searchsorted(words, r)
    seen = np.zeros(len(words), dtype=bool)
    cycles = []
    for row in range(len(words)):
        if seen[row]:
            continue
        seq = [row]
        nxt = shifts[row]
        while nxt != row:
            seq.append(nxt)
            nxt = shifts[nxt]
        seen[seq] = True
        cycles.append(seq)
    return cycles


def fourier(cycles, nwords, j, L):
    r, c, v = [], [], []
    col = 0
    for seq in cycles:
        if len(seq) * j % L:
            continue
        r.extend(seq)
        c.extend([col] * len(seq))
        v.extend(np.exp(-2j * np.pi * j * np.arange(len(seq)) / L) / np.sqrt(len(seq)))
        col += 1
    w = coo_matrix((v, (r, c)), shape=(nwords, col)).tocsr()
    if j == 0 or 2 * j == L:
        w = w.real
    return w


def low_modes(matrix, n=10):
    if matrix.shape[0] < 100:
        vals, vecs = eigh(matrix.toarray())
        return vals[:n], vecs[:, :n]
    v0 = np.random.default_rng(318144).standard_normal(matrix.shape[0])
    vals, vecs = eigsh(matrix, k=n + 2, which="SA", tol=2e-13, ncv=48, v0=v0)
    # Some real-symmetric spectra have exact degeneracies. Orthogonalize their
    # invariant span before measuring summed weight, irrespective of solver phase.
    basis = np.linalg.qr(vecs)[0]
    vals, vec = eigh(basis.conj().T @ (matrix @ basis))
    return vals[:n], (basis @ vec)[:, :n]


def measure_ring(route, tables):
    n = len(route)
    L = n // 2
    words = uniform(n)
    h = reservoir.compile_words(n, words, tables)
    assert connected_components(h, directed=False)[0] == 1
    if n == 10:
        all_words, full = reservoir.integer_sector(n, L, tables)
        ix = np.searchsorted(all_words, words)
        assert (h - full[np.ix_(ix, ix)]).nnz == 0
        prior = json.loads((HERE / "reservoir.json").read_text())["rings"][-1]
        benchmark = prior["sector_minima"][L]
    cycles = translations(words, n)
    rows = []
    for j in range(L // 2 + 1):
        k = 2 * np.pi * j / L
        W = fourier(cycles, len(words), j, L)
        block = (W.getH() @ h @ W) / 8
        vals, v = low_modes(block)
        full = W @ v
        error = np.max(np.linalg.norm(h @ full / 8 - full * vals, axis=0))
        assert error < 1e-8
        if j == 0:
            base = vals[0]
            g = full[:, 0]
            if n == 10:
                assert abs(base - benchmark) < 1e-10
        gap = vals - base
        bands = []
        for parity in [0, 1]:
            f = (
                sum(
                    np.exp(-1j * k * c) * ((words >> (2 * c + parity)) & 1)
                    for c in range(L)
                )
                / L
            )
            src = f * g
            elastic = abs(g.conj() @ src) ** 2
            connected = src - g * (g.conj() @ src)
            total = np.linalg.norm(connected) ** 2
            np.testing.assert_allclose(
                np.linalg.norm(W.getH() @ connected) ** 2, total, atol=1e-12
            )
            weights = abs(full.conj().T @ connected) ** 2
            # Store the combined weight of equal energies; individual vectors
            # inside a degenerate eigenspace have no physically fixed identity.
            groups = []
            for e, w in zip(gap, weights, strict=True):
                if not groups or abs(e - groups[-1][0]) > 1e-8:
                    groups.append([float(e), float(w)])
                else:
                    groups[-1][1] += w
            lines = [[e, w / total] for e, w in groups] if total > 1e-14 else []
            candidate = next((e for e, weight in lines if weight >= 0.02), None)
            bands.append(
                dict(
                    parity=parity,
                    elastic=elastic,
                    inelastic=total,
                    energy_and_fraction=lines,
                    residual_weight=max(0.0, 1 - sum(w for _, w in lines)),
                    first_visible=candidate,
                )
            )
        rows.append(
            dict(
                index=j,
                k=k,
                dimension=block.shape[0],
                lowest_energies_above_ground=gap.tolist(),
                matrix_residual=error,
                channels=bands,
            )
        )
    return dict(
        n=n, ports=route, uniform_words=len(words), sector_ground=base, response=rows
    )


def rounded(v):
    if isinstance(v, dict):
        return {k: rounded(x) for k, x in v.items()}
    if isinstance(v, list):
        return [rounded(x) for x in v]
    if isinstance(v, (float, np.floating)):
        return round(float(v), 10) + 0.0
    return v


def run():
    assert sha256((HERE / "NEUTRAL_RESPONSE.md").read_bytes()).hexdigest() == LOCK
    board = geometry.federation(2)
    graph = nx.Graph(board.seams)
    _, reps = geometry.cycles(graph)
    rings = reservoir.boards()
    routes = [rings[0][0], rings[-1][0], reps[18], reps[20]]
    assert all(geometry.qualify(board, graph, r)["qualified"] for r in routes)
    tables = reservoir.native_contacts()
    return dict(
        protocol_sha256=LOCK, rings=[rounded(measure_ring(r, tables)) for r in routes]
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    out = json.dumps(run(), sort_keys=True, indent=2) + "\n"
    target = HERE / "neutral_response.json"
    if opt.write:
        target.write_text(out)
    elif target.read_text() != out:
        raise AssertionError("Neutral response replay mismatch")
    print("Verified homogeneous neutral density response at every registered momentum.")
