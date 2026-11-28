"""Predict release of a separate native double from two ordinary units and their source."""

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from types import SimpleNamespace

import numpy as np
from scipy import sparse

from . import pair_contact as contact
from . import pair_motion as motion

HERE = motion.HERE
BATCH = 2**48
BATCHES = 4096
TRIALS = BATCH * BATCHES


def operators(board, inside, *, normalize=True):
    h, _, i, j, index = motion.hamiltonian(board, [])
    h = h.copy()
    h.setdiag(0)
    h.setdiag(-np.asarray(h.sum(axis=1)).ravel())
    occupancy = inside[i].astype(int) + inside[j]
    N = sparse.diags(occupancy, dtype=float, format="csr")
    comm = h @ N - N @ h
    kinetic = h @ comm - comm @ h
    adj = board.adjacency().toarray().astype(int)
    u, v = np.where(adj & inside[:, None] & ~inside[None, :])
    diagonal = np.zeros(len(i), dtype=int)
    rows, cols, coeff = [], [], []
    for a, b in zip(u, v, strict=True):
        empty = (i != a) & (j != a) & (i != b) & (j != b)
        gradient = adj[a] - adj[b]
        diagonal += 4 * empty * (gradient[i] + gradient[j])
        k = np.flatnonzero(
            (np.arange(board.n) != a) & (np.arange(board.n) != b) & (gradient != 0)
        )
        rows.extend(index[a, k])
        cols.extend(index[b, k])
        coeff.extend(gradient[k])
    cross = sparse.csr_matrix((coeff, (rows, cols)), shape=h.shape)
    contrast = sparse.diags(diagonal, dtype=float, format="csr") + cross + cross.T
    if normalize:
        return kinetic / (board.n - 2), contrast / (board.n - 2)
    return kinetic, contrast


def check_generator():
    # Ring with two chords, irregular degrees and a nonsymmetric boundary.
    size = 7
    seams = [(i, (i + 1) % size) for i in range(size)] + [(1, 4), (0, 3)]
    adj = np.zeros((size, size), dtype=int)
    for a, b in seams:
        adj[a, b] = adj[b, a] = 1
    board = SimpleNamespace(n=size, adjacency=lambda: sparse.csr_matrix(adj))
    ordinary = list(combinations(range(size), 2))
    states = [(d, S) for S in ordinary for d in range(size) if d not in S]
    index = {w: n for n, w in enumerate(states)}
    W = np.zeros((len(states), len(ordinary)), dtype=int)
    exact_refresh = []
    _, _, H8, H0 = contact.stock.kernels()
    for kappa, raw in ((0, H0), (8, H8)):
        H = np.zeros((len(states), len(states)), dtype=int)
        for n, (d, S) in enumerate(states):
            loads = np.zeros(size, dtype=int)
            loads[d] = 2
            loads[list(S)] = 1
            W[n, ordinary.index(S)] = 1
            for a, b in seams:
                word = loads[a] * 3 + loads[b]
                for outcome in np.flatnonzero(raw[:, word]):
                    assert outcome < 9
                    ports = loads.copy()
                    ports[a], ports[b] = divmod(int(outcome), 3)
                    state = (
                        np.flatnonzero(ports == 2).item(),
                        tuple(np.flatnonzero(ports == 1)),
                    )
                    H[index[state], n] += raw[outcome, word]
            H[n, n] += 4 * kappa * ((6 - adj.sum(axis=1)) @ (loads**2) - 36)
        np.testing.assert_array_equal(H, H.T)
        for sites in ((0,), (0, 2, 3), (1, 5, 6)):
            inside = np.isin(np.arange(size), sites)
            P = np.diag([int(inside[d]) for d, _ in states])
            comm = H @ P - P @ H
            direct = W.T @ (comm @ H - H @ comm) @ W
            a, b = operators(board, inside, normalize=False)
            # Before dividing by the number of empty ports these entries are
            # dyadic. The integer matrix comparison is exact for every column.
            rebuilt = (16 * (a + kappa * b)).toarray()
            np.testing.assert_array_equal(direct, rebuilt.astype(np.int64))
            np.testing.assert_array_equal(rebuilt, rebuilt.astype(np.int64))
            exact_refresh.append(int(np.count_nonzero(direct)))
    return exact_refresh


def instrument(operator, waves):
    upper = sparse.triu(operator, 1, format="csr")
    degree = (
        abs(operator.diagonal())
        + np.asarray((abs(upper) + abs(upper.T)).sum(axis=1)).ravel()
    )
    norm = max(degree)
    if norm == 0:
        return 0, np.tile([0.0, 0.0, 1.0], (waves.shape[1], 1))
    # Route each pair coherently into one arm per supported partner and
    # recombine the two paths for each unordered pair; the unused port is counted.
    active = degree @ (abs(waves) ** 2) / norm
    value = np.sum(waves.conj() * (operator @ waves), axis=0).real / norm
    probs = np.column_stack(((active + value) / 2, (active - value) / 2, 1 - active))
    assert min(probs.ravel()) >= -1e-12
    return float(norm), np.maximum(probs, 0)


def log_factor(source):
    a, b = contact.overlap_bands(source)["connected_charge_nats"]
    lower, upper = contact.overlap_bands(source)["reference_log_odds"]
    assert a > 0 and upper < 0
    return np.array([-a / lower, -b / upper])


def stats(bank, norm):
    pos, neg, miss = bank
    assert pos + neg + miss == TRIALS and min(bank) >= 0
    return norm * (
        (pos - neg) / TRIALS + np.array([-1, 1]) * np.sqrt(2 * contact.RISK / TRIALS)
    )


def mix_bands(x, y):
    z = np.multiply.outer(x, y)
    return np.array([np.min(z), np.max(z)])


def decode(records, source, gain):
    factor = log_factor(source)
    output = {}
    for t, regions in records.items():
        output[t] = {}
        for key, path in regions.items():
            quantities = {}
            for tag, (norm, banks) in path.items():
                vals = {label: stats(bank, norm) for label, bank in banks.items()}
                m = vals["00"] + vals["11"] - vals["10"][::-1] - vals["01"][::-1]
                quantities[tag] = dict(
                    from_charge=mix_bands(factor, m),
                    held=vals["product"] - vals["native"][::-1],
                )
            output[t][key] = {
                k: (
                    quantities["kinetic"][k]
                    + mix_bands(gain, quantities["contrast"][k])
                ).tolist()
                for k in ("from_charge", "held")
            }
    return output


def run(old, tape, saved=None):
    gain = motion.timing.lab.ordinary_gain(tape) / 4
    if saved is not None:
        return decode(saved["ordinary_banks"], old["source"], gain)
    cross = check_generator()
    edges, table, _, _ = contact.source_counts()
    nat = table.ravel() / table.sum()
    independent = np.outer(table.sum(1), table.sum(0)).ravel() / table.sum() ** 2
    board = motion.cold.federation(0)
    h, z, _, _, _ = motion.hamiltonian(board, edges)
    cases = list(product(motion.VALUES, repeat=2))
    psi = np.ones(h.shape[0]) / np.sqrt(h.shape[0])
    waves, cert = motion.chebyshev(h, z, cases, psi, motion.TIMES)
    rad = motion.cold.distances(board, sorted({n for e in edges for n in e}))
    rng = np.random.default_rng(686203)
    raw = {str(t): {} for t in motion.TIMES}
    native = {str(t): {} for t in motion.TIMES}
    for r in motion.RADII:
        k, v = operators(board, rad <= r)
        for sharp in (False, True):
            key = f"{r}_{'read' if sharp else 'undisturbed'}"
            for t, w in zip(motion.TIMES, waves, strict=True):
                raw[str(t)][key] = {}
                target = 0.0
                for tag, operator, coefficient in (
                    ("kinetic", k, 1),
                    ("contrast", v, motion.K),
                ):
                    mat = (
                        sparse.diags(operator.diagonal(), format="csr")
                        if sharp
                        else operator
                    )
                    norm, prob = instrument(mat, w)
                    conditional = motion.binary_responses(
                        {
                            ",".join(map(str, c)): p
                            for c, p in zip(cases, prob, strict=True)
                        }
                    )
                    probs = dict(
                        zip(("00", "01", "10", "11"), conditional, strict=True)
                    )
                    probs.update(
                        native=nat @ conditional, product=independent @ conditional
                    )
                    bank = {
                        label: rng.multinomial(BATCH, p / p.sum(), size=BATCHES)
                        .sum(0)
                        .tolist()
                        for label, p in probs.items()
                    }
                    raw[str(t)][key][tag] = [norm, bank]
                    diff = probs["product"] - probs["native"]
                    target += coefficient * norm * (diff[0] - diff[1])
                native[str(t)][key] = float(target)
    decoded = decode(raw, old["source"], gain)
    for t, panels in decoded.items():
        for key, forecasts in panels.items():
            for lo, hi in forecasts.values():
                assert lo <= native[t][key] <= hi
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "PAIR_PROBE.md",
                HERE / "pair_contact.json",
                HERE / "mediated_release.json",
            )
        },
        full_matrix_nonzero_counts=cross,
        clock=cert,
        ordinary_banks=raw,
        decoded=decoded,
        model_release=native,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "pair_probe.json"
    old = json.loads((HERE / "pair_contact.json").read_bytes())
    tape = json.loads((HERE / "mediated_release.json").read_bytes())
    if args.counts_only:
        out = json.loads(path.read_bytes())
        assert out["decoded"] == run(old, tape, out)
    else:
        out = run(old, tape)
        if args.write:
            path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
        else:
            motion.timing.lab.lab.cold.apparatus.check(
                json.loads(path.read_bytes()), out
            )
    print(out["decoded"]["64"])
