"""Binary port-interference and head-analyzer records of the two terms in the shift Ward identity."""

import argparse
import json
from fractions import Fraction

import networkx as nx
import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import bridge
from . import measure


def selection_norm(op, unit):
    """Recover O from port-pair X or Y and two single-head Pauli outcomes."""
    eye = np.eye(2)
    pauli = (
        eye,
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]),
        np.diag([1, -1]),
    )
    basis = np.array([np.kron(a, b) for a in pauli for b in pauli])
    nz = op.tocoo()
    upper = sorted(
        {(i // 4, j // 4) for i, j in zip(nz.row, nz.col) if i // 4 <= j // 4}
    )
    cost = 0
    for i, j in upper:
        local = op[4 * i : 4 * i + 4, 4 * j : 4 * j + 4].toarray() / unit
        weights = np.einsum("kij,ji->k", basis, local) / 4
        np.testing.assert_array_equal(np.einsum("k,kij->ij", weights, basis), local)
        cost += np.sum(abs(weights.real) + abs(weights.imag))
    return float(cost)


def big_binomial(rng, n, p):
    block = 2**60
    assert n % block == 0
    counts = rng.binomial(block, p, size=n // block)
    return sum(map(int, counts))


def run():
    rng = np.random.default_rng(83192638)
    atlas = measure.Atlas(1)
    board = atlas.board
    patch = measure.chosen(atlas)
    route = bridge.native.iso.pairing.formation.region(board)
    graph = nx.Graph(board.seams)
    graph.remove_edges_from((route[:2], route[2:]))
    assert nx.is_connected(graph)
    pos = bridge.beam.indices(
        sum((list(r) for r in bridge.beam.screens(board, route)), [])
    )
    states = bridge.action.preparations()
    angles = {
        ("ordinary" if row["source"] == "ordinary" else "pairing"): row[
            "angle_per_unit"
        ]
        for row in json.loads(
            (bridge.HERE.parent / "signed_transfers.json").read_text()
        )["tests"]
    }
    n = 2**76
    rows = {}
    for q in (-1, 1):
        effects = {}
        norms = {}
        for odd in (False, True):
            h, ops = measure.op_table(board, route, patch, q, odd)
            start = np.zeros((h.shape[0], len(pos)))
            start[pos] = np.eye(len(pos))
            ht = -2j * h
            w = expm_multiply(ht, start, traceA=ht.diagonal().sum())
            for name in ("boundary", "marked_force", "free_rate"):
                op, unit = ops[name]
                cost = selection_norm(op, unit)
                norms[name] = max(cost, norms.get(name, 0))
                m = [a.conj().T @ (op @ a) / unit for a in np.hsplit(w, 10)]
                effects.setdefault(name, []).append(np.mean(m, axis=0))
        reads = {}
        for name, m in effects.items():
            scale = norms[name]
            expected = bridge.evaluate(m, states)
            cnt = {
                key: big_binomial(rng, n, 0.5 + mu / (2 * scale))
                for key, mu in expected.items()
            }
            est = {key: float(Fraction(2 * c - n, n)) * scale for key, c in cnt.items()}
            jets = {}
            for key, angle in angles.items():
                val = (
                    est[key]
                    - est["equilibrium"]
                    + angle * (est["plus"] - est["minus"]) / 2
                )
                width = (
                    scale
                    * np.sqrt(2 * np.log(2 * 2 * 3 * 5 / 0.01) / float(n))
                    * (2 + abs(angle))
                )
                ref = (
                    expected[key]
                    - expected["equilibrium"]
                    + angle * (expected["plus"] - expected["minus"]) / 2
                )
                assert abs(val - ref) < width and abs(val) > width
                jets[key] = dict(estimate=val, halfwidth=width)
            reads[name] = dict(plus_counts=cnt, range=scale, slopes=jets)
        rows[str(q)] = reads
    return dict(
        seed=83192638,
        simulated_trials_per_setting=n,
        simultaneous_coverage=0.99,
        signs=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "counts.json").write_text(data)
    else:
        print(data, end="")
