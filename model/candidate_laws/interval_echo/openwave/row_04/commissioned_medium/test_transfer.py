"""Check the exposure-37 counterexamples directly from the two frozen data files."""

from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent


def test_long_flight_from_raw_edges():
    record = json.loads((HERE / "results.json").read_text())
    donor = json.loads((HERE.parents[1] / "row_30_formed/results.json").read_text())
    board_path = HERE.parents[4] / "executable_model/carriers/FEDERATION_LEVEL0.npz"
    with np.load(board_path, allow_pickle=False) as board:
        edges = sorted(tuple(sorted(edge)) for edge in board["edges"])
    ports = 240
    pairs = np.array(list(combinations(range(ports), 2)))
    index = np.full((ports, ports), -1)
    a, b = pairs.T
    index[a, b] = index[b, a] = np.arange(len(pairs))
    reference = 6 - np.bincount(np.ravel(edges), minlength=ports)
    home = [e for e in edges if e[0] // 12 == e[1] // 12 == 13]
    assert [list(e) for e in home] + ["quiet"] == record["outcome_order"]
    launch = np.zeros((len(pairs), 2))
    launch[[index[156, 157], index[156, 166]], [0, 1]] = 1
    failure = Fraction(364, 729)
    signed = (
        np.array([-13, 14]) / 27 * np.sqrt(float(sum(failure**i for i in range(8))))
    )
    records = {s["bank"]: s for s in record["studies"]}
    for material in (s for s in donor["studies"] if s["cells"] == 20):
        study = records[material["seed"]]
        packed = bytes.fromhex(
            material["formation"]["coin_pairs_packed_little_bit_first_hex"]
        )
        # Decode each two-bit word by address and add both directions of each allowed move.
        rows, cols, live = [], [], set()
        for i, (p, q) in enumerate(edges):
            if (packed[i // 4] >> (2 * (i % 4))) & 3 == 3:
                continue
            live.add((p, q))
            other = np.delete(np.arange(ports), [p, q])
            rows.append(index[p, other])
            cols.append(index[q, other])
        x, y = np.concatenate(rows), np.concatenate(cols)
        moves = sparse.csr_matrix(
            (np.ones(2 * len(x)), (np.r_[x, y], np.r_[y, x])),
            shape=(len(pairs), len(pairs)),
        )
        count = np.asarray(moves.sum(axis=1)).ravel()
        h = (sparse.diags(count, dtype=float) - moves) / 2
        h += sparse.diags(0.5 * (count + reference[a] + reference[b]))
        states = expm_multiply(-37j * h, launch)
        np.testing.assert_allclose(
            states.conj().T @ states, np.eye(2), atol=2e-11, rtol=0
        )
        at_home = states[index[tuple(np.array(home).T)]]
        assert (156, 157) in live
        acquired = signed if (156, 166) in live else np.array([1.0, 0.0])
        (entry,) = [
            r
            for r in study["appointments"]
            if r["exposure"] == 37 and r["kappa"] == 0.5
        ]
        for arm, initial in (("filed_input", signed), ("acquired", acquired)):
            bins = abs(at_home @ initial) ** 2
            np.testing.assert_allclose(
                np.r_[bins, 1 - bins.sum()], entry["bins"][arm], atol=2e-10, rtol=0
            )
        best = np.linalg.svd(at_home, compute_uv=False)[0] ** 2
        np.testing.assert_allclose(
            best, entry["best_normalized_XZ_home"], atol=2e-10, rtol=0
        )
        assert best < 0.94
