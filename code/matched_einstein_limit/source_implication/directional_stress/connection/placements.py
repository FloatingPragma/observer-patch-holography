"""Translate the mobile comparison and the two actual marked contacts together."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .. import bridge
from ..translation import measure as translation

HERE = Path(__file__).resolve().parent


def physical_twice(board, route, q, odd):
    cells = bridge.operators(board, route, q, odd)
    full = bridge.masked(cells, range(board.cells)) + bridge.bridges(
        board, board.external
    )
    shift = full - bridge.native.generator(board, route, q, odd)
    diff = int(shift.diagonal()[0])
    assert (shift - diff * sparse.eye(full.shape[0], dtype=int)).count_nonzero() == 0
    return full


def traces(h, f):
    values = []
    op = f.tocoo()
    cols = np.unique(op.col)
    v = np.zeros((h.shape[0], len(cols)), dtype=np.int64)
    v[cols] = np.eye(len(cols), dtype=np.int64)
    assert int(abs(h).sum(axis=1).max()) ** 7 < 2**62
    for i in range(8):
        values.append(
            sum(
                int(a) * int(b)
                for a, b in zip(op.data, v[op.row, np.searchsorted(cols, op.col)])
            )
        )
        if i < 7:
            v = h @ v
    return values


def address_conservation():
    from model.candidate_laws.interval_echo.openwave.photon_identity.native_emission import (
        measure as native,
    )

    h = native.four_h(8)
    blank = np.arange(h.shape[0]) < native.LAW.N_SRC
    assert not np.any(h[np.ix_(blank, ~blank)])
    assert not np.any(h[np.ix_(~blank, blank)])
    return {"primitive_columns": int(h.shape[0]), "blank_written_transition_entries": 0}


def run():
    atlas = translation.Atlas(1)
    board = atlas.board
    patch = translation.chosen(atlas)
    route = tuple(bridge.native.iso.pairing.formation.region(board))
    s = translation.map_matrix(board, patch)
    target = tuple(int(s[:, 4 * p].nonzero()[0][0] // 4) for p in route)
    variants = {
        "fixed": route,
        "left_moved": target[:2] + route[2:],
        "right_moved": route[:2] + target[2:],
        "both_moved": target,
    }
    domain = s.T @ s
    image = s @ s.T
    identity = sparse.eye(s.shape[0], dtype=np.int64, format="csr")
    rings = bridge.beam.screens(board, route)
    src = bridge.beam.indices(sum((list(z) for z in rings), []))
    start = np.zeros((s.shape[0], len(src)))
    start[src] = np.eye(len(src))
    shifted = s @ start
    np.testing.assert_array_equal(shifted.T @ shifted, np.eye(len(src)))
    states = bridge.tangents()
    rows = {}
    for q in (-1, 1):
        terms = {}
        sizes = []
        broken = []
        for odd in (False, True):
            a = physical_twice(board, route, q, odd)
            trial = {}
            centered = bridge.native.generator(board, route, q, odd)
            for label, marks in variants.items():
                b = physical_twice(board, marks, q, odd)
                residual = image @ b @ image @ s - s @ domain @ a @ domain
                trial[label] = residual
                assert (residual.count_nonzero() == 0) == (label == "both_moved")
            assert (
                trial["fixed"] - trial["left_moved"] - trial["right_moved"]
            ).count_nonzero() == 0
            sizes.append({k: int(v.multiply(v).sum()) for k, v in trial.items()})
            moment = {
                k: traces(centered, x + x.T)
                for k, x in trial.items()
                if k in ("left_moved", "right_moved")
            }
            assert moment["left_moved"][-1] + moment["right_moved"][-1] == 64
            broken.append(moment)
            b = physical_twice(board, target, q, odd)
            cut = (identity - image) @ b @ s - s @ a @ (identity - domain)
            assert (b @ s - s @ a - cut).count_nonzero() == 0

            def propagate(h, state, t):
                ht = -0.5j * t * h
                return expm_multiply(ht, state, traceA=ht.diagonal().sum())

            wa = propagate(a, start, 4)
            wb = propagate(b, shifted, 4)

            def cross(wa, wb, mat):
                shifted = mat @ wa
                h = (
                    sum(
                        y.conj().T @ x
                        for x, y in zip(np.hsplit(shifted, 10), np.hsplit(wb, 10))
                    )
                    / 10
                )
                return (h + h.conj().T) / 2, (h - h.conj().T) / (2j)

            re, im = cross(wa, wb, s)
            cre, cim = cross(wa, wb, cut)
            step = 1 / 2048
            rates = []
            for sign in (-1, 1):
                aa = propagate(a, wa, sign * step)
                bb = propagate(b, wb, sign * step)
                rates.append(cross(aa, bb, s))
            # 2<Q_AB> is Im Z, and its quadrature is Re Z.
            for k, m in (
                ("real_overlap", re),
                ("imag_overlap", im),
                ("imag_rate_from_boundary", cre / 2),
                ("real_rate_from_boundary", -cim / 2),
                (
                    "imag_rate_from_free_flights",
                    (rates[1][1] - rates[0][1]) / (2 * step),
                ),
                (
                    "real_rate_from_free_flights",
                    (rates[1][0] - rates[0][0]) / (2 * step),
                ),
            ):
                terms.setdefault(k, []).append(m)
        values = {k: bridge.evaluate(m, states) for k, m in terms.items()}
        for quad in ("real", "imag"):
            for name in states:
                assert (
                    abs(
                        values[quad + "_rate_from_free_flights"][name]
                        - values[quad + "_rate_from_boundary"][name]
                    )
                    < 0.0001
                )
        rows[str(q)] = dict(
            squared_hilbert_schmidt_residual_by_parity=sizes,
            independent_head_residual_moments=broken,
            coherent_reference=values,
        )
    return dict(
        actual_addresses=route,
        translated_addresses=target,
        free_address_grammar=address_conservation(),
        diagonal_address_ensemble_interference_exact=0,
        coherent_comparison_readout="Z=Tr(U_B^dagger S U_A rho_AB) with coincident cross-branch amplitudes normalized to one at t=0",
        charges=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "placements.json").write_text(data)
    else:
        print(data, end="")
