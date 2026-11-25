"""Exact compression of the integral translation operators on the stationary four-tick flight."""

import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv

from ...action_response.four_ticks import raw_effect
from .. import bridge
from ..traction.certificate import four_tick_columns
from . import measure


def compress(w, mat, unit):
    a, b = w
    co = mat.tocoo()
    ar, ai = np.zeros_like(a), np.zeros_like(b)
    np.add.at(ar, co.row, co.data[:, None] * a[co.col])
    np.add.at(ai, co.row, co.data[:, None] * b[co.col])
    real = sp.zeros(4)
    imag = sp.zeros(4)
    for x, y, hx, hy in zip(*(np.hsplit(item, 10) for item in (a, b, ar, ai))):
        real += sp.Matrix(x.T @ hx + y.T @ hy)
        imag += sp.Matrix(x.T @ hy - y.T @ hx)
    divisor = 2 * sp.I if unit == 2j else unit
    m = (real + sp.I * imag) / (10 * divisor * bridge.native.prior.FRAC**2)
    m = m.applyfunc(sp.expand)
    assert m.H == m
    return m


def run():
    iv.dps = 30
    atlas = measure.Atlas(1)
    board = atlas.board
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    pos = bridge.beam.indices(list(rings[0]) + list(rings[1]))
    source = bridge.action.front.inputs()
    co = bridge.native.prior.weights()[7]
    expected = json.loads((measure.HERE / "results.json").read_text())["signs"]
    data = {}
    for q in (-1, 1):
        materials = {}
        for odd in (False, True):
            h, ops = measure.op_table(board, route, measure.chosen(atlas), q, odd)
            idx = np.unique(np.concatenate([m.nonzero()[0] for m, u in ops.values()]))
            w = four_tick_columns(h, idx, pos, co)
            for k, (m, u) in ops.items():
                assert abs(m).sum(axis=1).max() / abs(u) < 200
                materials.setdefault(k, []).append(compress(w, m[idx, :][:, idx], u))
        row = {}
        for k, x in materials.items():
            # <=960 compressed port blocks, each with operator error <1e-19.
            # ||O||<200 gives <1.3e-15 per expectation and <1.3e-12 per matched jet.
            val = bridge.action.source_jets(x, source)
            val = {v: t + iv.mpf(["-2e-12", "2e-12"]) for v, t in val.items()}
            val["reference"] = raw_effect(
                x, source[0], sp.Rational("6e-7"), sp.Rational("2e-15")
            )
            assert all(expected[str(q)][k][v] in t for v, t in val.items())
            row[k] = {v: str(t) for v, t in val.items()}
            if k in ("marked_force", "boundary", "free_rate"):
                assert all(0 not in t for t in val.values())
        data[str(q)] = row
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "certificate.json").write_text(data)
    else:
        print(data, end="")
