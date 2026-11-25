"""Bounded fixed-point force and work responses at the stationary source."""

import argparse
import json
from fractions import Fraction
from math import isqrt

import numpy as np
import sympy as sp
from mpmath import iv

from ...action_response.four_ticks import parse_interval, raw_effect
from .. import bridge
from . import measure, strain


def four_tick_columns(h, ends, starts, coefficient):
    """Selected rows, columns and time of native_precision.kernel_panels."""
    prior = bridge.native.prior
    degree = prior.DEGREE
    row = int(abs(h).sum(axis=1).max())
    eps = Fraction(
        3 * (1 + isqrt(h.shape[0])) * degree * (degree + 1) * (degree + 2),
        6 * prior.UNIT,
    )
    eps += Fraction(3 * (degree + 1), prior.FRAC) + Fraction(1, 10**59)
    assert 4 * eps < Fraction(1, 10**19) and row < 64
    high = np.zeros((h.shape[0], len(starts)), dtype=np.int64)
    low = high.copy()
    high[starts] = (prior.UNIT // prior.BASE) * np.eye(len(starts), dtype=np.int64)
    out = np.zeros((2, len(ends), len(starts)), dtype=object)
    prev_hi = prev_lo = None
    for k in range(degree + 1):
        assert row * (int(np.abs(high).max()) + 1) < 2**62 and row * prior.BASE < 2**62
        val = high[ends].astype(object) * prior.BASE + low[ends]
        out[k % 2] += (1 if k % 4 in (0, 3) else -1) * coefficient[k] * val
        if k < degree:
            nh, nl = prior.divided_product(h, high, low, 6 if k == 0 else 5)
            if k:
                nl -= prev_lo
                nh += (nl >> prior.LIMB) - prev_hi
                nl &= prior.MASK
            prev_hi, prev_lo, high, low = high, low, nh, nl
    return out // prior.UNIT


def compress(w, matrix, unit):
    a, b = w
    nz = np.flatnonzero(np.any(matrix, axis=1))
    h = matrix[np.ix_(nz, nz)].astype(object)
    a, b = a[nz], b[nz]
    real, imag = sp.zeros(4), sp.zeros(4)
    for x, y in zip(np.hsplit(a, 10), np.hsplit(b, 10)):
        hx, hy = h @ x, h @ y
        real += sp.Matrix(x.T @ hx + y.T @ hy)
        imag += sp.Matrix(x.T @ hy - y.T @ hx)
    effect = (real + sp.I * imag) / (10 * unit * bridge.native.prior.FRAC**2)
    assert effect.H == effect
    return effect


def round_rulers(route, charts, source):
    out = {}
    for slot, (name, chart) in enumerate(charts.items()):
        half = []
        for odd in (False, True):
            t = measure.marked_block(route, slot, 1, odd)[1, 2]
            t = (
                np.kron(t, np.eye(2, dtype=int))
                if slot == 0
                else np.kron(np.eye(2, dtype=int), t)
            )
            half.append(sp.Matrix(t @ t.T + t.T @ t) / 2)
        coeff = raw_effect(half, source[0], sp.Rational("6e-7"), sp.S.Zero)
        x = measure.UNIT / iv.sqrt(coeff)
        val = chart["marked_length_numerator"]
        assert x.a > val - 0.5 and x.b < val + 0.5
        out[name] = str(x / measure.UNIT)
    return out


def run():
    iv.dps = 30
    board = bridge.native.iso.pairing.formation.federation(1)
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    starts = bridge.beam.indices(list(rings[0]) + list(rings[1]))
    charts = measure.coordinates(board, route)
    coefficient = bridge.native.prior.weights()[7]
    source = bridge.action.front.inputs()
    ratio = json.loads(
        (bridge.HERE.parent / "action_response/four_ticks.json").read_text()
    )
    rho = parse_interval(ratio["ratio_of_volume_adjustment_angles"])
    rows = {}
    for q in (-1, 1):
        materials = {}
        for odd in (False, True):
            h = bridge.native.generator(board, route, q, odd)
            ports, families = measure.setup(board, route, q, odd, charts)
            w = four_tick_columns(h, bridge.beam.indices(ports), starts, coefficient)
            ops = {
                "_".join(key): (items[0], 16 * measure.UNIT)
                for key, items in families.items()
            }
            ops.update(
                {k: (rec, 2) for k, rec in strain.observables(ports, route).items()}
            )
            for key, (tensor, unit) in ops.items():
                assert np.max(np.sum(abs(tensor), axis=1)) < 20000 * unit
                materials.setdefault(key, []).append(compress(w, tensor, unit))
        reads = {}
        raw = {}
        nom = json.loads((measure.HERE / "results.json").read_text())["force"][str(q)]
        work = json.loads((measure.HERE / "strain.json").read_text())["charges"][
            str(q)
        ]["means"]
        for name, hs in materials.items():
            # 24 output-port blocks have operator error <1e-19 each, hence the
            # full rectangular wave error <5e-19. F norm<20000 gives <2.1e-14
            # bias per expectation and <2.1e-11 per matched jet.
            jets = bridge.action.source_jets(hs, source)
            raw[name] = {k: v + iv.mpf(["-3e-11", "3e-11"]) for k, v in jets.items()}
            reads[name] = {k: str(v) for k, v in raw[name].items()}
            val = raw_effect(hs, source[0], sp.Rational("6e-7"), sp.Rational("3e-14"))
            reads[name]["reference"] = str(val)
            nominal = nom[name]["commutator"] if name in nom else work[name]
            for key, interval in reads[name].items():
                assert nominal[key] in parse_interval(interval)
            raw_o = raw_effect(
                hs, source[1]["ordinary"], sp.Rational("12e-7"), sp.Rational("6e-14")
            )
            raw_p = raw_effect(
                hs, source[1]["unpair"], sp.Rational("18e-7"), sp.Rational("6e-14")
            )
            raw[name]["without_rotation"] = raw_o - rho * raw_p
        a, b = raw["full_seam"], raw["exchange_only"]
        determinant = (
            a["without_rotation"] * b["pairing"] - a["pairing"] * b["without_rotation"]
        )
        assert 0 not in determinant
        for k in ("ordinary", "pairing"):
            assert (raw["right_front"][k] * raw["right_hops"][k]).b < 0
        held = (
            1 + determinant / (a["pairing"] * b["ordinary"])
            if 0 not in b["ordinary"]
            else None
        )
        rows[str(q)] = dict(
            observables=reads,
            work_determinant=str(determinant),
            whole_seam_ordinary_over_pairing=str(
                rho + a["without_rotation"] / a["pairing"]
            ),
            ordinary_calibrated_work_ratio=str(held) if held is not None else None,
        )
    return dict(ruler_length=round_rulers(route, charts, source), charges=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "certificate.json").write_text(data)
    else:
        print(data, end="")
