"""Fixed-point four-tick propagation and source enclosure for the tagged flux."""

import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv

from ..action_response.four_ticks import parse_interval, raw_effect
from . import bridge


def run():
    iv.dps = 30
    board = bridge.native.iso.pairing.formation.federation(1)
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    nodes = list(route[1:3]) + list(rings[0]) + list(rings[1])
    idx = bridge.beam.indices(nodes)
    coefficients = bridge.native.prior.weights()
    source = bridge.action.front.inputs()
    ratio = json.loads(
        (bridge.HERE.parent / "action_response/four_ticks.json").read_text()
    )
    volume_ratio = parse_interval(ratio["ratio_of_volume_adjustment_angles"])
    rows = {}
    for q in (-1, 1):
        effects = []
        for odd in (False, True):
            h = bridge.native.generator(board, route, q, odd)
            cells = bridge.operators(board, route, q, odd)
            k = (
                bridge.current(board, cells, route[1:3])[idx, :][:, idx]
                .toarray()
                .astype(object)
            )
            assert np.max(np.abs(k).sum(axis=1)) <= 8 * 500
            raw = bridge.native.kernel_panels(h, idx, coefficients)[:, 7]
            real = sp.zeros(4)
            imag = sp.zeros(4)
            for launch in range(2, len(nodes)):
                a, b = raw[:, :, launch].reshape(2, len(idx), 4)
                real += sp.Matrix(b.T @ k @ a - a.T @ k @ b)
                imag += sp.Matrix(a.T @ k @ a + b.T @ k @ b)
            effect = (real + sp.I * imag) / (10 * 8 * bridge.native.prior.FRAC**2)
            assert effect.H == effect
            effects.append(effect)
        # 12 four-row blocks, error <= 1e-19/block. Compressed U has norm <= 1,
        # ||J|| <= 500, so the effect error < 4e-16 and the jet error < 4e-13.
        allowance = iv.mpf(["-1e-12", "1e-12"])
        jets = {
            k: v + allowance
            for k, v in bridge.action.source_jets(effects, source).items()
        }
        raw_o = raw_effect(
            effects, source[1]["ordinary"], sp.Rational("12e-7"), sp.Rational("1e-15")
        )
        raw_p = raw_effect(
            effects, source[1]["unpair"], sp.Rational("18e-7"), sp.Rational("1e-15")
        )
        reference = raw_effect(
            effects, source[0], sp.Rational("6e-7"), sp.Rational("1e-15")
        )
        direct = jets["ordinary"] / jets["pairing"]
        shared = volume_ratio + (raw_o - volume_ratio * raw_p) / jets["pairing"]
        tighter = iv.mpf([max(direct.a, shared.a), min(direct.b, shared.b)])
        nominal = json.loads((bridge.HERE / "bridge.json").read_text())["charges"][
            str(q)
        ]["equal_tagged_launch"]
        for name, interval in jets.items():
            assert nominal[name] in interval and interval.b < 0
        assert nominal["ordinary"] / nominal["pairing"] in tighter
        rows[str(q)] = dict(
            reference=str(reference),
            ordinary=str(jets["ordinary"]),
            pairing=str(jets["pairing"]),
            ordinary_over_pairing=str(tighter),
            pairing_over_ordinary=str(1 / tighter),
        )
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (bridge.HERE / "certificate.json").write_text(output)
    else:
        print(output, end="")
