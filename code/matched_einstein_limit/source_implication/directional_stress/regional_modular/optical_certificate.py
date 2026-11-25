"""Enclose exchange-clock focusing on the exact fixed-radar equal-traction tangent."""

import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv

from ...jacobi_screen import certificate as optical
from ...jacobi_screen import measure as jacobi
from . import measure

c = measure.c


def hull(values):
    return iv.mpf([min(x.a for x in values), max(x.b for x in values)])


def probe(operators, ref, dirs, source_error, pooled_charges=1):
    states, ref_radius, odds = ref
    p = 0.1343
    nominal = np.array(
        [np.array(b, dtype=complex) * w for b, w in zip(states, (1 - p, p))]
    )
    denominator = pooled_charges * optical.prior.FRAC**2
    base = (
        np.array(operators[:, 0, 0], dtype=float) / denominator
        + 1j * np.array(operators[:, 0, 1], dtype=float) / denominator
    )
    out, inc = optical.translated_basis(optical.beam.probabilities(base, nominal))
    baseline, variation = [], []
    for k in range(3):
        denom = 2**k * denominator * optical.beam.UNIT**2
        sector = []
        jets = []
        for a in out.T:
            for b in inc.T:
                weight = np.outer(a, b)
                assert sum(abs(v) for v in weight.ravel()) <= 100 * optical.beam.UNIT**2
                re, im = (
                    (operators[:, k] * weight[None, None, :, :, None, None])
                    .sum((2, 3))
                    .swapaxes(0, 1)
                )
                blocks = [
                    (sp.Matrix(x) + sp.I * sp.Matrix(y)) / denom for x, y in zip(re, im)
                ]
                err = sp.Rational(100 * 10 * 64**k, 10**19)
                centers = [c.ri(c.inner(v, s)) for v, s in zip(blocks, states)]
                effective = [v - v.trace() * sp.eye(4) / 4 for v in blocks]
                norms = [iv.sqrt(c.ri(c.inner(v, v))) for v in effective]
                er = max(n.b for n in norms) * c.ri(ref_radius) + c.ri(err)
                sector.append((centers, er))
                slope = c.ri(sum(c.inner(v, d) for v, d in zip(effective, dirs)))
                # Tr norm of D is at most sqrt(8)||D||_HS <3 on this source.
                wl = iv.sqrt(sum(v * v for v in norms)) * c.ri(source_error) + 3 * c.ri(
                    err
                )
                jets.append(slope + iv.mpf([-1, 1]) * wl)
        baseline.append(sector)
        variation.append(np.array(jets, dtype=object).reshape(3, 3))
    vals = []
    for index in range(8):
        lo = odds.a + (odds.b - odds.a) * index / 8
        hi = odds.a + (odds.b - odds.a) * (index + 1) / 8
        pp = iv.mpf([lo.a, hi.b])
        m = [
            np.array(
                [e + pp * (o - e) + iv.mpf([-1, 1]) * er for (e, o), er in row],
                dtype=object,
            ).reshape(3, 3)
            for row in baseline
        ]
        j = jacobi.matrix_jets(m)
        det, b, cc = jacobi.optical_matrices(j)
        assert 0 not in det
        db, dc = jacobi.deformation_response(j, jacobi.response_jets(m, variation))
        vals.append([b.trace(), -cc.trace(), db.trace(), -dc.trace()])
    return [hull(group) for group in zip(*vals)]


def common_clock(panels):
    from ...action_response.four_ticks import parse_interval

    result = {}
    for tag in ("forward", "reverse"):
        a, b = [panels[q][tag] for q in ("-1", "1")]
        focusing = parse_interval(a["source_focusing_slope"]) - parse_interval(
            b["source_focusing_slope"]
        )
        expansion = parse_interval(a["source_expansion_slope"]) - parse_interval(
            b["source_expansion_slope"]
        )
        assert 0 not in expansion
        result[tag] = str(-focusing / expansion)
    assert parse_interval(result["forward"]).b < parse_interval(result["reverse"]).a
    return result


def run(arm=None, sample=None, *, compare_charges=True):
    iv.dps = 32
    board = measure.bridge.native.iso.pairing.formation.federation(1)
    route = measure.bridge.native.iso.pairing.formation.region(board)
    rings = measure.bridge.beam.screens(board, route)
    src = measure.source()
    if arm is None:
        arm = json.loads((measure.HERE / "control.json").read_text())
    coords = np.array(arm["center_coordinates"]).reshape(2, 15)
    frame = [c.prior.sym_matrix(x) for x in c.iso.basis(False)]
    d = [
        sum((sp.Rational(v) * op for v, op in zip(row, frame)), sp.zeros(4))
        for row in coords
    ]
    err = sp.Rational(arm["direction_error_with_serialization_upper"])
    weights = measure.bridge.native.prior.weights()[[7]]
    if sample is None:
        sample = json.loads((measure.HERE / "screen_ward.json").read_text())
    out = {}
    banks = []
    for q in (-1, 1):
        tables = {"forward": [], "reverse": []}
        for odd in (False, True):
            h = measure.bridge.native.generator(board, route, q, odd)
            amplitudes = np.array(
                optical.amplitude_jets(h, *rings, weights), dtype=object
            )
            for name, w in (
                ("forward", amplitudes),
                ("reverse", amplitudes.transpose(0, 1, 2, 4, 3, 6, 5)),
            ):
                tables[name].append(optical.probability_jets(w)[:, :, 0])
        banks.append(tables)
        if not compare_charges:
            continue
        out[str(q)] = {}
        for name, row in tables.items():
            z = probe(np.array(row, dtype=object), src, d, err)
            direct = sample["charges"][str(q)]["directional_images"][name][
                "material_terminal_ward_read"
            ]
            assert all(v in x for v, x in zip(direct, z))
            out[str(q)][name] = dict(zip(sample["order"], map(str, z)))
    pooled = {}
    for name in ("forward", "reverse"):
        total = np.array(banks[0][name], dtype=object) + np.array(
            banks[1][name], dtype=object
        )
        z = probe(total, src, d, err, pooled_charges=2)
        direct = sample["fair_untagged_images"][name]["material_terminal_ward_read"]
        assert all(v in x for v, x in zip(direct, z))
        pooled[name] = dict(zip(sample["order"], map(str, z)))
    return dict(
        flight_time=4,
        weights=sample["launches"],
        fixed_radar=True,
        work_derivative=arm["work_derivative_exact"],
        charges=out,
        fair_untagged_images=pooled,
        required_common_nonaffinity_for_equal_work=common_clock(out)
        if compare_charges
        else None,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "optical_certificate.json").write_text(data)
    else:
        print(data, end="")
