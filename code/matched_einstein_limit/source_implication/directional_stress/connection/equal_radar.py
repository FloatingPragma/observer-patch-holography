"""Full marked-force response along the exact twelve-range-invisible source chord."""

import argparse
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from ...full_algebra import certificate as cert
from ...full_algebra.reuse import invisible
from .. import bridge
from ..traction.certificate import four_tick_columns
from ..translation import certificate as material
from ..translation.counts import big_binomial, selection_norm
from . import late_read, measure

HERE = measure.HERE
TRANSFER_NORM_ERROR = sp.Rational("2e-15")


def tangent_source():
    path = Path(invisible.__file__).with_suffix(".json")
    source = json.loads(path.read_text())
    assert source["range_rows"] == 12 and source["exact_range_kernel_dimension"] == 6
    # The producer serializes each exactly corrected rational coordinate as float.
    # Round-trip binary64 coordinates below 2 incur <2^-52 absolute error each.
    coords = sp.Matrix([sp.Rational(x) for x in source["tangent_center_even_then_odd"]])
    assert len(coords) == 18 and max(abs(x) for x in coords) < 2
    error = sp.Rational(source["true_kernel_correction_HS_upper"]) + sp.Rational(
        5, 2**52
    )
    assert error < sp.Rational("4.1e-7")
    basis = [cert.prior.sym_matrix(h) for h in invisible.station.iso.basis(True)]
    assert sp.Matrix([[cert.inner(a, b) for b in basis] for a in basis]) == sp.eye(9)
    heads = [
        sum((v * h for v, h in zip(coords[i : i + 9], basis)), sp.zeros(4))
        for i in (0, 9)
    ]
    return source, coords, error, basis, heads


def delayed_read(board, route, patch, launches, charge, states, coefficient):
    """Read two positive sources through calibrated late local hits and both phase quadratures."""
    mp.mp.dps = 52
    s = measure.translation.map_matrix(board, patch)
    images = [int(s[:, 4 * p].nonzero()[0][0]) // 4 for p in route]
    ports = sorted(set(launches + list(route) + images))
    pairs = [
        (p, t)
        for i, p in enumerate(route)
        for t in images[2 * (i // 2) : 2 * (i // 2) + 2]
    ]
    paulis = [
        np.eye(2),
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]),
        np.diag([1, -1]),
    ]
    analyzers = [late_read.mm(np.kron(a, b)) for a in paulis for b in paulis]
    sums = [mp.mpf(0), mp.mpf(0)]
    detected = {}
    for odd in (False, True):
        h, ops = measure.translation.op_table(board, route, patch, charge, odd)
        numerator, unit = ops["marked_force"]
        ids = bridge.beam.indices(ports)
        wave = four_tick_columns(h, ids, ids, coefficient)

        def block(p, t):
            i, j = ports.index(p), ports.index(t)
            return late_read.mp_wave(*wave[:, 4 * i : 4 * i + 4, 4 * j : 4 * j + 4])

        records = {}
        for p, t in pairs:
            a, b = block(p, p), block(p, t)
            ca, cb, _ = late_read.calibrate(a, b, analyzers, False)
            force = (
                late_read.mm(numerator[4 * t : 4 * t + 4, 4 * p : 4 * p + 4].toarray())
                / unit
            )
            dual = (cb.H**-1) * force * (ca**-1)
            early_p = [block(p, k) for k in launches]
            early_t = [block(t, k) for k in launches]
            pair = []
            for i, group in enumerate(states):
                read, _ = late_read.decode_local(
                    a,
                    b,
                    early_p,
                    early_t,
                    late_read.mm(group[int(odd)]),
                    dual,
                    analyzers,
                    False,
                )
                original = read / 2
                sums[i] += original
                pair.append(mp.nstr(original, 25))
            records[f"{p}_{t}"] = pair
        detected[str(int(odd))] = records
    return sums, detected


def run():
    cert.iv.dps = 32
    source, coords, error, corr, heads = tangent_source()
    atlas = measure.translation.Atlas(1)
    board = atlas.board
    route = bridge.native.iso.pairing.formation.region(board)
    launches = sum(bridge.beam.screens(board, route), [])
    patch = measure.translation.chosen(atlas)
    coefficient = bridge.native.prior.weights()[7]
    norm = cert.iv.sqrt(cert.ri(coords.dot(coords)))
    step = sp.Rational(source["physical_step"])
    # Native calibration states around the two published conditional fixed-point
    # centers. The exact force difference is independent of the reference mixture.
    trial = sp.Rational(1, 100)
    p = sp.Rational(".1343")
    frame = [cert.prior.sym_matrix(h) for h in invisible.station.iso.basis(False)]
    info = source["fixed_state_centers"]
    conditional = [
        sp.eye(4) / 4
        + sum(
            (sp.Rational(n, info["denominator"]) * a for n, a in zip(row, frame)),
            sp.zeros(4),
        )
        for row in info["even_then_odd"]
    ]
    states = [
        [
            weight * center + sign * trial * block
            for weight, center, block in zip((1 - p, p), conditional, heads)
        ]
        for sign in (-1, 1)
    ]
    for state in states:
        assert all(cert.global_survey.positive_definite(m) for m in state)
        # The two admissible native Pauli lotteries have positive identity slack.
        for block, weight in ((state[0] - state[1], 1 - 2 * p), (2 * state[1], 2 * p)):
            assert sum(abs(2 * cert.inner(block, op)) for op in frame) < weight
    rows = {}
    rng = np.random.default_rng(28112026)
    trials = 2**82
    for q in (-1, 1):
        blocks = []
        scales = []
        for odd in (False, True):
            h, ops = measure.translation.op_table(board, route, patch, q, odd)
            force, unit = ops["marked_force"]
            scales.append(selection_norm(force, unit))
            assert unit == 4 and abs(force).sum(axis=1).max() / unit < 200
            idx = np.unique(force.nonzero()[0])
            w = four_tick_columns(h, idx, bridge.beam.indices(launches), coefficient)
            blocks.append(material.compress(w, force[idx, :][:, idx], unit))
        read = sp.Matrix([cert.inner(op, h) for op in blocks for h in corr])
        # Per-block operator norm error <2e-15 by the regional propagation
        # certificate, hence <=sqrt(8)*2e-15 HS over both four-dimensional blocks.
        sensitivity = cert.iv.sqrt(cert.ri(read.dot(read)))
        radius = sensitivity * cert.ri(error) + 3 * cert.ri(TRANSFER_NORM_ERROR) * (
            norm + cert.ri(error)
        )
        slope = cert.ri(read.dot(coords)) + cert.iv.mpf([-1, 1]) * radius
        centers = [
            cert.ri(cert.inner(op, sigma)) for op, sigma in zip(blocks, conditional)
        ]
        weight = cert.iv.mpf(source["matched_odd_weight"])
        radial = []
        for op in blocks:
            centered = op - op.trace() * sp.eye(4) / 4
            radial.append(
                cert.iv.sqrt(cert.ri(cert.inner(centered, centered)))
                * cert.ri(sp.Rational(source["conditional_HS_error_upper"]))
                + cert.ri(TRANSFER_NORM_ERROR)
            )
        baseline = (
            centers[0]
            + weight * (centers[1] - centers[0])
            + cert.iv.mpf([-1, 1]) * max(r.b for r in radial)
        )
        values, detections = delayed_read(
            board, route, patch, launches, q, states, coefficient
        )
        sampled = (values[1] - values[0]) / (2 * mp.mpf(str(float(trial))))
        assert float(sampled) in slope and 0 not in slope
        nominal = [
            float(sum(cert.inner(op, s) for op, s in zip(blocks, arm)))
            for arm in states
        ]
        assert max(abs(c - float(v)) for c, v in zip(nominal, values)) < 2e-15
        scale = max(scales)
        counts = [big_binomial(rng, trials, 0.5 + y / (2 * scale)) for y in nominal]
        observed = 2 * scale * float(Fraction(counts[1] - counts[0], trials))
        # Four arm means simultaneously. The numerical probability conversion
        # error at this scale adds less than 2e-13 to their signed difference.
        halfwidth = 2 * scale * np.sqrt(2 * np.log(8 / 0.01) / float(trials)) + 2e-13
        assert (
            observed > halfwidth
            and abs(observed - (nominal[1] - nominal[0])) < halfwidth
        )
        assert 2 * halfwidth < float((2 * cert.ri(trial) * slope).a)
        rows[str(q)] = dict(
            force_derivative=str(slope),
            stationary_reference_force=str(baseline),
            force_change_on_certified_chord=str(cert.ri(step) * slope),
            force_difference_between_native_endpoints=str(2 * cert.ri(trial) * slope),
            nominal_late_local_derivative=mp.nstr(sampled, 30),
            total_projected_force_norm=str(sensitivity),
            calibration_arm_force=[mp.nstr(v, 30) for v in values],
            direct_native_counts=dict(
                scale=scale,
                plus=counts,
                trials_per_arm=trials,
                difference=observed,
                simultaneous_99pct_halfwidth=halfwidth,
            ),
            pair_contributions=detections,
        )
    return dict(
        native_ticks=4,
        launches=launches,
        stationary_weight=source["matched_odd_weight"],
        certified_step=source["physical_step"],
        source_coordinate_count=18,
        source_correction_upper=float(np.nextafter(float(error), np.inf)),
        range_rows_identical=12,
        entropy_derivative_nats=source["entropy_derivative_nats"],
        nominal_tomography_calibration_odd_weight=str(p),
        nominal_tomography_half_step=str(trial),
        charges=rows,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "equal_radar.json").write_text(data)
    else:
        print(data, end="")
