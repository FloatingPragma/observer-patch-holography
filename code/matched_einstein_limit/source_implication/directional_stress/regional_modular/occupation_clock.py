"""Test running material clocks on the same source as the two null screens."""

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv

from ...action_response.four_ticks import parse_interval
from ...full_algebra.reuse import invisible
from . import control, energy_odd, measure
from .clock_area import measure as optics

c = measure.c
HERE = Path(__file__).parent


def occupation_rates():
    board = measure.bridge.native.iso.pairing.formation.federation(1)
    route = measure.bridge.native.iso.pairing.formation.region(board)
    rings = measure.bridge.beam.screens(board, route)
    weights = measure.bridge.native.prior.weights()[7]
    frame = [c.prior.sym_matrix(a) for a in c.iso.basis(False)]
    ref = measure.source()
    rows, probabilities, rates = [], [], []
    for charge in (-1, 1):
        panel, _ = measure.local_effects(board, route, rings, charge, weights)
        d = panel["reference"][0]
        dot = panel["total_rate"][0]
        rows.append([c.inner(op, b) for op in dot for b in frame])
        probabilities.append(measure.mean(d, ref, sp.Rational("3e-15")))
        rates.append(measure.mean(dot, ref, sp.Rational("3e-15")))
    return sp.Matrix(rows), probabilities, rates


def source_slope(row, arm):
    d = sp.Matrix([sp.Rational(float(v)) for v in arm["center_coordinates"]])
    radius = invisible.frobenius(row) * c.ri(
        arm["direction_error_with_serialization_upper"]
    ) + c.ri("1e-14")
    return c.ri((row * d)[0]) + iv.mpf([-1, 1]) * radius


def specified_clock(row, probabilities, rates, receipt):
    d = sum(probabilities) / 2
    dot = sum(rates) / 2
    assert d.a > 0
    results = {}
    for key in ("fixed_work_control", "fixed_work_second_control"):
        arm = receipt[key]
        optics = arm["screen"]["odd_forward_minus_reverse"]
        g = parse_interval(optics["source_focusing_slope"])
        theta = parse_interval(optics["expansion"])
        h = parse_interval(optics["source_expansion_slope"])
        delta_dot = source_slope((row.row(0) + row.row(1)) / 2, arm)
        # lambda(t) is accumulated regional occupation. delta d=0 at the appointment,
        # while delta dot(d) can change the acceleration of the material clock.
        affine_slope = (g + (dot * h + theta * delta_dot) / d) / d**2
        assert 0 not in affine_slope
        results[key] = dict(
            occupation_time_rate_derivative=str(delta_dot),
            odd_focusing_derivative_in_dwell_time=str(affine_slope),
        )
    return dict(
        clock_speed=str(d), clock_acceleration=str(dot), source_controls=results
    )


def frozen_rate_control(square, err, slot, component):
    target = sp.zeros(30, 1)
    target[26 + slot] = 1
    center, scale, correction, defect, norm = control.solve(square, err, target)
    assert (c.ri(".20") - 16 * c.ri(".01") * norm).a > 0
    arm = dict(
        work_derivative_exact="0",
        center_coordinates=[float(x) for x in center],
        direction_error_with_serialization_upper=float(
            np.nextafter(float(correction.b) + 2e-14, np.inf)
        ),
    )
    return dict(
        **arm,
        target_head_coordinate=component,
        target_derivative_exact=f"1/{scale}",
        augmented_right_inverse_defect_upper=str(defect.b),
    )


def run():
    iv.dps = 32
    rate, probabilities, rates = occupation_rates()
    receipt = json.loads((HERE / "energy_odd.json").read_text())
    dwell = specified_clock(rate, probabilities, rates, receipt)
    rows, err = control.constraint_rows()
    currents, _ = energy_odd.energy_rows()
    rows = rows.col_join(currents).col_join(rate)
    error = iv.sqrt(
        err**2 + c.ri(16 * (sp.Rational("1e-15") ** 2 + sp.Rational("3e-15") ** 2))
    )
    assert rows.shape == (26, 30)
    # A pivoted QR of the material rows selected these four fixed head axes
    # before computing the optics. The inverse defect proves full rank
    # for the actual instrument, independently of that numerical selection.
    indices = (9, 4, 20, 8)
    square = rows.col_join(sp.eye(30)[indices, :])
    arms = [
        frozen_rate_control(square, error, i, axis) for i, axis in enumerate(indices)
    ]
    raw = optics.image_records()
    ref = measure.source()
    models = [optics.moment_model(pair, 2, ref[0]) for pair in raw[0] + raw[1]]
    controls = {
        i: (
            [c.ri(sp.Rational(x)) for x in arm["center_coordinates"]],
            c.ri(sp.Rational(arm["direction_error_with_serialization_upper"])),
        )
        for i, arm in enumerate(arms)
    }
    _, changes = optics.odd_pair(models, ref, controls)
    slopes = [changes[i][1] for i in range(len(arms))]
    for arm, (h, g) in zip(arms, slopes, strict=True):
        arm["odd_expansion_derivative"] = str(h)
        arm["odd_focusing_derivative"] = str(g)
    minors = {}
    for i, j in combinations(range(4), 2):
        hi, gi = slopes[i]
        hj, gj = slopes[j]
        det = hi * gj - hj * gi
        if 0 not in det:
            minors[f"{i},{j}"] = str(det)
    assert minors, slopes
    return dict(
        held_probability_and_material_rows=26,
        reference_occupations=list(map(str, probabilities)),
        reference_occupation_rates=list(map(str, rates)),
        dwell_clock=dwell,
        fixed_occupation_and_rate_controls=arms,
        nonzero_two_control_minors=minors,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(out)
    else:
        print(out, end="")
