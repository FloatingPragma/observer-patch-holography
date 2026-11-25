"""Calibrate entropy and material-force responses using complete return-time counts."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from .. import certificate as cert
from .. import station
from . import clock_history, invisible
from ...directional_stress.connection.force_gradients import force_effects, force_gradients

HERE = Path(__file__).resolve().parent


def entropy_gradient(record, *, log_coordinates=cert.prior.entropy_coordinates):
    """Modular score at the independent stationary reference, on both correlation blocks."""
    data = record["fixed_state_centers"]
    axes = [cert.prior.sym_matrix(op) for op in station.iso.basis(False)]
    products = [cert.prior.sym_matrix(op) for op in station.iso.basis(True)]
    entries, errors = [], []
    for coordinates in data["even_then_odd"]:
        center = sp.eye(4) / 4
        for coefficient, axis in zip(coordinates, axes, strict=True):
            center += sp.Rational(coefficient, data["denominator"]) * axis
        gradient, log_error = log_coordinates(center, 0, sp.Rational(".18"), products)
        entries.extend(gradient)
        errors.extend([2 * log_error] * len(products))
    # The scalar log of the parity weight pairs to zero on each traceless tangent.
    reference_error = cert.ri(record["conditional_HS_error_upper"])
    uncertainty = invisible.frobenius(sp.Matrix(errors)) + cert.iv.sqrt(2) * reference_error / cert.ri(".18")
    return sp.Matrix(entries), uncertainty


def upper(bound):
    return float(np.nextafter(float(bound.b), np.inf))


def freeze(matrix, inverse_transpose, input_error, target, target_error, direction, correction):
    """Freeze dyadic bin weights, bound their error on all tangents, and test a held change."""
    denominator = 2**52
    numerators = [round(item * denominator) for item in inverse_transpose * target]
    weights = sp.Matrix(numerators) / denominator
    response = matrix.T * weights
    clock_error = invisible.frobenius(weights) * input_error
    norm = invisible.frobenius
    residual = norm(response - target) + clock_error + target_error
    # Evaluate the actual clock observable independently of the target's error.
    center = response.dot(direction)
    radius = norm(response) * correction + clock_error * (norm(direction) + correction)
    measured = cert.ri(center) + cert.iv.mpf([-1, 1]) * radius
    assert measured.a > 0
    return {
        "bin_score_numerators": [str(n) for n in numerators],
        "bin_score_denominator": str(denominator),
        "bin_score_l2_norm_upper": upper(norm(weights)),
        "target_error_per_unit_HS_change_upper": upper(residual),
        "clock_score_on_held_tangent": str(measured),
    }


def absolute_forces():
    """Predict absolute force on every normalized state, including unknown parity and marginals."""
    full = json.loads((HERE / "full_clock.json").read_text())
    affine = full["affine_clock"]
    denominator = int(affine["denominator"])
    panel = sp.Matrix([[int(v) for v in row] for row in affine["numerators"]]) / denominator
    offset = sp.Matrix([int(v) for v in affine["offset_numerators"]]) / denominator
    inverse = panel.T.inv()
    eps = cert.ri(sp.Rational(clock_history.timing_error()[0]))
    axes = [cert.prior.sym_matrix(op) for op in station.iso.basis(False)]
    scores = {}
    for name, blocks in zip(("vacancy", "excess"), force_effects(), strict=True):
        even, odd = blocks
        target = sp.Matrix([cert.inner(op, a) for op in blocks for a in axes]
                           + [(odd.trace() - even.trace()) / 4])
        baseline = (odd.trace() + even.trace()) / 8
        scale = 2**52
        numerators = [round(v * scale) for v in inverse * target]
        w = sp.Matrix(numerators) / scale
        intercept = round((baseline - offset.dot(w)) * scale)
        # For a normalized physical state, ||x||² = ||rho-I/8||_HS² + z²/2 <= 1.
        # Each event offset costs eps, and each calibrated basis effect costs 2 eps.
        residual = invisible.frobenius(panel.T * w - target)
        residual += cert.ri(abs(sp.Rational(intercept, scale) + offset.dot(w) - baseline))
        residual += (62 + cert.iv.sqrt(31)) * eps * invisible.frobenius(w)
        residual += cert.ri(sp.Rational("2e-15"))  # uniform material expectation error
        bound = upper(residual)
        assert bound < 1e-9
        scores[name] = {
            "weight_numerators": [str(n) for n in numerators],
            "intercept_numerator": str(intercept),
            "denominator": str(scale),
            "clock_weight_l2_norm_upper": upper(invisible.frobenius(w)),
            "absolute_error_all_normalized_heads_upper": bound,
        }
    return {"event_order": full["kept_clock_bins"], "signed_force": scores}


def run():
    cert.iv.dps = 32
    history = json.loads((HERE / "clock_history.json").read_text())
    record = json.loads((HERE / "invisible.json").read_text())
    central = history["pooled_timestamp_matrix"]
    panel = sp.Matrix([[int(i) for i in row] for row in central["numerators"]]) / int(central["denominator"])
    inverse_transpose = panel.T.inv()
    matrix_error = 18 * cert.ri(central["entry_error_upper"])
    direction = sp.Matrix([sp.Rational(x) for x in record["tangent_center_even_then_odd"]])
    correction = cert.ri(sp.Rational(record["true_kernel_correction_HS_upper"]) + sp.Rational(1, 10**14))

    def calibrate(target, error):
        return freeze(panel, inverse_transpose, matrix_error, target, error, direction, correction)

    entropy = calibrate(*entropy_gradient(record))
    assert entropy["target_error_per_unit_HS_change_upper"] < 2.1e-5
    force = {}
    gradients, row_error = force_gradients()
    for i, name in enumerate(("vacancy", "excess")):
        force[name] = calibrate(gradients.row(i).T, cert.ri(row_error))
        assert force[name]["target_error_per_unit_HS_change_upper"] < 1e-10
    return {
        "reference": "volume-matched head state stationary under the uniform 24-setting reuse channel",
        "event_order": history["kept_clock_bins"],
        "normalization": "original launch denominator at each directed station pair, charge pooled equally",
        "entropy_differential": entropy,
        "independent_material_force": force,
        "whole_algebra_absolute_force": absolute_forces(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    content = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(content)
    else:
        print(content, end="")
