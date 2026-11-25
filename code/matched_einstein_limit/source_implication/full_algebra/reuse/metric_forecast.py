"""Forecast the six separately surveyed lengths from one directed relay and return."""

import argparse
import json
from pathlib import Path

from fractions import Fraction

import numpy as np
import sympy as sp

from .. import certificate as cert, native_precision, station
from . import (
    clock_scores,
    duration_geometry as geometry,
    force_forecast,
    invisible,
    measure,
)

HERE = Path(__file__).resolve().parent


def targets():
    numerators, unit, error = native_precision.exact_effects()
    directed = dict(
        zip(measure.ROUTES, cert.rational_effects(numerators, unit), strict=True)
    )
    axes = [cert.prior.sym_matrix(op) for op in station.iso.basis(False)]
    for edge in geometry.EDGES:
        blocks = [
            (f + r) / 2
            for f, r in zip(directed[edge], directed[edge[::-1]], strict=True)
        ]
        even, odd = blocks
        slopes = [cert.inner(block, axis) for block in blocks for axis in axes]
        slopes.append((odd.trace() - even.trace()) / 4)
        yield sp.Matrix(slopes), (odd.trace() + even.trace()) / 8, sp.Rational(error)


def fit(menu, denom, observables):
    matrix = sp.Matrix([[int(v) for v in row] for row in menu["matrix"][:36]]) / denom
    offsets = sp.Matrix([int(v) for v in menu["offsets"][:36]]) / denom
    # Solve the full-rank least-norm calibration on exact native moments; the
    # frozen dyadic score is reproducible across BLAS implementations.
    inverse = matrix * (matrix.T * matrix).inv(method="DM")
    errors = [sp.Rational(str(v)) for v in menu["event_effect_error_upper"][:36]]
    readers = []
    scale = 2**48
    for target, origin, target_error in observables:
        proposal = inverse * target
        nums = [round(x * scale) for x in proposal]
        w = sp.Matrix(nums) / scale
        intercept = round((origin - w.dot(offsets)) * scale)
        residual = invisible.frobenius(matrix.T * w - target)
        residual += cert.ri(
            abs(sp.Rational(intercept, scale) + offsets.dot(w) - origin) + target_error
        )
        # For normalized states ||x||_2 <= 1. Each event has offset error <=eps
        # and 31 slopes with individual error <=2eps, including the parity axis.
        residual += (1 + 2 * cert.iv.sqrt(31)) * cert.ri(
            sum(abs(s) * e for s, e in zip(w, errors, strict=True))
        )
        bound = clock_scores.upper(residual)
        assert bound < 1e-5
        readers.append(
            dict(
                denominator=str(scale),
                score_numerators=[str(n) for n in nums],
                intercept_numerator=str(intercept),
                uniform_error_upper=bound,
            )
        )
    return readers


def direct_effects(labels):
    """Independently propagate all launches and retain the original denominator."""
    board = station.iso.pairing.formation.federation(1)
    ports = station.iso.pairing.formation.region(board)
    local = np.zeros((2, len(labels), 2, 4, 4), complex)
    survey = np.tile(3 * np.eye(4, dtype=complex), (6, 2, 1, 1))
    for qi, charge in enumerate((-1, 1)):
        for parity in (0, 1):
            h = station.iso.ranging.generator(board, ports, charge, 8, bool(parity))
            panels = station.iso.ranging.kernels(h, ports, 0.5, 12)
            hits = {
                r: station.iso.ranging.first_hits(panels, *r) for r in measure.ROUTES
            }
            for j, edge in enumerate(geometry.EDGES):
                for route in (edge, edge[::-1]):
                    for t in range(2, 12):
                        for a in range(1, t):
                            amp = hits[route[::-1]][t - a] @ hits[route][a]
                            effect = amp.conj().T @ amp
                            survey[j, parity] -= (12 - t) * effect / 16
                            if route == (2, 1) and t <= 9:
                                local[qi, labels.index([a, t]), parity] = effect
    return local, survey


def run():
    cert.iv.dps = 32
    archive = json.loads((HERE / "edge_clock_history.json").read_text())
    labels = archive["event_labels"][:36]
    observables = list(targets())
    readers = [
        fit(menu, int(archive["affine_denominator"]), observables)
        for menu in archive["native_affine_clock"]
    ]
    local, survey = direct_effects(labels)
    observations = {}
    axes = [cert.prior.sym_matrix(op) for op in station.iso.basis(False)]
    for name, rho in force_forecast.preparations().items():
        # All three native Pauli lotteries have dyadic entries. Remove only
        # floating representation of the preparation, before any observation.
        def dyadic(value):
            return sp.Rational(Fraction(float(value)).limit_denominator(128))

        blocks = [
            sp.Matrix(
                4, 4, lambda i, j: dyadic(op[i, j].real) + sp.I * dyadic(op[i, j].imag)
            )
            for op in rho
        ]
        assert (
            max(np.max(abs(np.array(b, complex) - op)) for b, op in zip(blocks, rho))
            < 1e-14
        )
        x = sp.Matrix(
            [cert.inner(b, a) for b in blocks for a in axes]
            + [blocks[1].trace() - sp.Rational(1, 2)]
        )
        exact_lengths = [origin + t.dot(x) for t, origin, _ in observables]
        lengths = np.einsum("ebij,bji->e", survey, rho).real
        probs = np.einsum("qabij,bji->qa", local, rho).real
        forecasts = []
        float_forecasts = []
        float_error = []
        for qi, (p, panel) in enumerate(zip(probs, readers, strict=True)):
            menu = archive["native_affine_clock"][qi]
            denom = int(archive["affine_denominator"])
            matrix = sp.Matrix([[int(v) for v in row] for row in menu["matrix"][:36]])
            offset = sp.Matrix([int(v) for v in menu["offsets"][:36]])
            native = (offset + matrix * x) / denom
            independent = []
            amplification = []
            estimated = []
            for target, certified_target, row in zip(
                lengths, exact_lengths, panel, strict=True
            ):
                w = np.array(
                    [int(v) / int(row["denominator"]) for v in row["score_numerators"]]
                )
                result = (
                    int(row["intercept_numerator"]) / int(row["denominator"]) + w @ p
                )
                exact_w = sp.Matrix([int(v) for v in row["score_numerators"]]) / int(
                    row["denominator"]
                )
                certified = sp.Rational(
                    int(row["intercept_numerator"]), int(row["denominator"])
                ) + exact_w.dot(native)
                assert abs(certified - certified_target) < row["uniform_error_upper"]
                assert abs(float(certified_target) - target) < 1e-12
                # Independent numerical propagation of an ill-conditioned reader
                # can magnify eigensolver errors. Record this discrepancy separately.
                delta_p = np.abs(np.asarray(native, dtype=float).ravel() - p)
                budget = abs(w) @ delta_p + np.finfo(float).eps * 64 * (
                    abs(w) @ abs(p) + 4
                )
                assert abs(result - float(certified)) <= budget
                assert abs(result - target) < 5e-8
                estimated.append(float(certified))
                independent.append(float(result))
                amplification.append(float(budget))
            float_forecasts.append(independent)
            float_error.append(amplification)
            forecasts.append(estimated)
        matrix = np.array(geometry.gram(lengths**2), dtype=float)
        observations[name] = dict(
            separate_lengths=lengths.tolist(),
            single_edge_predictions=forecasts,
            independent_float_predictions=float_forecasts,
            float_probability_and_roundoff_budgets=float_error,
            separate_volume=float(np.sqrt(np.linalg.det(matrix)) / 6),
        )
    return dict(
        launch_port=38,
        relay_port=1,
        joint_events=labels,
        species=[-1, 1],
        edges=geometry.EDGES,
        absolute_readers=readers,
        held_population_comparisons=observations,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    content = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(content)
    else:
        print(content, end="")
