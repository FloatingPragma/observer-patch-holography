"""Geometric bounds and opposite stochastic time order on the native four-station radar."""

import argparse
from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
import sympy as sp

from .. import certificate as cert
from ...directional_stress.connection import force_entropy
from . import clock_history, clock_scores, invisible, measure

HERE = Path(__file__).resolve().parent
EDGES = tuple(combinations(range(4), 2))
EDGE = (1, 2)
# Geometry-independent lower bounds on the original-denominator timeout probability.
TIMEOUT = tuple(sp.Rational(v, 1000) for v in (630, 993, 999, 836, 994, 725))
# Proposal only. The physical null-space projection and all signs are verified below.
HELD = (
    -23449464,
    3280032,
    1048050,
    -65900256,
    -3942453,
    159475,
    -1772164,
    -1468353,
    14273,
    5056496,
    11396474,
    -810798,
    70085021,
    -1142406,
    148836,
    1532880,
    1692371,
    -14880,
)


def effects():
    """Cumulative return effects for each route, species, sector and deadline 2..12."""
    summed = {}
    for q, b, i, t, matrix in clock_history.native_events(last_tick=12):
        end = (sp.Matrix(matrix[0]) + sp.I * sp.Matrix(matrix[1])) / clock_history.F
        end = (end + end.H) / 2
        key = q, b, measure.ROUTES[i]
        records = summed.setdefault(key, [])
        assert len(records) == t - 2
        records.append(end + (records[-1] if records else sp.zeros(4)))
    return summed


def gram(squared):
    distances = sp.zeros(4)
    for (i, j), value in zip(EDGES, squared, strict=True):
        distances[i, j] = distances[j, i] = value
    return sp.Matrix(
        3,
        3,
        lambda i, j: (
            distances[0, i + 1] + distances[0, j + 1] - distances[i + 1, j + 1]
        )
        / 2,
    )


def metric_bound(cumulative, bin_error):
    for edge, lower in zip(EDGES, TIMEOUT, strict=True):
        for q, b in product(range(2), repeat=2):
            complete = (
                cumulative[q, b, edge][-1] + cumulative[q, b, edge[::-1]][-1]
            ) / 2
            assert cert.global_survey.positive_definite(
                (1 - lower - 11 * bin_error) * sp.eye(4) - complete
            )
    eigen_floor = sp.Rational("1.55")
    determinants = []
    for corner in product((0, 1), repeat=6):
        squared = [9 if bit else (3 * lo) ** 2 for bit, lo in zip(corner, TIMEOUT)]
        matrix = gram(squared)
        assert cert.global_survey.positive_definite(matrix - eigen_floor * sp.eye(3))
        determinants.append(matrix.det())
    assert min(determinants) > (6 * sp.Rational("1.57")) ** 2
    return {
        "edges": EDGES,
        "minimum_timeout_probability": [float(v) for v in TIMEOUT],
        "minimum_anchored_gram_eigenvalue_lower": str(eigen_floor),
        "smallest_corner_determinant": str(min(determinants)),
        "volume_lower": "1.57",
    }


def traceless_log_coordinates(center, input_error, lower, products):
    """Exact rational series for -log(center), modulo the identity."""
    assert input_error == 0 and all(op.trace() == 0 for op in products)
    eye = sp.eye(4)
    # Y = I - 4C, so the traceless part of -log(C) is sum Y^n/n.
    y = eye - 4 * center
    bound = sp.Rational(1, 2)
    assert cert.global_survey.positive_definite(center - lower * eye)
    assert all(cert.global_survey.positive_definite(bound * eye + sign * y)
               for sign in (-1, 1))
    terms = 40
    matrix = y.to_DM().to_field()
    domain = matrix.domain
    power = matrix
    series = matrix
    for n in range(2, terms + 1):
        power = power * matrix
        series = series + power.scalarmul(domain.convert(sp.Rational(1, n)))
    tail = bound ** (terms + 1) / ((terms + 1) * (1 - bound))
    score = series.to_Matrix()
    return sp.Matrix([cert.inner(score, op) for op in products]), tail


def physical_control(corr):
    ranges, range_error = force_entropy.range_gradients(corr)
    force = json.loads((force_entropy.HERE / "force_entropy.json").read_text())
    rows = sp.Matrix(force["force_gradient_numerators"]) / int(
        force["force_gradient_denominator"]
    )
    d, correction, _, _ = force_entropy.corrected_kernel_direction(
        ranges,
        range_error,
        rows,
        sp.Rational(force["force_gradient_error_per_row"]),
        proposal=sp.Matrix(HELD) / 10**8,
    )
    source = json.loads((HERE / "invisible.json").read_text())
    g, ge = clock_scores.entropy_gradient(source, log_coordinates=traceless_log_coordinates)
    norm = invisible.frobenius(d) + correction
    radius = invisible.frobenius(g) * correction + ge * norm
    slope = cert.ri(g.dot(d)) + cert.iv.mpf([-1, 1]) * radius
    assert slope.a > cert.ri("8e-6") and norm.b < cert.ri("1.00001")
    step = cert.ri("1e-7")
    # The source has minimum eigenvalue >.024 on this chord, by the previously
    # verified parity-weight and conditional-state enclosures on this same board.
    increase = step * slope.a - step**2 * norm.b**2 / (2 * cert.ri(".024"))
    assert increase.a > cert.ri("5e-13")
    assert cert.ri(source["lottery_slack_lower"]) - 12 * cert.ri(".01") * norm.b > 0
    return (
        d,
        correction,
        norm,
        {
            "entropy_slope": str(slope),
            "finite_entropy_gain_at_1e-7_lower": str(increase.a),
            "kernel_correction_upper": str(correction.b),
        },
    )


def timing_order(cumulative, bin_error, corr, d, correction, norm):
    bounds = {}
    identity = sp.eye(4)
    xz = 2 * corr[2]
    for q, name in enumerate(("vacancy", "excess")):
        for route in (EDGE, EDGE[::-1]):
            simple_bounds, held_bounds = [], []
            cdfs = [[], []]
            for i in range(11):
                block = [cumulative[q, b, route][i] for b in range(2)]
                err = (i + 1) * bin_error
                sign = 1 if q == 0 else -1
                row = sp.Matrix(
                    [cert.inner(op, basis) for op in block for basis in corr]
                )
                simple = cert.ri_band(sign * row[11], 2 * err)
                # The combined two-block error has HS norm at most sqrt(8)*err.
                radius = (
                    invisible.frobenius(row) * correction
                    + cert.iv.sqrt(8) * cert.ri(err) * norm
                )
                held = cert.ri(sign * row.dot(d)) + cert.iv.mpf([-1, 1]) * radius
                assert simple.a > cert.ri("5.13e-6") and held.a > cert.ri("1e-8")
                simple_bounds.append(str(simple))
                held_bounds.append(str(held))
                for arm, orientation in enumerate((-1, 1)):
                    p = cert.inner(block[0], 3 * identity / 16) + cert.inner(
                        block[1], (identity + orientation * xz) / 16
                    )
                    cdfs[arm].append(float(p))
            bounds[name + "_" + str(route[0]) + "_" + str(route[1])] = {
                "simple_signed_slopes": simple_bounds,
                "held_signed_slopes": held_bounds,
                "cdf_centers": cdfs,
            }
    return bounds


def counts(timing, bin_error):
    rng = np.random.default_rng(2551126)
    trials = 2**46
    count_bound = 2 * np.sqrt(np.log(176 / 0.01) / (2 * trials))
    out = {}
    for name, row in timing.items():
        observed = []
        for cdf in row["cdf_centers"]:
            rates = np.diff([0.0, *cdf, 1.0])
            assert rates.min() >= 0
            observed.append(rng.multinomial(trials, rates))
        # The native enclosures also cover conversion to binary64 for the sampler.
        half_width = count_bound + 2 * (11 * float(bin_error) + 1e-13)
        diff = (observed[1][:-1] - observed[0][:-1]).cumsum() / trials
        oriented = diff if name.startswith("vacancy") else -diff
        assert min(oriented) > half_width
        out[name] = {
            "minus_counts": observed[0].tolist(),
            "plus_counts": observed[1].tolist(),
            "smallest_signed_difference": float(min(oriented)),
        }
    return {
        "trials_per_arm_and_setting": trials,
        "simultaneous_99pct_half_width": float(half_width),
        "records": out,
    }


def run():
    cert.iv.dps = 32
    cumulative = effects()
    error = sp.Rational(clock_history.timing_error(last_tick=12)[0])
    corr = [cert.prior.sym_matrix(op) for op in cert.iso.basis(True)]
    d, correction, norm, source = physical_control(corr)
    timing = timing_order(cumulative, error, corr, d, correction, norm)
    # The squared-contrast-off control makes the right even contact scalar on
    # its coin. The other marked contact acts on the left coin only. Every
    # clock effect therefore factors through I on the right, for either parity.
    for q in (-1, 1):
        blocks = invisible.station.iso.ranging.parity_blocks(q, 0, 1)
        np.testing.assert_array_equal(blocks, blocks[..., 0, 0, None, None] * np.eye(2))
    return {
        "ports": [0, 1, 38, 39],
        "deadlines": list(range(2, 13)),
        "bin_operator_error": str(Fraction(error)),
        "common_metric_box": metric_bound(cumulative, error),
        "matched_source": source,
        "opposite_time_orders": timing,
        "simulated_duration_counts": counts(timing, error),
        "zero_contrast_correlation_clock_slope": 0,
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
