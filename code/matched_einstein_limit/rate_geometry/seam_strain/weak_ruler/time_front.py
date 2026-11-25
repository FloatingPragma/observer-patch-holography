"""Reciprocal leading returns at different ages of the same autonomous pump."""

import argparse
import json
from functools import lru_cache
from math import factorial, isqrt

import numpy as np
import sympy as sp

from . import flight, late_read, observability
from .. import live_source
from ..source import certificate as cert

# The choices are made on front sensitivities to the nine Pauli coordinates.
# They use neither the regional logarithm nor a cut current.
REPEATED = (1, 2, 3)
SCALE = 1 << 48
LATER_STEP = sp.Rational(1, 8)


def rotate(wave):
    """Multiply the packed real and imaginary columns by minus i."""
    re, im = np.hsplit(wave, 2)
    return np.column_stack((im, -re))


def retained_layout(wave):
    bank = wave.reshape(2, -1, 4, 16).transpose(1, 0, 2, 3).reshape(-1, 8, 16)
    return [a.astype(object) for a in np.split(bank, 2, axis=-1)]


def compressed(left, operator, right):
    """Exact left-adjoint times (identity on position tensor operator) times right."""
    ore, oim = (operator @ a for a in right)
    lr, li = (a.reshape(-1, 8).T for a in left)
    rr, ri = (a.reshape(-1, 8) for a in (ore, oim))
    return sp.Matrix(lr @ rr + li @ ri) + sp.I * sp.Matrix(lr @ ri - li @ rr)


def axes_response(effect):
    rail = sp.diag(2, 1) / 3
    return sp.Matrix(
        [
            [
                sp.re(sp.trace(effect * sp.kronecker_product(rail, axis)))
                for axis in observability.AXES
            ]
        ]
    )


def operators():
    return [
        (
            np.array([int(v) for v in op], dtype=object).reshape(8, 8),
            observability.norm_inf(op),
        )
        for op in late_read.centered_front_operators()
    ]


def scaled_action(h, wave, divisor):
    """Round a scaled exchange-generator action on each real component."""
    peak = int(abs(wave).max())
    assert int(abs(h).sum(axis=1).max()) * peak + divisor < 2**62
    return (rotate(h @ wave) + divisor // 2) // divisor


@lru_cache(maxsize=1)
def calibrated_waves():
    h = flight.hardware(-1)[0]
    hi, lo, column_error, _ = live_source.pump_center(True)
    # The exact two-limb center divided by 2^82 rounds to the 48-bit grid.
    scale_shift = (1 << 82) // SCALE
    base = np.asarray(
        ((hi.astype(object) * (1 << 36) + lo + scale_shift // 2) // scale_shift),
        dtype=np.int64,
    )
    rounding = sp.Rational(1 + isqrt(2 * h.shape[0] * 8), 2 * SCALE)
    base_error = 3 * column_error + rounding
    h_bound = sp.Rational(int(abs(h).sum(axis=1).max()), 40)
    waves = [base]
    errors = [base_error]
    for _ in range(2):
        waves.append(scaled_action(h, waves[-1], 40))
        errors.append(h_bound * errors[-1] + rounding)

    # Taylor in pump age about four. Unitarity bounds the integrated
    # remainder by L^(m+1)/(m+1)! for a normalized true starting isometry.
    duration = int(1 / LATER_STEP)
    length = h_bound / duration
    degree = 32
    step = base
    later = base.copy()
    accumulated = 0
    step_error = 0
    for j in range(1, degree + 1):
        step = scaled_action(h, step, 40 * duration * j)
        assert int(abs(later).max()) + int(abs(step).max()) < 2**62
        later += step
        step_error = length * step_error / j + rounding
        accumulated += step_error
    # The truncated exponential on the initial rounding error has norm
    # at most 1+tail, by comparison with its unitary exponential.
    tail = length ** (degree + 1) / factorial(degree + 1)
    later_error = accumulated + tail + (1 + tail) * base_error
    return (
        [retained_layout(w) for w in waves],
        errors,
        h_bound,
        retained_layout(later),
        later_error,
    )


def product_error(a, b, errors, bound):
    # The true j-th wave has operator norm <= ||H||^j.
    return errors[a] * bound**b + errors[b] * bound**a + errors[a] * errors[b]


def read(operator, den, left, right, symmetric=False):
    eff = compressed(left, operator, right)
    if symmetric:
        eff += eff.H
    return axes_response(eff / (SCALE**2 * den))


def rounded(row, error):
    center = late_read.rounded_matrix(row)
    return center, (row - center).applyfunc(abs) + sp.ones(*row.shape) * error


@lru_cache(maxsize=1)
def measurements():
    waves, errors, bound, later, later_error = calibrated_waves()
    z = np.diag([1] * 4 + [-1] * 4).astype(object)
    dz = read(z, 1, waves[1], waves[0], True)
    ddz = read(z, 1, waves[2], waves[0], True) + 2 * read(z, 1, waves[1], waves[1])
    zrows = [
        rounded(read(z, 1, waves[0], waves[0]), product_error(0, 0, errors, bound)),
        rounded(dz, 2 * product_error(1, 0, errors, bound)),
        rounded(
            ddz,
            2
            * (product_error(2, 0, errors, bound) + product_error(1, 1, errors, bound)),
        ),
    ]
    base, rate, delayed = [], [], []
    for i, (op, den) in enumerate(operators()):
        base.append(
            rounded(
                read(op, den, waves[0], waves[0]), product_error(0, 0, errors, bound)
            )
        )
        if i in REPEATED:
            rate.append(
                rounded(
                    read(op, den, waves[1], waves[0], True),
                    2 * product_error(1, 0, errors, bound),
                )
            )
            delayed.append(
                rounded(read(op, den, later, later), 2 * later_error + later_error**2)
            )
    return zrows, base, rate, delayed


def stack(rows):
    return tuple(sp.Matrix.vstack(*parts) for parts in zip(*rows, strict=True))


def rank_certificate(rows):
    center, widths = stack(rows)
    inverse = center.inv()
    defect = observability.norm_inf(inverse.applyfunc(abs) * widths)
    assert defect < 1
    return dict(
        rank=center.rows,
        inverse_defect_upper=cert.enclosed(defect, 0)[1],
        true_inverse_infinity_norm_upper=cert.enclosed(
            observability.norm_inf(inverse) / (1 - defect), 0
        )[1],
    )


def matched_front_readings(rows):
    """Delayed raw contrasts for the two native lotteries with identical age-four fronts."""
    record = json.loads((flight.HERE / "observability.json").read_text())
    direction = sp.Matrix([sp.Rational(v) for v in record["blind_direction_center"]])
    # The center is dyadic. The published correction bounds every component.
    error = sp.Rational(record["maximum_coordinate_correction_upper"]) + sp.Rational(
        1, 10**15
    )
    scale = record["native_lotteries"]["common_scale"]
    results = []
    for center, widths in rows:
        value = center.dot(direction)
        radius = sum(abs(v) for v in center) * error + sum(
            width * (abs(v) + error) for width, v in zip(widths, direction, strict=True)
        )
        results.append(cert.enclosed(value / scale, radius / scale))
    return results


def run():
    from . import jet_lottery

    zrows, base, rate, delayed = measurements()
    zmid, zwidth = stack(zrows)
    return dict(
        pump_time=4,
        later_pump_time=float(4 + LATER_STEP),
        units="native exchange time",
        axes=list(late_read.AXES),
        initial_links=list(observability.SELECTED),
        repeated_links=[observability.SELECTED[i] for i in REPEATED],
        with_front_derivatives=rank_certificate(base + rate),
        with_delayed_fronts=rank_certificate(base + delayed),
        nearest_strain_response_centers=np.array(zmid, float).tolist(),
        nearest_entry_errors_upper=[
            cert.enclosed(max(zwidth.row(i)), 0)[1] for i in range(3)
        ],
        matched_age_four_lotteries=dict(
            front_velocity_contrasts=matched_front_readings(rate),
            delayed_front_contrasts=matched_front_readings(delayed),
            normalization="orientation-averaged, identity-subtracted front coefficient divided by the operator row norm, per unit common admixture",
        ),
        matched_preparations=jet_lottery.witness(zmid, zwidth),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (flight.HERE / "time_front.json").write_text(output)
    else:
        print(output, end="")
