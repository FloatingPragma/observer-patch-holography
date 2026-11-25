"""Certify which registered four-station range menus determine the head energy and state."""
import argparse
from fractions import Fraction
import json
from math import isqrt
from pathlib import Path

import numpy as np

from ..full_algebra import native_precision as native
from ..full_algebra import station


CAPS = (2, 4, 6)
SCALE = 2**50


def integer_directions():
    """Four times each rational HS-orthogonal direction, including (+I,-I)/4."""
    directions = 4 * station.directions()
    directions[-1] = np.stack([np.eye(4), -np.eye(4)])
    assert np.array_equal(directions.real, directions.real.round())
    assert np.array_equal(directions.imag, directions.imag.round())
    return directions.real.astype(int).astype(object), directions.imag.astype(int).astype(object)


def resolved_range_matrices():
    """Dyadic traceless design and tracial effects, retaining the probe charge."""
    board = station.iso.pairing.formation.federation(1)
    route = station.iso.pairing.formation.region(board)
    ports = np.array([4*p+i for p in route for i in range(4)])
    coefficients = native.prior.weights()
    d_re, d_im = integer_directions()
    rows = {(q, cap): np.zeros((12, 31), dtype=object) for q in (-1, 1) for cap in CAPS}
    traces = {key: np.zeros(12, dtype=object) for key in rows}
    for parity in (0, 1):
        for charge in (-1, 1):
            h = native.generator(board, route, charge, bool(parity))
            kernels = native.kernel_panels(h, ports, coefficients)
            for cap in CAPS:
                effects = native.prior.direction_effects(kernels, 2*cap)
                for i, (re, im) in enumerate(effects):
                    # Real trace against a Hermitian direction ignores the anti-Hermitian part.
                    rows[charge, cap][i] += np.sum(
                        re.T*d_re[:, parity] - im.T*d_im[:, parity], axis=(1, 2)
                    )
                    traces[charge, cap][i] += re.trace()
    # Each nonconstant effect has denominator 4*FRAC; directions above were quadrupled.
    return rows, traces, 16*native.prior.FRAC, 2*native.prior.error_effect()


def range_matrices():
    rows, _, unit, error = resolved_range_matrices()
    pooled = {cap: rows[-1, cap]+rows[1, cap] for cap in CAPS}
    return pooled, 2*unit, error


def dyadic(array):
    return np.vectorize(lambda x: int(round(x*SCALE)), otypes=[object])(array)


def max_row_sum(matrix):
    return int(max(abs(matrix).sum(axis=1)))


def sqrt_upper(value):
    """Round upward to 1e-10 by integer arithmetic."""
    scale = 10**10
    return Fraction(isqrt(value.numerator*scale**2//value.denominator)+1, scale)


def rational(numerator, denominator):
    return Fraction(int(numerator), int(denominator))


def hidden_energy(matrix, unit, error):
    """Pin the odd probability derivative to one and solve all selected ranges exactly."""
    first = matrix[:, :-1]
    pivot = matrix[:, -1]
    f = first.astype(float)/unit
    target = pivot.astype(float)/unit
    numeric_inverse = np.linalg.pinv(f, rcond=1e-16)
    q = dyadic(numeric_inverse)
    v = dyadic(numeric_inverse @ target)
    k = len(first)
    identity = np.eye(k, dtype=int).astype(object)
    remainder = first @ q - identity * (unit*SCALE)
    # ||(A0 Q - I)||_inf for the true A0, with each uncertain matrix entry <= error.
    right_error = rational(max_row_sum(remainder), unit*SCALE)
    right_error += error * rational(abs(q).sum(), SCALE)
    assert right_error < 1
    defect = rational(abs(first @ v - pivot*SCALE).max(), unit*SCALE)
    defect += error * (1 + rational(abs(v).sum(), SCALE))
    # The inverse Neumann series corrects the range residual, leaving v[30]=-1 exactly.
    q_norm = sqrt_upper(rational(sum(abs(q).sum(axis=1)**2), SCALE**2))
    correction = q_norm * defect/(1-right_error)
    hs = sqrt_upper(rational(sum(v**2), SCALE**2) + Fraction(1, 2)) + correction
    lower = 1/hs
    return dict(right_inverse_defect_upper=float(right_error), witness_HS_upper=float(hs),
                normalized_energy_lower=float(lower), correction_HS_upper=float(correction))


def tomography(matrix, unit, error):
    """A verified left inverse supplies an exact linear dual for every head observable."""
    q = dyadic(np.linalg.pinv(matrix.astype(float)/unit, rcond=1e-16))
    m, n = matrix.shape
    identity = np.eye(n, dtype=int).astype(object)
    preliminary = q @ matrix - identity*(unit*SCALE)
    q_norm = rational(max_row_sum(q), SCALE)
    defect = rational(max_row_sum(preliminary), unit*SCALE) + n*q_norm*error
    assert defect < 1
    inverse_norm = q_norm/(1-defect)
    e_norm = rational(abs(q[-1]).sum(), SCALE)
    e_defect = rational(abs(preliminary[-1]).sum(), unit*SCALE) + n*e_norm*error
    # This bounds the 1-norm of the exact dual row for parked energy.
    energy_norm = e_norm + e_defect*inverse_norm
    return dict(range_rows=m, coordinates=n, left_inverse_defect_upper=float(defect),
                reconstruction_inf_norm_upper=float(inverse_norm),
                energy_dual_one_norm_upper=float(energy_norm))


def parked_terms():
    """Integer 99-state check on the blank and both two-coin reversal sectors."""
    emission = station.iso.ranging.formation.emission
    indices = emission.columns((1, 1))
    basis = np.zeros((emission.LAW.DIM, 5), dtype=int)
    basis[indices[0], 0] = 1
    for odd in (0, 1):
        for coin in (0, 1):
            for sign in (-1, 1):
                row = indices[1 + emission.LAW.LABELS.index((coin, sign))]
                basis[row, 1 + 2*odd + coin] = sign**odd
    exchange4 = emission.four_h(0)
    contrast32 = emission.four_h(8) - exchange4
    exchange = np.array([0, 0, 0, 1, 1])
    contrast = np.array([0, 1, 1, 1, 1])
    np.testing.assert_array_equal(exchange4 @ basis, 4 * basis * exchange)
    np.testing.assert_array_equal(contrast32 @ basis, 32 * basis * contrast)
    return dict(local_column_order=["blank", "even0", "even1", "odd0", "odd1"],
                exchange=exchange.tolist(), contrast=contrast.tolist())


def verify():
    station.iso.ranging.ordinary_write()
    rows, unit, error = range_matrices()
    hidden = {}
    for caps, bound in (((6,), .0249), ((2, 4), .00125), ((2, 6), .000584), ((4, 6), .0002317)):
        data = hidden_energy(np.vstack([rows[c] for c in caps]), unit, error)
        assert data["normalized_energy_lower"] > bound
        hidden["+".join(map(str,caps))] = data
    complete = tomography(np.vstack([rows[c] for c in CAPS]), unit, error)
    assert complete["left_inverse_defect_upper"] < .1
    return dict(uniform_design_entry_error_upper=float(error), invisible=hidden,
                complete=complete, parked_terms=parked_terms())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = json.dumps(verify(), indent=2, allow_nan=False)+"\n"
    if args.write:
        Path(__file__).with_name("results.json").write_text(result)
    else:
        print(result, end="")
