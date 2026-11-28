"""Bound the exact time-four sharp-head covariance using integer arithmetic."""

import argparse
from fractions import Fraction as F
import json
from math import factorial, isqrt
from pathlib import Path

import numpy as np
from scipy import sparse

from model.experiments.interval_echo_einstein_station.mobile_collar_source import (
    measure as formation,
)

SCALE = 2**45
ORDER = 13


def interval(lo, hi):
    """Outward endpoints in billionths, computed as integers."""
    unit = 10**9
    return [int((lo * unit) // 1), int(-((-hi * unit) // 1))]


def certify(sign):
    board = formation.federation(1)
    route = formation.region(board)
    h = formation.hamiltonian(board, route, sign, 8, center=False)
    # All local blocks come from integer embeddings of 4H divided by eight.
    # Without the floating trace subtraction 8H is exactly representable.
    scaled = 8 * h
    integer = scaled.astype(np.int64)
    assert (integer != scaled).nnz == 0
    diagonal = integer.diagonal()
    offset = int((diagonal.min() + diagonal.max()) // 2)
    b = integer - offset * sparse.eye(h.shape[0], dtype=np.int64, format="csr")
    assert (b != b.T).nnz == 0
    bound = int(abs(b).sum(axis=1).max())  # Symmetry bounds the spectral norm as well.
    steps = 1 << (bound - 1).bit_length()
    x = F(bound, 2 * steps)  # Four time units divided by 8 and by the step count.
    assert x <= F(1, 2)
    dimension = b.shape[0]
    # Each complex component is rounded to the nearest signed integer grid point.
    rounding = F(isqrt(2 * dimension) + 1, 2 * SCALE)
    tail = x ** (ORDER + 1) / factorial(ORDER + 1)
    # Unitary Taylor's integral remainder is at most 'tail' times the
    # input norm. A term-rounding error enters later terms with total norm
    # multiplier at most exp(x)<2. Telescoping uses unitary propagation.
    error = steps * (2 * tail + 2 * ORDER * rounding)
    assert error < 1  # Induction gives input norm <= 1+error < 2 in every step.
    r = np.zeros(dimension, dtype=np.int64)
    i = r.copy()
    lhs, rhs = formation.starting_coins(route, sign)
    r[4 * route[0] + 2 * lhs + rhs] = SCALE
    for _ in range(steps):
        tr, ti = r.copy(), i.copy()
        for degree in range(1, ORDER + 1):
            den = 2 * steps * degree
            limit = max(abs(tr).max().item(), abs(ti).max().item())
            assert bound * limit + den < 2**63
            tr, ti = ((b @ ti + den // 2) // den, (-(b @ tr) + den // 2) // den)
            assert (
                max(abs(r).max().item(), abs(i).max().item())
                + max(abs(tr).max().item(), abs(ti).max().item())
                < 2**63
            )
            r += tr
            i += ti
    sums = [
        sum(int(a) ** 2 + int(bi) ** 2 for a, bi in zip(r[j::4], i[j::4], strict=True))
        for j in range(4)
    ]
    p = [F(value, SCALE * SCALE) for value in sums]
    # For every coordinate projector P, |v*Pv-u*Pu| <= error*(2+error).
    eps = error * (2 + error)
    joint, left, right = p[3], p[2] + p[3], p[1] + p[3]
    assert min(left, right) > eps
    lo = (joint - eps) - (left + eps) * (right + eps)
    hi = (joint + eps) - (left - eps) * (right - eps)
    assert lo > 0
    return dict(
        source_sign=sign,
        cells=board.cells,
        dimension=dimension,
        time=4,
        eight_h_scalar_removed=offset,
        maximum_absolute_row_sum=bound,
        steps=steps,
        taylor_order=ORDER,
        integer_scale=SCALE,
        state_error_ceiling_billionths=interval(0, error)[1],
        covariance_billionths=interval(lo, hi),
        native_energy_increment_billionths=interval(17 * lo, 17 * hi),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = [certify(sign) for sign in (-1, 1)]
    path = Path(__file__).with_name("formation_certificate.json")
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert json.loads(path.read_text()) == result
