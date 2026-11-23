"""Continuous-window home bound for every port-dephased release on the 20-cell board."""

from fractions import Fraction
from math import ceil, factorial, isqrt

import numpy as np
from scipy import sparse

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .gateway_transport_bound import exchange_step_counts


SCALE = 1 << 26
ORDER = 12
# Half-integer centers of the 23 unit intervals covering [15,38].
CENTERS_TWICE = range(31, 76, 2)


def _exchange(board):
    step = exchange_step_counts(board)
    if board.cells != 20:
        raise ArithmeticError("The confirmation uses the 20-cell carrier")
    identity = sparse.eye(board.n, format="csr")
    if (2*candidate.generator(board) + step != (6+12*candidate.KAPPA)*identity).nnz:
        raise ArithmeticError("Exchange half-step differs from the native one-load evolution")
    # In units of the off-diagonal Hamiltonian, the crossing block at each home
    # is one half of a partial permutation. Thus ||i[H,P_home]|| <= 1/2.
    entries = step.tocoo()
    crossing = entries.row//12 != entries.col//12
    if (np.any(entries.data[crossing] != 1)
            or np.any(np.bincount(entries.row[crossing], minlength=board.n) > 1)):
        raise ArithmeticError("The crossing seams must be a unit partial matching")
    return step


def _half_step_error(size):
    """Column-norm error in exp(i M/4) acting on a column of norm at most two."""
    radius = Fraction(3, 2)  # M = 6I-L has operator norm at most 6.
    half_unit = Fraction(isqrt(2*size-1)+1, 2*SCALE)
    term_error = Fraction(0)
    sum_error = Fraction(0)
    for k in range(1, ORDER+1):
        term_error = radius*term_error/k + half_unit
        sum_error += term_error
    tail = radius**(ORDER+1)/factorial(ORDER+1)/(1-radius/(ORDER+2))
    return sum_error + 2*tail


def _half_step(step, re, im):
    """Apply the order-12 Taylor polynomial, rounding each term to the dyadic grid."""
    term_re, term_im = re, im
    out_re, out_im = re.copy(), im.copy()
    for k in range(1, ORDER+1):
        # Adding half the denominator before floor division rounds ties upward.
        term_re, term_im = ((-(step @ term_im)+2*k)//(4*k),
                            ((step @ term_re)+2*k)//(4*k))
        out_re += term_re
        out_im += term_im
    return out_re, out_im


def dephased_home_ceiling():
    """Bound every sharp-start home probability at every real exposure in [15,38]."""
    board = federation(0)
    step = _exchange(board)
    error = max(CENTERS_TWICE)*_half_step_error(board.n)
    # Exact time steps are unitary. While total column error is <1 the input
    # norm to each approximate step is <=2, so its local error adds by induction.
    if error >= 1:
        raise ArithmeticError("Error budget exceeds the input norm bound")
    # A Taylor term has norm at most 2*(3/2)**k/k! + term_error < 4.
    # Sparse multiplies use six counts per row. Each terminal home sums the
    # squares of twelve complex entries of modulus at most two.
    if max(24*SCALE+4*ORDER, 4*(ORDER+1)*SCALE, 24*(2*SCALE)**2) > np.iinfo(np.int64).max:
        raise ArithmeticError("Taylor accumulation exceeds the integer range")
    re = np.eye(board.n, dtype=np.int64)*SCALE
    im = np.zeros_like(re)
    largest = 0
    for half_ticks in range(1, max(CENTERS_TWICE)+1):
        re, im = _half_step(step, re, im)
        if half_ticks in CENTERS_TWICE:
            home = (re*re+im*im).reshape(board.cells, 12, board.n).sum(axis=1)
            largest = max(largest, int(home.max()))
    midpoint = Fraction(largest, SCALE**2)
    # A norm-epsilon column perturbation costs <= 2 epsilon + epsilon^2 in
    # any projector probability. A center is within 1/2 of every exposure
    # in its interval, giving another (1/2)*(1/2) from the commutator.
    ceiling = midpoint + error*(2+error) + Fraction(1, 4)
    return midpoint, error, ceiling


if __name__ == "__main__":
    midpoint, error, ceiling = dephased_home_ceiling()
    print(f"Largest grid home weight {midpoint}")
    print(f"Column norm error <= {ceil(error*10**12)}/10^12")
    print(f"All 240 sharp starts, all 20 homes, exposure [15,38]: "
          f"joint home probability <= {ceil(ceiling*10**12)}/10^12")
