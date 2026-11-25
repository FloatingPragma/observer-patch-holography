"""Rational error certificate for all-history comparisons on the formed 80-cell source."""

from fractions import Fraction
from math import factorial, isqrt

import numpy as np
import sympy as sp
from scipy import sparse

from ...mobile_collar_source import measure as formation
from .measure import parity_blocks


def integer_flight(board, route, charge=1):
    """Twice the even-head generator for the specified charge, modulo an integer scalar."""
    arms = [tuple(sorted(route[:2])), tuple(sorted(route[2:]))]
    blank = np.array([e for e in board.seams if e not in arms])
    a, b = blank.T
    degree = np.bincount(blank.ravel(), minlength=board.n)
    diagonal = 17*degree + 16*(2*charge + charge**2)*(6-board.degree)
    adjacency = sparse.csr_matrix((np.ones(2*len(a), dtype=np.int64),
                                   (np.r_[a,b], np.r_[b,a])), shape=(board.n, board.n))
    h = sparse.kron(sparse.diags(diagonal, dtype=np.int64)-adjacency, np.eye(4,dtype=np.int64),
                    format="csr")
    blocks = 2*parity_blocks(charge,8,1)
    np.testing.assert_array_equal(blocks, blocks.astype(np.int64))
    blocks = blocks.astype(np.int64)
    for slot, ends in enumerate(arms):
        for i in range(2):
            for j in range(2):
                block = blocks[i+1,j+1] - int(i == j)*blocks[0,0]
                local = [np.eye(2, dtype=np.int64), np.eye(2,dtype=np.int64)]
                local[slot] = block
                dyad = sparse.csr_matrix(([1],([ends[i]],[ends[j]])), shape=(board.n,board.n))
                h += sparse.kron(dyad, np.kron(*local), format="csr")
    center = (int(h.diagonal().min()) + int(h.diagonal().max())) // 2
    return h - center*sparse.eye(h.shape[0],dtype=np.int64,format="csr")


def certificate():
    board = formation.federation(1)
    route = formation.region(board)
    h = integer_flight(board, route)
    row = int(abs(h).sum(axis=1).max())
    assert row <= 64
    dimension = h.shape[0]
    precision, steps, order = 2**54, 128, 18
    left, right = formation.starting_coins(route, 1)
    real = np.zeros(dimension, dtype=np.int64)
    imag = np.zeros_like(real)
    real[4*route[0]+2*left+right] = precision
    limit = np.iinfo(np.int64).max
    for _ in range(steps):
        # Bound the norm of the input to this step and certify all int64 operations.
        rr, ii = real.astype(object), imag.astype(object)
        assert rr @ rr + ii @ ii <= (2 * precision)**2
        x, y = real.copy(), imag.copy()
        for k in range(1, order+1):
            m = int(max(abs(x).max(),abs(y).max()))
            assert m * row < limit
            x, y = (h @ y) // (64*k), -(h @ x) // (64*k)
            assert int(abs(real).max()) + int(abs(x).max()) < limit
            assert int(abs(imag).max()) + int(abs(y).max()) < limit
            real += x
            imag += y
    # Unitary exact steps do not amplify the accumulated local errors. Each
    # rounded Taylor term costs at most sqrt(2d)/precision, amplified by at
    # most exp(row/64)<3 in the remaining terms of its own step.
    rounding = Fraction(3 * order * (1+isqrt(2*dimension)), precision)
    taylor = 2 * Fraction(row,64)**(order+1) / factorial(order+1)
    error = steps * (rounding + taylor)
    assert error < Fraction(1,10**9)
    re = real.reshape(-1,4).astype(object)
    im = imag.reshape(-1,4).astype(object)
    formed = (sp.Matrix(re.T@re+im.T@im) + sp.I*sp.Matrix(im.T@re-re.T@im))/precision**2
    a = sp.Matrix(2,2,lambda i,j: sum(formed[2*i+k,2*j+k] for k in range(2)))
    b = sp.Matrix(2,2,lambda i,j: sum(formed[2*k+i,2*k+j] for k in range(2)))
    # The sixth-order term of the charge-blind directed return reads this contrast.
    assert b[1,1]-b[0,0] > sp.Rational(95,100)+3*error
    product = sp.kronecker_product(a,b)
    low, high = sp.Rational(1,10**7), 8
    # With normalized true marginals, the four diagonal correlation changes
    # alternate in sign. This entry certifies a nonzero two-arm front effect.
    assert formed[0,0]-product[0,0] > sp.Rational(17,1000)+10*error
    for target in (formed-(3*error+sp.Rational(1,10**9))*sp.eye(4),
                   formed-low*product-(3+7*low)*error*sp.eye(4),
                   high*product-formed-(3+7*high)*error*sp.eye(4)):
        assert target == target.H
        assert all(target[:n,:n].det() > 0 for n in range(1,5))
    return dict(max_abs_row_sum=row, dimension=dimension, fractional_bits=54,
                steps=steps, degree=order, state_vector_error_upper=float(error),
                formed_minimum_certified=1e-9, probability_ratio_floor="1/10000000",
                probability_ratio_ceiling=8, diagonal_covariance_lower="17/1000",
                right_coin_one_minus_zero_lower="95/100")


if __name__ == "__main__":
    print(certificate())
