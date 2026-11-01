"""Computed static response, admitted only by exact residual arithmetic.

A floating-point linear solve proposes an integer matrix X. Its accuracy
is never trusted: all residuals, square-root enclosures and accepted error
bounds below use integers and rational arithmetic. No proposed matrix hash
is fixed, so harmless platform differences can pass the same certificate.
"""
from fractions import Fraction as F
from math import isqrt

import numpy as np
from scipy.linalg import solve

from ..interacting_vacuum import ground
from . import response as r


BITS = 40
SCALE = 1 << BITS
SQRT_SCALE = 1 << 80
SQRT_LOWER = isqrt(3*SQRT_SCALE**2)
C_ERROR = F(1, 10**7)
INVERSE_ERROR = F(1, 10**6)
C_DIAGONAL, C_ADJACENCY = F('12.9058381883'), F('-4.9026412320')
INVERSE_DIAGONAL, INVERSE_ADJACENCY = F('.1366376936'), F('.0519056245')


def quadratic_frobenius_bound(a, b, denominator):
    """Outward norm bound of (a+sqrt3 b)/denominator, exact entry intervals."""
    if type(denominator) is not int or denominator <= 0:
        raise ValueError('positive integer denominator required')
    if np.shape(a) != np.shape(b): raise ValueError('matching numerator shapes required')
    square = 0
    for x, y in zip(np.asarray(a).flat, np.asarray(b).flat, strict=True):
        if not isinstance(x, (int,np.integer)) or not isinstance(y, (int,np.integer)):
            raise ValueError('integer numerators required')
        x, y = int(x), int(y)
        lo, hi = x*SQRT_SCALE+y*SQRT_LOWER, x*SQRT_SCALE+y*(SQRT_LOWER+1)
        bound = max(abs(lo), abs(hi))
        square += bound*bound
    root = isqrt(square)
    if root*root != square: root += 1
    return F(root, denominator*SQRT_SCALE)


def propose():
    data = r.dilation('ground','down')
    v = (data.vectors[0]+np.sqrt(3)*data.vectors[1])/36
    k = data.h16.toarray()/16+np.diag(data.sqrt3_shift*np.sqrt(3)/16)
    value = solve(k, v, assume_a='pos')*SCALE
    if not np.all(np.isfinite(value)) or np.max(np.abs(value)) >= 2**60:
        raise ArithmeticError('untrusted proposal outside integer range')
    return np.rint(value).astype(np.int64)


def verify(candidate):
    """Verify ||V†K^-1V-(c0 I+c1 A)||<1e-7, then its inverse."""
    value = np.asarray(candidate)
    if value.shape != (560,8) or value.dtype != np.int64 or r._maxabs(value) > 2**50:
        raise ValueError('bounded560x8 integer proposal required')
    data = r.dilation('ground','down'); a,b = data.vectors
    h = data.h16
    row_bound = int(np.max(np.asarray(abs(h).sum(axis=1))))
    # Keep each actual int64 operation exact, including36 times H16 X.
    if 36*(row_bound+12)*r._maxabs(value)+16*SCALE*max(r._maxabs(a),r._maxabs(b)) > r.LIMIT:
        raise ArithmeticError('integer residual would overflow')
    ra = 16*SCALE*a-36*(h@value)
    rb = 16*SCALE*b-36*data.sqrt3_shift[:,None]*value
    residual = quadratic_frobenius_bound(ra,rb,576*SCALE)
    inverse_vector_error = residual/r.GAP
    # V†X is known overQsqrt3. Compare with the fixed rational matrix C0.
    ta,tb = r._matmul(a.T,value),r._matmul(b.T,value)
    adjacency=ground.band_data()[0]
    fixed=(int(C_DIAGONAL*10**10)*np.eye(8,dtype=object)
           +int(C_ADJACENCY*10**10)*adjacency.astype(object))
    da=ta.astype(object)*10**10-fixed*(36*SCALE)
    db=tb.astype(object)*10**10
    compression_error=quadratic_frobenius_bound(da,db,36*SCALE*10**10)
    if inverse_vector_error+compression_error >= C_ERROR:
        raise ArithmeticError('static response proposal failed its exact residual budget')
    # C>=I/2 because0<K<2. Thus||C^-1||<=2. The proposed rational inverse
    # Y=y0I+y1A has exact product with C0 because A²=3I.
    scalar=1-C_DIAGONAL*INVERSE_DIAGONAL-3*C_ADJACENCY*INVERSE_ADJACENCY
    edge=-(C_DIAGONAL*INVERSE_ADJACENCY+C_ADJACENCY*INVERSE_DIAGONAL)
    inverse_residual=abs(scalar)+2*abs(edge)
    inverse_norm=abs(INVERSE_DIAGONAL)+2*abs(INVERSE_ADJACENCY)
    if 2*(inverse_residual+C_ERROR*inverse_norm) >= INVERSE_ERROR:
        raise ArithmeticError('fixed inverse response failed its exact perturbation budget')
    return {'static_positive_response':'C=V^*K^(-1)V=-G_down(0)',
            'C_approximation':{'I':str(C_DIAGONAL),'A':str(C_ADJACENCY)},
            'C_operator_error_strict':str(C_ERROR),
            'inverse_C_approximation':{'I':str(INVERSE_DIAGONAL),'A':str(INVERSE_ADJACENCY)},
            'inverse_C_operator_error_strict':str(INVERSE_ERROR),
            'static_memory_approximation':{'I':str(INVERSE_DIAGONAL-F(1,2)),
                                           'A':str(INVERSE_ADJACENCY-F(1,16))},
            'static_memory_operator_error_strict':str(INVERSE_ERROR),
            'proof':'Exact residual ofKX=V gives||X-K^-1V||<=||residual||/delta. V is an isometry. Fixed rational C0 and Y are checked by outward Qsqrt3 Frobenius bounds and A²=3I; C>=I/2 gives||C^-1||<=2. Sigma(0)=C^-1-M1.',
            'proposal_fixed_point_bits':BITS,
            'floating_point_solve_admits_the_result':False,
            'exact_equality_to_affine_adjacency_asserted':False,
            'ground_down_channel_only':True}


def certificate():
    return verify(propose())
