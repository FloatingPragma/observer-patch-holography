"""Read and certify a fixed size-three frame; no runtime eigensolver selection."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import math
import sys

sys.dont_write_bytecode = True
import numpy as np
import providers

HERE = Path(__file__).resolve().parent
RANK = 646
DENOMINATOR = 1 << 52


def exact_product(left, right):
    """Integer GEMM through binary64 with every product/partial sum exact.

    All operands are integers of magnitude below2^53, and the absolute
    dot-product sum is strictly below2^53. Every product and every addition
    is therefore an exactly representable binary64 integer, in any order.
    """
    assert left.shape[1] == right.shape[0]
    maximum_left, maximum_right = int(np.max(abs(left))), int(np.max(abs(right)))
    assert left.shape[1]*maximum_left*maximum_right < 1 << 53
    assert max(maximum_left, maximum_right) < 1 << 53
    assert np.finfo(np.float64).nmant == 52
    product = left.astype(np.float64)@right.astype(np.float64)
    assert np.array_equal(product, np.rint(product))
    return product.astype(np.int64).astype(object)


def sqrt_upper(value):
    value = F(value); n = value.numerator*value.denominator
    r = math.isqrt(n)
    return F(r+(r*r < n), value.denominator)


def sqrt_lower(value):
    value = F(value)
    return F(math.isqrt(value.numerator*value.denominator), value.denominator)


@lru_cache(None)
def geometry():
    data = json.loads((HERE/'GEOMETRY.json').read_text())
    edges = data['edges']; B = np.zeros((324, 972), np.int64)
    for e, (a, b) in enumerate(edges): B[a, e] = 1; B[b, e] = -1
    A = np.zeros((972, 702), np.int64)
    for e, f, value in data['face_entries']: A[e, f] = value
    periods = np.zeros((3, 972), np.int64)
    for axis, edge, value in data['period_times3_entries']: periods[axis, edge] = value
    winding = np.array(data['winding_columns'], np.int64).T
    assert not np.any(B@A) and not np.any(periods@A)
    assert not np.any(B@winding) and np.array_equal(periods@winding, 3*np.eye(3, dtype=int))
    return B, A, A@A.T, periods, winding


def certificate_from_data(frame_data):
    K = geometry()[2]; den = frame_data['denominator']
    Q = np.array(frame_data['frame_integers'], np.int64)
    lam = np.array(frame_data['curvature_numerators'], np.int64)
    assert den == DENOMINATOR and Q.shape == (972, RANK) and len(lam) == RANK
    assert min(lam) > 0 and all(a <= b for a, b in zip(lam, lam[1:]))
    assert int(np.max(np.sum(abs(K), axis=1)))*int(np.max(abs(Q))) < 1 << 63
    # Three signed 18-bit limbs; exact reconstruction after guarded int64 sums.
    base = 1 << 18; limbs = []; tail = Q.copy()
    for _ in range(2):
        high = tail//base; limbs.append(tail-high*base); tail = high
    limbs.append(tail)
    gram = np.zeros((RANK, RANK), dtype=object)
    for i, left in enumerate(limbs):
        for j, right in enumerate(limbs):
            assert 972*int(np.max(abs(left)))*int(np.max(abs(right))) < 1 << 53
            gram += exact_product(left.T, right)*base**(i+j)
    gram -= den**2*np.eye(RANK, dtype=object)
    KQ = sum((exact_product(K, limb)*base**i for i, limb in enumerate(limbs)), np.zeros(Q.shape, dtype=object))
    residual = KQ*den-Q.astype(object)*lam.astype(object)[None, :]
    gram2 = F(sum(int(x)**2 for x in gram.flat), den**4)
    residual2 = F(sum(int(x)**2 for x in residual.flat), den**4)
    epsilon, r = sqrt_upper(gram2), sqrt_upper(residual2)
    minimum, maximum = F(int(lam[0]), den), F(int(lam[-1]), den)
    shift = r/minimum; gamma = epsilon+2*(1+epsilon)*shift+shift**2
    eta = r+(F(1504, 17)+maximum)*(shift+gamma)
    assert gamma < F(1, 2) and eta < minimum/2
    frequency = F(5, 4)*eta/sqrt_lower(minimum)
    return {'denominator':den, 'modes':RANK, 'Gram_squared_Frobenius_exact':str(gram2),
        'residual_squared_Frobenius_exact':str(residual2), 'Gram_error_upper':str(epsilon),
        'residual_upper':str(r), 'polar_Gram_error_upper':str(gamma),
        'curvature_operator_error_upper':str(eta), 'per_frequency_error_upper':str(frequency),
        'cap1_H_error_including_vacuum':str((RANK//2+1)*frequency),
        'minimum_curvature':str(minimum), 'maximum_curvature':str(maximum),
        'exact_frame':'U=polar(K Q Lambda^-1); full coexact range, BU=0 exactly',
        'integer_product_domain':'18-bit integer limbs, exact binary64 products and partial sums below2^53; arbitrary-precision recombination',
        'harmonic_modes':3, 'physical_cap1_states':647, 'padding_states':377}


@lru_cache(None)
def data():
    c = json.loads((HERE/'CERTIFICATE.json').read_text())
    assert providers.sha(HERE/'FRAME.json') == c['frame_sha256']
    assert providers.sha(HERE/'GEOMETRY.json') == c['geometry_sha256']
    return json.loads((HERE/'FRAME.json').read_text()), c


def spectrum():
    d, c = data()
    return [F(n, d['denominator']) for n in d['curvature_numerators']], c
