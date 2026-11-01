"""Exact dyadic matrix powers with a conditional unitary error theorem.

The caller supplies a proved unitary target and entry enclosures for it.
This module never infers unitarity from a floating residual.  Every matrix
product is a Gaussian-integer product followed by nearest dyadic rounding.
If A,B approximate unitary targets with operator errors a,b, the rounded
product error is at most a+b+ab+d*2**(-bits).  Bounds themselves are rounded
outward to the same dyadic lattice, keeping fractional precision fixed.
For controlled near-unitary inputs this avoids growing integer sizes even
when the exponent contains millions of steps; no such size bound is
claimed for arbitrary inputs whose error bound grows without control.
"""
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json

import numpy as np

from ..certified_reference import interval as iv


def _exact(value):
    if type(value) not in (int, F):
        raise ValueError('exact integer or Fraction required')
    return F(value)


def _bits(bits):
    if type(bits) is not int or not 8 <= bits <= 512:
        raise ValueError('integer precision between8 and512 required')
    return bits


def _ceil_dyadic(value, bits):
    value = _exact(value)
    scale = 1 << bits
    return F(-((-value.numerator*scale)//value.denominator), scale)


def _nearest_fraction(value, bits):
    value = _exact(value)
    numerator = value.numerator*(1 << bits)
    # Ties go toward positive infinity; either choice has error at most1/2 ulp.
    return (2*numerator+value.denominator)//(2*value.denominator)


def _integers(array):
    result = np.array(array, dtype=object, copy=True)
    if result.ndim != 2 or result.shape[0] != result.shape[1] or not result.shape[0]:
        raise ValueError('nonempty square matrix required')
    if any(isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, np.integer))
           for x in result.flat):
        raise ValueError('integer matrix numerators required; no floating input')
    result = np.array([int(x) for x in result.flat], dtype=object).reshape(result.shape)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class CertifiedUnitary:
    """Dyadic approximation, conditional on the caller's unitary-target proof."""

    real: np.ndarray
    imag: np.ndarray
    bits: int
    error: F
    multiplications: int = 0

    def __post_init__(self):
        _bits(self.bits)
        real, imag = _integers(self.real), _integers(self.imag)
        error = _exact(self.error)
        if real.shape != imag.shape or error < 0:
            raise ValueError('matching matrices and nonnegative error required')
        if type(self.multiplications) is not int or self.multiplications < 0:
            raise ValueError('nonnegative integer multiplication count required')
        object.__setattr__(self, 'real', real)
        object.__setattr__(self, 'imag', imag)
        object.__setattr__(self, 'error', error)

    @property
    def dimension(self):
        return self.real.shape[0]

    def diagnostic_matrix(self):
        """Floating display only; never used by a certificate calculation."""
        return (np.asarray(self.real, dtype=float)+1j*np.asarray(self.imag, dtype=float))/float(1 << self.bits)

    def certificate(self):
        payload = {'bits': self.bits, 'real': self.real.tolist(), 'imag': self.imag.tolist()}
        return {'dimension': self.dimension, 'precision_bits': self.bits,
                'operator_error_upper': str(self.error),
                'rounded_matrix_multiplications': self.multiplications,
                'matrix_sha256': hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
                'scope': 'operator norm relative to the supplied proved unitary target; no host matrix arithmetic or residual admission'}


def _endpoints(entry):
    if type(entry) is iv.Interval:
        return F(entry.lo, iv.SCALE), F(entry.hi, iv.SCALE)
    if type(entry) in (int, F):
        value = _exact(entry)
        return value, value
    if type(entry) in (tuple, list) and len(entry) == 2:
        low, high = map(_exact, entry)
        if low <= high:
            return low, high
    raise ValueError('ordered exact rational endpoints or a certified interval required')


def from_intervals(real, imag, *, bits=128):
    """Enclose a separately proved unitary target and quantize its center.

    Inputs are square lists of exact rationals, (lower,upper) rational pairs,
    or certified_reference.interval.Interval objects.  The maximum complex
    entry error is bounded by |real error|+|imag error|.  Frobenius norm is
    then at most dimension times that maximum; operator norm is no larger.
    """
    _bits(bits)
    try:
        d = len(real)
        if not d or len(imag) != d or any(len(row) != d for row in list(real)+list(imag)):
            raise ValueError('matching nonempty square interval grids required')
    except TypeError as exc:
        raise ValueError('square interval grids required') from exc
    re, im = np.zeros((d, d), dtype=object), np.zeros((d, d), dtype=object)
    maximum = F()
    for i in range(d):
        for j in range(d):
            errors = []
            for target, grid in ((re, real), (im, imag)):
                low, high = _endpoints(grid[i][j])
                numerator = _nearest_fraction((low+high)/2, bits)
                target[i, j] = numerator
                center = F(numerator, 1 << bits)
                errors.append(max(abs(center-low), abs(center-high)))
            maximum = max(maximum, sum(errors, F()))
    return CertifiedUnitary(re, im, bits, _ceil_dyadic(d*maximum, bits))


def identity(dimension, *, bits=128):
    _bits(bits)
    if type(dimension) is not int or dimension < 1:
        raise ValueError('positive integer dimension required')
    return CertifiedUnitary(np.eye(dimension, dtype=object)*(1 << bits),
                            np.zeros((dimension, dimension), dtype=object), bits, F())


def multiply(left, right):
    """Return a certified approximation to the ordered unitary product UV."""
    if not isinstance(left, CertifiedUnitary) or not isinstance(right, CertifiedUnitary):
        raise ValueError('certified dyadic inputs required')
    if left.bits != right.bits or left.dimension != right.dimension:
        raise ValueError('equal dimensions and dyadic precisions required')
    # dtype=object multiplication is exact Python-integer arithmetic.
    real = left.real@right.real-left.imag@right.imag
    imag = left.real@right.imag+left.imag@right.real
    scale = 1 << left.bits
    rounded = lambda array: (2*array+scale)//(2*scale)
    # Each of2d^2 real components has error<=2^(-bits-1), hence
    # ||rounding||_op <= ||rounding||_F <= d*2^(-bits)/sqrt2 <= d*2^(-bits).
    rounding = F(left.dimension, scale)
    error = _ceil_dyadic(left.error+right.error+left.error*right.error+rounding, left.bits)
    return CertifiedUnitary(rounded(real), rounded(imag), left.bits, error,
                            left.multiplications+right.multiplications+1)


def power(matrix, exponent):
    """Binary power; the error tracks target U**exponent, not a BLAS result.

    The multiplication count reports operations actually performed in this
    call (plus construction of the supplied base), not the expanded tree.
    """
    if not isinstance(matrix, CertifiedUnitary):
        raise ValueError('certified dyadic matrix required')
    if type(exponent) is not int or exponent < 0:
        raise ValueError('nonnegative integer exponent required')
    result = identity(matrix.dimension, bits=matrix.bits)
    base = matrix
    remaining, count = exponent, matrix.multiplications
    while remaining:
        if remaining & 1:
            result = multiply(result, base)
            count += 1
        remaining >>= 1
        if remaining:
            base = multiply(base, base)
            count += 1
    return CertifiedUnitary(result.real, result.imag, result.bits, result.error, count)


def expectation(matrix, vector_real, vector_imag=None):
    """Exact center/error of <v,Uv>/<v,v> for a nonzero integer vector.

    Normalization is evaluated exactly, avoiding a numerical square root.
    Both real and imaginary ideal values lie within `operator_error_upper`
    of their respective reported rational centers.
    """
    if not isinstance(matrix, CertifiedUnitary):
        raise ValueError('certified dyadic matrix required')
    vr = np.asarray(vector_real, dtype=object)
    vi = np.zeros(vr.shape, dtype=object) if vector_imag is None else np.asarray(vector_imag, dtype=object)
    if vr.shape != (matrix.dimension,) or vi.shape != vr.shape:
        raise ValueError('matching vectors of the matrix dimension required')
    if any(isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, np.integer)) for x in list(vr)+list(vi)):
        raise ValueError('exact integer vector required')
    vr = np.array(list(map(int, vr)), dtype=object)
    vi = np.array(list(map(int, vi)), dtype=object)
    norm = int(vr@vr+vi@vi)
    if not norm:
        raise ValueError('nonzero vector required')
    yr = matrix.real@vr-matrix.imag@vi
    yi = matrix.real@vi+matrix.imag@vr
    re = F(int(vr@yr+vi@yi), norm*(1 << matrix.bits))
    im = F(int(vr@yi-vi@yr), norm*(1 << matrix.bits))
    return {'real_center': str(re), 'imag_center': str(im),
            'operator_error_upper': str(matrix.error), 'integer_vector_norm_squared': norm}
