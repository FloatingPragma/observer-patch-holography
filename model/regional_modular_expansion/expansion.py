"""Exact noncommutative regional cumulants and a rational logarithm tail.

The fixed plaquette is the existing complete-Fock Weyl/contact action. No
particle-number sector, boundary coupling or eigenvalue is removed. The
general bound requires a unital CP partial trace, B >= ||H|| and
abs(beta)*B < log(2); the exact rational checker admits nonnegative beta
and B with 0 <= beta*B <= 1/2.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import math

import numpy as np
import sympy as s
from scipy import sparse

from ..field_stress.energy import GaussianInteger as GI
from ..field_symmetry.symmetry import operators
from ..full_fock_field.field import FieldGeometry


ORDER = 8
BETA = F(1, 32)
GEOMETRY = FieldGeometry((2, 2, 1))
KEPT = (0, 1, 2, 3)


def rational(value):
    if type(value) not in (int, F):
        raise ValueError('an original exact integer or Fraction is required')
    return F(value)


def _sp(value):
    value = rational(value)
    return s.Rational(value.numerator, value.denominator)


def _f(value):
    value = s.Rational(value)
    return F(int(value.p), int(value.q))


def centered(matrix):
    return matrix-s.eye(matrix.rows)*s.trace(matrix)/matrix.rows


def exact_trace(matrix, modes=8, kept=KEPT):
    """Normalized signed CAR trace, independently routed in integer arithmetic.

    The matrix numerator is a complete 2**modes square Gaussian integer.
    The discarded dimension is divided only after summing every basis state.
    """
    if (type(modes) is not int or not 1 <= modes <= 8 or type(kept) is not tuple
            or not kept or tuple(sorted(set(kept))) != kept
            or any(type(j) is not int or not 0 <= j < modes for j in kept)):
        raise ValueError('a bounded complete CAR carrier and ordered kept modes are required')
    if matrix.real.shape != (1 << modes,)*2:
        raise ValueError('the complete Fock matrix is required')
    outside = tuple(j for j in range(modes) if j not in kept)
    labels = []
    for mask in range(1 << modes):
        local = sum(((mask >> j) & 1) << k for k, j in enumerate(kept))
        env = sum(((mask >> j) & 1) << k for k, j in enumerate(outside))
        inversions = sum((mask >> a & 1)*(mask >> b & 1)
                         for a in kept for b in outside if b < a)
        labels.append((local, env, (-1)**inversions))
    result = s.zeros(1 << len(kept))
    for part, factor in ((matrix.real, 1), (matrix.imag, s.I)):
        entries = part.tocoo()
        for i, j, value in zip(entries.row, entries.col, entries.data, strict=True):
            a, env_a, sign_a = labels[int(i)]
            b, env_b, sign_b = labels[int(j)]
            if env_a == env_b:
                result[a, b] += factor*int(value)*sign_a*sign_b
    return s.ImmutableMatrix(result/(1 << len(outside)))


def log_coefficients(moments):
    """Coefficients of -log(E exp(-beta H)), in their noncommuting order."""
    order = len(moments)
    if not 1 <= order <= ORDER or any(m.shape != moments[0].shape for m in moments):
        raise ValueError('one through eight common-shaped exact moments required')
    dimension = moments[0].rows
    a = [s.zeros(dimension)]+[(-1)**n*moments[n-1]/math.factorial(n)
                            for n in range(1, order+1)]
    power = [s.eye(dimension)]+[s.zeros(dimension) for _ in range(order)]
    result = [s.zeros(dimension) for _ in range(order+1)]
    for count in range(1, order+1):
        power = [s.zeros(dimension)]+[
            sum((power[n-j]*a[j] for j in range(1, n+1)), s.zeros(dimension))
            for n in range(1, order+1)]
        for n in range(count, order+1):
            result[n] += s.Rational((-1)**count, count)*power[n]
    if any(matrix != matrix.conjugate().T for matrix in result):
        raise ArithmeticError('regional logarithm coefficient is not Hermitian')
    return tuple(s.ImmutableMatrix(matrix) for matrix in result[1:])


def majorant_coefficients(order=ORDER):
    """Positive exact coefficients of f(x)=-log(2-exp(x))."""
    if type(order) is not int or not 1 <= order <= ORDER:
        raise ValueError('integer order from one through eight required')
    moments = tuple(s.Matrix([[(-1)**n]]) for n in range(1, order+1))
    # The matrix helper computes -log(1+sum a_n x^n). Here the desired
    # series is -log(1-(exp(x)-1)), so use moments -(-1)**n instead.
    moments = tuple(-m for m in moments)
    return tuple(_f(m[0, 0]) for m in log_coefficients(moments))


def _ceil(value, denominator=10**22):
    return F(-(-value.numerator*denominator//value.denominator), denominator)


def logarithm_tail(beta, bound, order=ORDER):
    """Outward rational remainder; neither exp nor log is evaluated in floats.

    exp(x) is bounded by 32 terms and a geometric tail. Then -log(1-q)
    uses 128 terms and its geometric tail. The positive coefficients through
    the requested order are subtracted exactly. Upward rational rounding is
    charged before each composition and in the returned bound.
    """
    beta, bound = rational(beta), rational(bound)
    if beta < 0 or bound < 0:
        raise ValueError('nonnegative beta and operator norm bound required')
    x = beta*bound
    if not 0 <= x <= F(1, 2):
        raise ValueError('the rational certificate requires beta*B <= 1/2')
    coefficients = majorant_coefficients(order)
    if x == 0:
        return F(0)
    terms, term = F(1), F(1)
    for n in range(1, 33):
        term *= x/n
        terms += term
    exp_upper = _ceil(terms+term*x/33/(1-x/34))
    q = exp_upper-1
    if not 0 < q < 1:
        raise ArithmeticError('the analytic logarithm disk was lost')
    log_upper = sum((q**n/n for n in range(1, 129)), F())
    log_upper += q**129/(129*(1-q))
    polynomial = sum((a*x**n for n, a in enumerate(coefficients, 1)), F())
    tail = _ceil(log_upper-polynomial)
    if tail <= 0:
        raise ArithmeticError('outward positive logarithm remainder failed')
    return tail


@lru_cache(maxsize=4)
def plaquette(interaction=1, boundary=True):
    """Exact coefficients on all 256 field states, for the retained half square.

    U=0 and removal of the two crossing bonds are explicit negative controls.
    Only a scalar trace shift is used; no chemical-potential term is added.
    """
    if type(interaction) is not int or interaction not in (0, 1) or type(boundary) is not bool:
        raise ValueError('U=0/1 and an explicit boundary Boolean are required')
    ops = operators(FieldGeometry((2, 2, 1), interaction=float(interaction)))
    scale = ops.scale
    identity = GI.from_exact(sparse.eye(256, format='csr'))
    h = ops.H-GI.from_exact(identity.matrix()*scale*interaction)
    if not boundary:
        # Literal source bond order (0,2),(0,1),(1,3),(2,3): crossing x bonds.
        h = h-ops.bonds[0]-ops.bonds[2]
    if not (h-h.adjoint()).is_zero():
        raise ArithmeticError('actual centered field action is not Hermitian')
    # Every entry is real or imaginary: its absolute value is an integer.
    if h.real.multiply(h.imag).nnz:
        raise ArithmeticError('exact row-norm certificate requires axis-pure entries')
    rows = np.asarray((abs(h.real)+abs(h.imag)).sum(axis=1)).ravel()
    norm_bound = F(int(max(rows)), scale)
    moments, power = [], identity
    for n in range(1, ORDER+1):
        power = power@h
        moments.append(s.ImmutableMatrix(exact_trace(power)/scale**n))
    coefficients = log_coefficients(tuple(moments))
    return {'H_numerator': h, 'scale': scale, 'bound': norm_bound,
            'moments': tuple(moments), 'coefficients': coefficients,
            'interaction': interaction, 'boundary': boundary}


def polynomial(data, beta=BETA, order=ORDER):
    beta = rational(beta)
    if type(order) is not int or not 1 <= order <= ORDER:
        raise ValueError('integer order from one through eight required')
    logarithm_tail(beta, data['bound'], order)
    return s.ImmutableMatrix(sum((_sp(beta)**n*c for n, c in enumerate(
        data['coefficients'][:order], 1)), s.zeros(16)))


def pauli_terms(matrix):
    """Exact four-mode tensor coefficients, highest occupation bit first."""
    if matrix.shape != (16, 16):
        raise ValueError('the complete four-mode regional matrix is required')
    sigma = (s.eye(2), s.Matrix([[0, 1], [1, 0]]),
             s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    result = {}
    for labels in product(range(4), repeat=4):
        operator = s.kronecker_product(*(sigma[j] for j in labels))
        coefficient = s.trace(operator*matrix)/16
        if coefficient:
            result[''.join('IXYZ'[j] for j in labels)] = str(coefficient)
    return result


def witness(data=None, beta=BETA):
    """Exact full-versus-empty regional witness rejects bare beta*H_A+scalar."""
    data = plaquette() if data is None else data
    beta = rational(beta)
    correction = polynomial(data, beta)-_sp(beta)*data['coefficients'][0]
    difference = _f(correction[15, 15]-correction[0, 0])
    remainder = logarithm_tail(beta, data['bound'])
    return {'polynomial_diagonal_difference': difference,
            'exact_difference_lower': difference-2*remainder,
            'exact_difference_upper': difference+2*remainder,
            'bare_generator_distance_lower': max(F(0), abs(difference)/2-remainder),
            'modular_operator_error_upper': remainder}
