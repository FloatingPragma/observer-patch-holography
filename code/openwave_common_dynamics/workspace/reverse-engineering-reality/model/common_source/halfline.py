"""Exact infinite coordinate-sign readout of finite Hermite polynomials.

A vector maps (n0, n1, gauge) to Gaussian-rational coefficients of the
UNNORMALIZED Bargmann basis (a0.dagger)^n0 (a1.dagger)^n1 |0,0,gauge>.
Its norm therefore includes n0! n1!. The observable is multiplication by
1_(x0 > 0) in the infinite oscillator, not sign of a finite Q matrix.

The complete postread branches are symbolically P_b |vector>: in coordinate
space multiply its Hermite-polynomial Gaussian wavefunction by 1_(x0 > 0)
or 1_(x0 < 0). These are infinite-dimensional vectors, with every surviving
coherence intact. The coefficient formulas below compute their exact norms;
they do not project those outputs back into finite occupation, replace them
by scalar probabilities, condition on success, or normalize a branch.
"""
from __future__ import annotations

from fractions import Fraction as Q
from functools import lru_cache
import math


MAX_OCCUPATION = 512
ARCTAN_TERMS = 64
ROOT_SCALE = 10**60
Gaussian = tuple[Q, Q]
Vector = dict[tuple[int, int, int], Gaussian]


def _occupation(value):
    if type(value) is not int or not 0 <= value <= MAX_OCCUPATION:
        raise ValueError(f"an integer occupation from 0 through {MAX_OCCUPATION} required")
    return value


def _vector(value):
    if type(value) is not dict:
        raise ValueError("a finite Gaussian-rational Bargmann coefficient dictionary required")
    for key, coefficient in value.items():
        if (type(key) is not tuple or len(key) != 3
                or any(type(n) is not int or not 0 <= n <= MAX_OCCUPATION for n in key[:2])
                or type(key[2]) is not int or not 0 <= key[2] < 4
                or type(coefficient) is not tuple or len(coefficient) != 2
                or any(type(part) is not Q for part in coefficient)):
            raise ValueError("bounded occupations, gauge 0..3 and exact Gaussian rationals required")
    return value


def _normalization(value):
    if type(value) is not Q or value <= 0:
        raise ValueError("a positive exact rational initial normalization required")
    return value


@lru_cache(maxsize=MAX_OCCUPATION + 1)
def _factorial(value):
    return math.factorial(value)


def norm_squared(vector: Vector) -> Q:
    """Exact norm of the supplied unnormalized polynomial, including factorials."""
    return sum((_factorial(n0) * _factorial(n1) * (real*real + imag*imag)
                for (n0, n1, _), (real, imag) in _vector(vector).items()), Q(0))


@lru_cache(maxsize=65536)
def _coefficient(m, n):
    if m == n:
        return Q(_factorial(m), 2), Q(0)
    if (m - n) % 2 == 0:
        return Q(0), Q(0)
    if m % 2:
        m, n = n, m
    j, k = m // 2, n // 2
    cross = Q((-1)**(j+k) * _factorial(m) * _factorial(n),
              2**(j+k) * _factorial(j) * _factorial(k) * (n-m))
    return Q(0), cross


def halfline_coefficient(m: int, n: int) -> tuple[Q, Q]:
    """Return (r,s) with <m_B|P_positive|n_B> = r+s/sqrt(2*pi).

    In the normalized Hermite basis, the even/odd entry is
    (-1)^(j+k) sqrt((2j)! (2k+1)!)/
    (sqrt(2*pi) 2^(j+k) j! k! (2k+1-2j)). Multiplication by sqrt(m! n!)
    for the Bargmann basis leaves the rational coefficient returned here.
    """
    return _coefficient(_occupation(m), _occupation(n))


def positive_weight_coefficients(vector: Vector, normalization=Q(1)) -> tuple[Q, Q]:
    """Exact r,s for the positive branch weight r+s/sqrt(2*pi).

    ``normalization`` is the squared norm of an initial polynomial used to
    specify a normalized preparation. It is explicit, never inferred from
    the current vector; an unnormalized retained branch stays unnormalized.
    Spectator occupation and gauge labels are orthogonal untouched factors.
    """
    vector, normalization = _vector(vector), _normalization(normalization)
    groups = {}
    for (n0, n1, gauge), coefficient in vector.items():
        groups.setdefault((n1, gauge), {})[n0] = coefficient
    cross = Q(0)
    for (n1, _), coefficients in groups.items():
        even = [(n, z) for n, z in coefficients.items() if n % 2 == 0]
        odd = [(n, z) for n, z in coefficients.items() if n % 2]
        for m, (a, b) in even:
            for n, (c, d) in odd:
                cross += 2 * _factorial(n1) * _coefficient(m, n)[1] * (a*c + b*d)
    return norm_squared(vector) / (2*normalization), cross / normalization


def _arctan_interval(denominator):
    """Alternating Taylor enclosure of atan(1/denominator), exactly rational."""
    total = sum((Q((-1)**k, (2*k+1) * denominator**(2*k+1))
                 for k in range(ARCTAN_TERMS)), Q(0))
    next_term = Q(1, (2*ARCTAN_TERMS+1) * denominator**(2*ARCTAN_TERMS+1))
    return ((total-next_term, total) if ARCTAN_TERMS % 2
            else (total, total+next_term))


@lru_cache(maxsize=1)
def pi_interval() -> tuple[Q, Q]:
    """Machin's pi = 16 atan(1/5)-4 atan(1/239), with signed tail bounds."""
    a, b = _arctan_interval(5)
    c, d = _arctan_interval(239)
    return 16*a-4*d, 16*b-4*c


def _sqrt_lower(value):
    return Q(math.isqrt(value.numerator * ROOT_SCALE**2 // value.denominator), ROOT_SCALE)


def _sqrt_upper(value):
    lower = _sqrt_lower(value)
    return lower if lower*lower == value else lower + Q(1, ROOT_SCALE)


@lru_cache(maxsize=1)
def inverse_sqrt_two_pi_interval() -> tuple[Q, Q]:
    """Rational enclosure from Machin tails and integer square roots only."""
    low, high = pi_interval()
    return _sqrt_lower(1/(2*high)), _sqrt_upper(1/(2*low))


def weight_interval(vector: Vector, outcome: int, normalization=Q(1)) -> tuple[Q, Q]:
    """Enclose one infinite sign branch, preserving its raw probability mass.

    Outcome one is x0>0 and outcome zero is x0<0. Intersecting with the
    proved interval [0, norm_squared(vector)/normalization] is safe even
    when the supplied polynomial has mass greater than one.
    """
    if type(outcome) is not int or outcome not in (0, 1):
        raise ValueError("a binary integer sign outcome required")
    constant, cross = positive_weight_coefficients(vector, normalization)
    cross *= 2*outcome-1
    low, high = inverse_sqrt_two_pi_interval()
    if cross < 0:
        low, high = high, low
    return max(Q(0), constant+cross*low), min(2*constant, constant+cross*high)


def positive_weight_interval(vector: Vector, normalization=Q(1)) -> tuple[Q, Q]:
    """Enclose the unnormalized positive-half-line weight of this polynomial."""
    return weight_interval(vector, 1, normalization)
