"""Integer interval arithmetic for the emitted binary-angle controls.

No floating trigonometric routine supplies the synthesis certificate.
Intervals use a fixed dyadic grid and directed integer rounding.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
from math import factorial

BITS = 120
Q = 1 << BITS
ANGLE_BITS = 41
AQ = 1 << ANGLE_BITS


def ceildiv(a, b):
    return -((-a)//b)


def exact(value):
    value = F(value)
    return value.numerator*Q//value.denominator, ceildiv(value.numerator*Q, value.denominator)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def neg(a):
    return -a[1], -a[0]


def mul(a, b):
    products = [x*y for x in a for y in b]
    return min(products)//Q, ceildiv(max(products), Q)


def divide(a, positive):
    if type(positive) is not int or positive <= 0:
        raise ValueError('positive integer divisor required')
    return a[0]//positive, ceildiv(a[1], positive)


@lru_cache(maxsize=32768)
def sincos(value):
    """Enclose sin/cos of an exact rational angle with absolute value at most4."""
    value = F(value)
    if abs(value) > 4:
        raise ValueError('the certified angle interval is [-4,4]')
    x = exact(value); minus_square = neg(mul(x, x))
    cosine = term_c = exact(1)
    sine = term_s = x
    for n in range(1, 49):
        term_c = divide(mul(term_c, minus_square), (2*n-1)*(2*n))
        term_s = divide(mul(term_s, minus_square), (2*n)*(2*n+1))
        cosine = add(cosine, term_c); sine = add(sine, term_s)
    # Taylor's remainder uses |sin|,|cos|<=1. The integer checks are exact.
    if not (4**98*Q < factorial(98) and 4**99*Q < factorial(99)):
        raise ArithmeticError('the fixed Taylor remainder is not below one grid unit')
    return (sine[0]-1, sine[1]+1), (cosine[0]-1, cosine[1]+1)


def atan_reciprocal(n, terms=90):
    total = sum((F((-1)**k, (2*k+1)*n**(2*k+1)) for k in range(terms)), F())
    next_term = F((-1)**terms, (2*terms+1)*n**(2*terms+1))
    return min(total, total+next_term), max(total, total+next_term)


@lru_cache(maxsize=1)
def pi_interval():
    """Machin's identity with alternating-series rational remainders."""
    a, b = atan_reciprocal(5), atan_reciprocal(239)
    lower, upper = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
    if not F(3) < lower < upper < F(22, 7):
        raise ArithmeticError('the Machin range check failed')
    return lower, upper


@lru_cache(maxsize=1)
def pi_grid():
    lower, upper = pi_interval()
    lo = (lower*AQ+F(1, 2)).numerator//(lower*AQ+F(1, 2)).denominator
    hi = (upper*AQ+F(1, 2)).numerator//(upper*AQ+F(1, 2)).denominator
    if lo != hi:
        raise ArithmeticError('the rational enclosure does not select one pi grid point')
    return lo


@lru_cache(maxsize=4096)
def split_angle(left, right):
    """A dyadic theta approximating atan2(sqrt(right),sqrt(left)) within1/AQ.

    The proof uses monotonicity of cos(2 theta) and a global square-root
    inverse bound at endpoints. Inputs are exact nonnegative squared norms.
    """
    left, right = F(left), F(right)
    if left < 0 or right < 0:
        raise ValueError('nonnegative exact subtree weights required')
    if right == 0:
        return 0
    upper = pi_interval()[1]*AQ/2
    high = ceildiv(upper.numerator, upper.denominator)
    if left == 0:
        return high
    wanted = (left-right)/(left+right)
    low = 0
    while high-low > 1:
        mid = (low+high)//2
        _, c = sincos(F(2*mid, AQ))
        if F(c[0], Q) > wanted:
            low = mid
        elif F(c[1], Q) < wanted:
            high = mid
        else:
            # |cos(2theta)-cos(2target)|<=width. On [0,pi/2],
            # |theta-target|<=2 sqrt(width), so this is <=1/AQ.
            if 4*(c[1]-c[0])*AQ*AQ > Q:
                raise ArithmeticError('trigonometric enclosure too wide for the angle grid')
            return mid
    return low


def phase_angle(real, imag):
    """Dyadic argument of a nonzero Gaussian rational, error at most2/AQ."""
    real, imag = F(real), F(imag)
    if real == imag == 0:
        return 0
    angle = split_angle(real*real, imag*imag)
    if real < 0:
        angle = pi_grid()-angle
    return -angle if imag < 0 else angle


def binary_exact(numerator, denominator):
    value = F(numerator, denominator)
    result = float(value)
    if F(result) != value:
        raise ArithmeticError('the synthesized angle is not exactly binary64 representable')
    return result


def scalar_error(value, interval):
    """Exact upper bound for a stored real number versus an enclosed target."""
    value = F(float(value))
    return max(abs(value-F(interval[0], Q)), abs(value-F(interval[1], Q)))
