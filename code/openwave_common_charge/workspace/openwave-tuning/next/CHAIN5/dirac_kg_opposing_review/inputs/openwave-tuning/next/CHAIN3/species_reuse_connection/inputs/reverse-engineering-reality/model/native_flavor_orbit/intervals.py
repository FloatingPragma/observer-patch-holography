"""Rational intervals for the declared transcendental source gates."""
from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial, isqrt

import sympy as s


@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    def __init__(self, lo, hi=None):
        lo, hi = F(lo), F(lo if hi is None else hi)
        if lo > hi:
            raise ValueError('reversed rational interval')
        object.__setattr__(self, 'lo', lo)
        object.__setattr__(self, 'hi', hi)

    def __add__(self, other):
        other = box(other)
        return I(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+-box(other)

    def __rsub__(self, other):
        return box(other)+-self

    def __mul__(self, other):
        other = box(other)
        values = tuple(a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi))
        return I(min(values), max(values))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ValueError('interval denominator contains zero')
        return I(1/self.hi, 1/self.lo)

    def __truediv__(self, other):
        return self*box(other).reciprocal()

    def __pow__(self, power):
        if type(power) is not int:
            raise ValueError('integer interval powers are required')
        if power < 0:
            return self.reciprocal()**(-power)
        if power == 0:
            return I(1)
        if power % 2 and self.lo < 0:
            return I(self.lo**power, self.hi**power)
        values = (self.lo**power, self.hi**power)
        return I(0 if self.lo <= 0 <= self.hi else min(values), max(values))

    def strings(self):
        return [str(self.lo), str(self.hi)]


def box(value):
    return value if isinstance(value, I) else I(value)


def sqrt_rational(value, digits=48):
    value = F(value)
    if value < 0:
        raise ValueError('a nonnegative radicand is required')
    scale = 10**digits
    integer = isqrt(value.numerator*scale*scale//value.denominator)
    lo, hi = F(integer, scale), F(integer+1, scale)
    if lo*lo == value:
        return I(lo)
    if not lo*lo < value < hi*hi:
        raise ArithmeticError('integer square-root enclosure failed')
    return I(lo, hi)


def atan_reciprocal(n, terms=36):
    if type(n) is not int or n < 2 or terms < 1:
        raise ValueError('an alternating arctangent domain is required')
    total = sum((F((-1)**k, (2*k+1)*n**(2*k+1)) for k in range(terms)), F(0))
    other = total+F((-1)**terms, (2*terms+1)*n**(2*terms+1))
    return I(min(total, other), max(total, other))


def pi_interval():
    return rounded_out(16*atan_reciprocal(5)-4*atan_reciprocal(239), 48)


def sine(value, terms=24):
    value = box(value)
    if not 0 <= value.lo <= value.hi <= 2:
        raise ValueError('the certified sine domain is [0,2]')
    total = sum(((-1)**k*value**(2*k+1)/factorial(2*k+1)
                 for k in range(terms)), I(0))
    remainder = value.hi**(2*terms+1)/factorial(2*terms+1)
    return total+I(-remainder, remainder)


def expression(value, substitutions=None):
    value = s.sympify(value)
    substitutions = {} if substitutions is None else substitutions
    if value in substitutions:
        return box(substitutions[value])
    if value.is_Rational:
        return I(F(int(value.p), int(value.q)))
    if value == s.pi:
        return pi_interval()
    if value.func is s.Add:
        return sum((expression(x, substitutions) for x in value.args), I(0))
    if value.func is s.Mul:
        result = I(1)
        for x in value.args:
            result *= expression(x, substitutions)
        return result
    if value.func is s.Pow:
        base, exponent = value.args
        if exponent.is_Integer:
            return expression(base, substitutions)**int(exponent)
        if exponent == s.Rational(1, 2) and base.is_Rational:
            return sqrt_rational(F(int(base.p), int(base.q)))
    if value.func is s.sin:
        return rounded_out(sine(expression(value.args[0], substitutions)), 38)
    raise ValueError('expression outside the rational source enclosure grammar: '+str(value))


def matrix(value, substitutions=None):
    return [[expression(value[i, j], substitutions) for j in range(value.cols)]
            for i in range(value.rows)]


def matmul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), I(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def shifted(a, x):
    return [[value-(x if i == j else 0) for j, value in enumerate(row)]
            for i, row in enumerate(a)]


def inertia_below(a, x):
    """Interval LDL congruence; every pivot must have a certified sign."""
    b = shifted(a, x)
    negative = 0
    for k in range(len(a)):
        pivot = b[k][k]
        if pivot.lo <= 0 <= pivot.hi:
            raise ValueError('an LDL pivot has uncertified sign')
        negative += pivot.hi < 0
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                b[i][j] -= b[i][k]*b[k][j]/pivot
    return negative


def rounded_out(value, digits=12, probability=False):
    scale = 10**digits
    lower = F(value.lo.numerator*scale//value.lo.denominator, scale)
    upper = F(-((-value.hi.numerator*scale)//value.hi.denominator), scale)
    if probability:
        lower, upper = max(F(0), lower), min(F(1), upper)
    return I(lower, upper)
