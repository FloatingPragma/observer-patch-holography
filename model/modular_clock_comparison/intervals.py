"""Small rational interval kernel; no floating decision enters a certificate.

Every operation rounds outwards to a rational decimal grid. Exponential
and logarithm enclosures use explicitly bounded positive series tails.
This is a bounded certificate kernel, not a general interval library.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt

SCALE = 10**30


def _q(x):
    if type(x) not in (int, F):
        raise ValueError('exact integer or rational endpoint required')
    return F(x)


def _floor(x):
    return F((x*SCALE).numerator//(x*SCALE).denominator, SCALE)


def _ceil(x):
    return -_floor(-x)


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __post_init__(self):
        lo, hi = _q(self.lo), _q(self.hi)
        if lo > hi:
            raise ValueError('ordered exact endpoints required')
        object.__setattr__(self, 'lo', _floor(lo))
        object.__setattr__(self, 'hi', _ceil(hi))

    @classmethod
    def point(cls, x):
        return cls(_q(x), _q(x))

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __add__(self, other):
        y = other if isinstance(other, Interval) else Interval.point(other)
        return Interval(self.lo+y.lo, self.hi+y.hi)

    __radd__ = __add__

    def __sub__(self, other):
        return self+-as_interval(other)

    def __rsub__(self, other):
        return as_interval(other)+-self

    def __mul__(self, other):
        y = as_interval(other)
        products = (self.lo*y.lo, self.lo*y.hi, self.hi*y.lo, self.hi*y.hi)
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ValueError('division interval contains zero')
        return Interval(1/self.hi, 1/self.lo)

    def __truediv__(self, other):
        return self*as_interval(other).reciprocal()

    def __rtruediv__(self, other):
        return as_interval(other)*self.reciprocal()

    def __pow__(self, n):
        if type(n) is not int or not 0 <= n <= 512:
            raise ValueError('bounded nonnegative integer power required')
        value = Interval.point(1)
        for _ in range(n):
            value *= self
        return value

    def json(self):
        return [str(self.lo), str(self.hi)]


def as_interval(x):
    return x if isinstance(x, Interval) else Interval.point(x)


def sqrt_integer(n):
    if type(n) is not int or not 0 <= n <= 100:
        raise ValueError('bounded nonnegative integer radicand required')
    root = isqrt(n*SCALE*SCALE)
    upper = root if root*root == n*SCALE*SCALE else root+1
    return Interval(F(root, SCALE), F(upper, SCALE))


def _exp_point(x):
    x = _q(x)
    if abs(x) > 4:
        raise ValueError('exponent outside the bounded certificate domain')
    if x < 0:
        return _exp_point(-x).reciprocal()
    term = total = Interval.point(1)
    for k in range(1, 97):
        term = term*x/k
        total += term
    # Every successive omitted term ratio is at most x/98.
    next_term = term*x/97
    tail = next_term/(1-x/F(98))
    return Interval(total.lo, total.hi+tail.hi)


def exp(x):
    x = as_interval(x)
    return Interval(_exp_point(x.lo).lo, _exp_point(x.hi).hi)


def _log_point(x):
    x = _q(x)
    if not F(1, 16) <= x <= 16:
        raise ValueError('positive logarithm argument in [1/16,16] required')
    if x < 1:
        return -_log_point(1/x)
    z = Interval.point((x-1)/(x+1))
    square, term, total = z*z, z, Interval.point(0)
    for k in range(160):
        total += 2*term/(2*k+1)
        term *= square
    # log(x)=2 atanh(z); all omitted denominators are at least 321.
    tail = 2*term/(321*(1-square))
    return Interval(total.lo, total.hi+tail.hi)


def log(x):
    x = as_interval(x)
    return Interval(_log_point(x.lo).lo, _log_point(x.hi).hi)


def distance_from_zero(x):
    x = as_interval(x)
    return max(F(), x.lo, -x.hi)
