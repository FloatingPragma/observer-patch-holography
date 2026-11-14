"""Small outward dyadic interval kernel with integer-only transcendental bounds.

The represented scale is fixed at 2**160. Sin/cos use finite Taylor sums
and a geometric absolute tail; square roots use integer square root. This
kernel neither calls a host transcendental routine nor trusts a residual.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, isqrt

BITS = 160
SCALE = 1 << BITS


def rational(x):
    if type(x) not in (int, F):
        raise ValueError('exact integer or Fraction required')
    return F(x)


def ceiling(a, b):
    return -((-a)//b)


@dataclass(frozen=True)
class Interval:
    lo: int
    hi: int

    def __post_init__(self):
        if type(self.lo) is not int or type(self.hi) is not int or self.lo > self.hi:
            raise ValueError('ordered integer dyadic endpoints required')

    @classmethod
    def point(cls, x):
        x = rational(x)
        return cls(x.numerator*SCALE//x.denominator,
                   ceiling(x.numerator*SCALE, x.denominator))

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __add__(self, other):
        other = convert(other)
        return Interval(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __sub__(self, other):
        return self+-convert(other)

    def __rsub__(self, other):
        return convert(other)+-self

    def __mul__(self, other):
        other = convert(other)
        products = (self.lo*other.lo, self.lo*other.hi,
                    self.hi*other.lo, self.hi*other.hi)
        return Interval(min(products)//SCALE, ceiling(max(products), SCALE))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = rational(other)
        if other <= 0:
            raise ValueError('strictly positive exact divisor required')
        return Interval(self.lo*other.denominator//other.numerator,
                        ceiling(self.hi*other.denominator, other.numerator))

    def absolute_upper(self):
        return F(max(abs(self.lo), abs(self.hi)), SCALE)

    def json(self):
        return [str(F(self.lo, SCALE)), str(F(self.hi, SCALE))]


def convert(x):
    return x if type(x) is Interval else Interval.point(x)


@lru_cache(maxsize=65536)
def sqrt_fraction(x):
    x = rational(x)
    if x < 0:
        raise ValueError('nonnegative exact radicand required')
    n = isqrt((x.numerator*SCALE*SCALE)//x.denominator)
    upper = n if n*n*x.denominator == x.numerator*SCALE*SCALE else n+1
    return Interval(n, upper)


@lru_cache(maxsize=65536)
def sine_cosine(x):
    x = rational(x)
    if abs(x) > 4:
        raise ValueError('the trigonometric certificate admits |angle| <= 4')
    value = Interval.point(x)
    square = value*value
    sine_term, cosine_term = value, Interval.point(1)
    sine, cosine = sine_term, cosine_term
    for n in range(1, 64):
        sine_term = -(sine_term*square)/((2*n)*(2*n+1))
        cosine_term = -(cosine_term*square)/((2*n-1)*(2*n))
        sine += sine_term
        cosine += cosine_term
    # Cosine's next degree is 128, sine's is 129. At |x|<=4 the
    # former dominates; every following same-parity ratio <=16/(129*130).
    tail = F(4**128, factorial(128))/(1-F(16, 129*130))
    error = ceiling(tail.numerator*SCALE, tail.denominator)
    return Interval(sine.lo-error, sine.hi+error), Interval(cosine.lo-error, cosine.hi+error)


def _atan_inverse(n):
    """Alternating atan(1/n), with the first omitted term as a remainder."""
    if type(n) is not int or n < 2:
        raise ValueError('integer inverse tangent denominator at least two required')
    total = sum((F((-1)**k, (2*k+1)*n**(2*k+1)) for k in range(80)), F())
    next_term = F(1, 161*n**161)
    # 80 terms end with a negative term, so the partial sum is a lower bound.
    lower, upper = Interval.point(total), Interval.point(total+next_term)
    return Interval(lower.lo, upper.hi)


@lru_cache(maxsize=1)
def pi():
    """Machin identity; exact alternating-series enclosures."""
    return 16*_atan_inverse(5)-4*_atan_inverse(239)


def pi_angle_error(actual, multiple):
    actual, multiple = rational(actual), rational(multiple)
    return (Interval.point(actual)-multiple*pi()).absolute_upper()
