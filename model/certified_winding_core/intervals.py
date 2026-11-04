"""Rational source intervals; the numerical eigensolver supplies no bounds."""
from fractions import Fraction as F
from math import isqrt

from ..winding_sector_search.algebra import winding


def sqrt_interval(value, bits=64):
    value = F(value)
    if value < 0 or type(bits) is not int or bits < 1:
        raise ValueError('nonnegative rational and positive precision required')
    scale_value = 1 << bits
    lower = isqrt((value.numerator*scale_value*scale_value)//value.denominator)
    lo = F(lower, scale_value)
    hi = lo if lo*lo == value else F(lower+1, scale_value)
    if not lo*lo <= value <= hi*hi:
        raise ArithmeticError('rational square-root enclosure failed')
    return lo, hi


def scale(a, x):
    values = a[0]*x, a[1]*x
    return min(values), max(values)


def divide(a, b):
    if b[0] <= 0:
        raise ArithmeticError('positive interval denominator required')
    values = [x/y for x in a for y in b]
    return min(values), max(values)


def direction(point, degree):
    """Enclose the literal rational-map texture over Q(sqrt(|point|^2))."""
    degree = winding(degree)
    if type(point) is not tuple or len(point) != 3 or any(type(x) is not int for x in point):
        raise ValueError('three integer coordinates required')
    sign, k = (1 if degree > 0 else -1), abs(degree)
    x, y, z = point
    if x == y == 0:
        return ((F(0), F(0)), (F(0), F(0)),
                (F(sign*((z > 0)-(z < 0))),)*2)
    real, imag = 1, 0
    for _ in range(k):
        real, imag = real*x-imag*y, real*y+imag*x
    norm = real*real+imag*imag
    root = sqrt_interval(F(x*x+y*y+z*z))
    b = root[0]-z, root[1]-z
    if b[0] <= 0:
        raise ArithmeticError('the nonpolar homogeneous coordinate was unresolved')
    power = b[0]**k, b[1]**k
    square = power[0]**2, power[1]**2
    denominator = norm+square[0], norm+square[1]
    nx = divide(scale(power, 2*real), denominator)
    ny = divide(scale(power, 2*imag), denominator)
    nz = divide((norm-square[1], norm-square[0]), denominator)
    return tuple(scale(v, sign) for v in (nx, ny, nz))


def dyadic_center(enclosure, bits):
    scale_value = 1 << bits
    integer = round((enclosure[0]+enclosure[1])*scale_value/2)
    center = F(integer, scale_value)
    error = max(abs(center-enclosure[0]), abs(center-enclosure[1]))
    # A common outward denominator prevents an unnecessary least-common-
    # multiple explosion when summing errors from different quadratic fields.
    error = ceil_dyadic(error, 48)
    return integer, error


def ceil_dyadic(value, bits):
    value = F(value)*(1 << bits)
    return F(-(-value.numerator//value.denominator), 1 << bits)
