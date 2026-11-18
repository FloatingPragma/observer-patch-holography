"""Rational enclosures, with explicit positive-series remainder bounds."""
from fractions import Fraction as F
from math import isqrt

GRID = 10**30


def outward(lo, hi):
    return F((lo*GRID).__floor__(), GRID), F((hi*GRID).__ceil__(), GRID)


def add(a, b):
    return outward(a[0]+b[0], a[1]+b[1])


def sub(a, b):
    return outward(a[0]-b[1], a[1]-b[0])


def mul(a, b):
    products = [x*y for x in a for y in b]
    return outward(min(products), max(products))


def div(a, b):
    if not 0 < b[0] <= b[1]:
        raise ValueError('strictly positive denominator required')
    return mul(a, (1/b[1], 1/b[0]))


def sqrt_point(x):
    if x < 0:
        raise ValueError('nonnegative square root required')
    q = isqrt(x.numerator*GRID**2//x.denominator)
    lo = F(q, GRID)
    return lo, lo if lo*lo == x else F(q+1, GRID)


def minus_log_point(x, branch=0):
    if branch != 0 or not F(1, 2) <= x <= 1:
        raise ValueError('positive real principal logarithm on [1/2,1] required')
    z = (1-x)/(1+x)
    # -log x = 2 atanh z. Every term is positive; the remaining
    # denominators are >=2N+1 and their powers form a geometric series.
    n = 40
    lo = 2*sum((z**(2*j+1)/(2*j+1) for j in range(n)), F())
    tail = 2*z**(2*n+1)/((2*n+1)*(1-z*z))
    return outward(lo, lo+tail)


def acos_point(x):
    if not F(1, 2) <= x <= 1:
        raise ValueError('principal phase on [0,pi/3] required')
    lo, hi = sqrt_point((1-x)/2)
    # acos x = 2 asin sqrt((1-x)/2). The positive asin coefficients
    # are <=1. After N terms its tail is <=y^(2N+1)/(1-y²).
    def lower(y):
        term, total = y, y
        for j in range(20-1):
            term *= F((2*j+1)**2, 2*(j+1)*(2*j+3))*y*y
            total += term
        return 2*total
    return outward(lower(lo), lower(hi)+2*hi**41/(1-hi*hi))


def spectral(transfer):
    lo, hi = transfer
    if not F(1, 2) <= lo <= hi <= 1:
        raise ValueError('positive admitted transfer spectrum required')
    e = outward(1260*minus_log_point(hi)[0], 1260*minus_log_point(lo)[1])
    theta = outward(acos_point(hi)[0], acos_point(lo)[1])
    return e, theta


def disjoint(a, b):
    return a[1] < b[0] or b[1] < a[0]


def trig_point(x, sine=False):
    if not 0 <= x <= 1:
        raise ValueError('Taylor trigonometric domain is [0,1]')
    from math import factorial
    shift = int(sine)
    upper = sum(((-1)**j*x**(2*j+shift)/factorial(2*j+shift)
                 for j in range(15)), F())
    lower = upper-x**(30+shift)/factorial(30+shift)
    return outward(lower, upper)


def sin(a):
    return outward(trig_point(a[0], True)[0], trig_point(a[1], True)[1])


def cos(a):
    return outward(trig_point(a[1])[0], trig_point(a[0])[1])


def sqrt(a):
    return outward(sqrt_point(a[0])[0], sqrt_point(a[1])[1])


def strings(a):
    return list(map(str, a))
