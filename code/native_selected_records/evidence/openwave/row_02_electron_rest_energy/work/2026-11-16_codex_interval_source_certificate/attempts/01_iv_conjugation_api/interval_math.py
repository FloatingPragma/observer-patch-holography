"""Outward arithmetic for finite source certificates; no source imports."""
import math
from mpmath import iv, mp

iv.dps = 60
mp.dps = 90


def require(condition, message):
    if not condition:
        raise ValueError(message)


def real(x=0):
    return iv.mpf(x)


def complex_(x=0, y=0):
    return iv.mpc(x, y)


def endpoints(x):
    return tuple(mp.mpf(a) for a in x._mpi_)


def upper(x):
    return math.nextafter(float(endpoints(x)[1]), math.inf)


def lower(x):
    return math.nextafter(float(endpoints(x)[0]), -math.inf)


def packed(x):
    lo, hi = x._mpi_
    return {'lo_dyadic': [lo[0], str(lo[1]), lo[2]],
            'hi_dyadic': [hi[0], str(hi[1]), hi[2]],
            'display': [mp.nstr(mp.mpf(lo), 24), mp.nstr(mp.mpf(hi), 24)]}


def width(x):
    return real(x.b) - real(x.a)


def sym_error(x, error):
    e = error.b
    return x + iv.mpf([-e, e])


def atan_small(x, count=110):
    require(abs(x) < real(1)/4, 'atan series domain')
    p = x
    total = real(0)
    for n in range(count):
        total += (-1 if n % 2 else 1)*p/(2*n+1)
        p *= x*x
    return sym_error(total, abs(p)/(2*count+1)/(1-x*x))


def pi_bound():
    return 16*atan_small(real(1)/5) - 4*atan_small(real(1)/239)


def root_of_unity(n):
    """exp(-2*pi*i/n), from series with explicit geometric tail bounds."""
    require(n >= 64, 'root series domain')
    x = 2*pi_bound()/n
    co = real(1)
    si = x
    tc, ts = real(1), x
    for j in range(1, 24):
        tc *= -x*x/((2*j-1)*(2*j))
        ts *= -x*x/((2*j)*(2*j+1))
        co += tc
        si += ts
    ec = abs(tc*x*x/(47*48))/(1-x*x)
    es = abs(ts*x*x/(48*49))/(1-x*x)
    return complex_(sym_error(co, ec), -sym_error(si, es))


def zeros(m, n=None):
    return [[complex_() for _ in range(m if n is None else n)] for _ in range(m)]


def eye(n):
    a = zeros(n)
    for i in range(n):
        a[i][i] = complex_(1)
    return a


def adj(a):
    return [[a[j][i].conjugate() for j in range(len(a))] for i in range(len(a[0]))]


def mm(a, b):
    require(len(a[0]) == len(b), 'matrix shape')
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), complex_())
             for j in range(len(b[0]))] for i in range(len(a))]


def add(a, b):
    return [[x+y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a, x):
    return [[z*x for z in row] for row in a]


def power(a, n):
    out = eye(len(a))
    while n:
        if n & 1:
            out = mm(out, a)
        n //= 2
        if n:
            a = mm(a, a)
    return out


def trace(a):
    return sum((a[i][i] for i in range(len(a))), complex_())


def sandwich(a, rho, b=None):
    return mm(mm(a, rho), adj(a if b is None else b))


def hermitian_trace(a):
    t = trace(a)
    require(0 in t.imag, 'Hermitian trace does not enclose real axis')
    return t.real


def difference_upper(x, saved):
    return upper(abs(x-real(saved)))


def diagonal(a):
    out = zeros(len(a))
    for i, x in enumerate(a):
        out[i][i] = x
    return out
