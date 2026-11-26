"""Two-tangent test for a common baseline field with variable scalar wavenumbers."""

from mpmath import iv

from .readout import quadratures


def cross(a, b):
    return [
        a[(j + 1) % 3] * b[(j + 2) % 3] - a[(j + 2) % 3] * b[(j + 1) % 3]
        for j in range(3)
    ]


def outward(value):
    import math

    return [
        math.nextafter(float(value.a), -math.inf),
        math.nextafter(float(value.b), math.inf),
    ]


def tangent(formed, destination):
    v = [(y / x).imag for x, y in zip(formed, destination, strict=True)]
    return [s - v[0] for s in v[1:]]


def rule_out(stocks):
    def directions(q):
        f, p, o = (stocks[q][s] for s in ("formed", "product", "ordinary"))
        return tangent(f, p), tangent(f, o)

    g1, o1 = directions(-1)
    g2, o2 = directions(1)
    constant = cross(g2, o2)
    a, b = cross(g1, o2), cross(g2, o1)
    linear = [-u - v for u, v in zip(a, b, strict=True)]
    square = cross(g1, o1)
    det = sum(x * y for x, y in zip(constant, cross(linear, square), strict=True))
    return dict(
        determinant=outward(det),
        columns=[[outward(v) for v in col] for col in (constant, linear, square)],
    )


def counts_certificate(receipt, systematic):
    stocks = {}
    for charge in (-1, 1):
        stocks[charge] = {}
        for name in ("formed", "product", "ordinary"):
            z, r = quadratures(
                receipt["counts"][str(charge)][name], receipt["shots"], systematic
            )
            # Each empirical complex disk is contained in this rectangle.
            stocks[charge][name] = [
                iv.mpc(
                    iv.mpf(u.real) + iv.mpf([-w - 1e-16, w + 1e-16]),
                    iv.mpf(u.imag) + iv.mpf([-w - 1e-16, w + 1e-16]),
                )
                for u, w in zip(z, r, strict=True)
            ]
    return rule_out(stocks)
