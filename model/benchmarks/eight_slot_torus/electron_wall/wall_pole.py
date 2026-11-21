"""Exact dyadic enclosures of the squared rest pole and the magnetic deficit."""
import argparse
import json
from fractions import Fraction


def continuants(length):
    """Coefficients in x of det(M M^T - x I) and its a derivative at a=1."""
    p0, p = [1], [1, -1]
    q0, q = [0], [2, 0]
    def recurrence(prev, curr):
        coeff = [0] * (len(curr) + 1)
        for i, value in enumerate(curr):
            coeff[i] += 5 * value
            coeff[i + 1] -= value
        for i, value in enumerate(prev):
            coeff[i] -= 4 * value
        return coeff

    for _ in range(1, length):
        next_p = recurrence(p0, p)
        next_q = recurrence(q0, q)
        for i, value in enumerate(p):
            next_q[i] += 2 * value
        for i, value in enumerate(p0):
            next_q[i] -= 8 * value
        p0, p, q0, q = p, next_p, q, next_q
    return p, q


def value_interval(coeff, lower, upper):
    """Outward exact Horner evaluation for 0 <= lower <= upper."""
    lo = hi = Fraction(0)
    for term in reversed(coeff):
        lo, hi = (lo * (upper if lo < 0 else lower) + term,
                  hi * (lower if hi < 0 else upper) + term)
    return lo, hi


def enclose(length, precision=48):
    if type(length) is not int or length < 2 or type(precision) is not int or precision < 8:
        raise ValueError("the isolated wall pole requires length >= 2 and at least eight guard bits")
    bits = 2 * length + precision
    unit = 1 << bits
    lo, hi = 0, unit
    while hi - lo > 1:
        mid = (lo + hi) // 2
        previous, current = 1, unit - mid
        for _ in range(1, length):
            previous, current = current, (5 * unit - mid) * current - 4 * unit ** 2 * previous
        # Interlacing gives a single simple root in (0,1).
        if current > 0:
            lo = mid
        else:
            hi = mid
    lower, upper = Fraction(lo, unit), Fraction(hi, unit)
    p, q = continuants(length)
    px = [i * value for i, value in enumerate(p) if i]
    dx_lo, dx_hi = value_interval(px, lower, upper)
    da_lo, da_hi = value_interval(q, lower, upper)
    if dx_hi >= 0 or da_lo <= 0:
        raise ArithmeticError("root interval needs more precision for an anomaly bound")
    # t = m s'(1) = -P_a / (2 P_x), deficit = 1 - g/2 = t/(1+t).
    t_lo, t_hi = da_lo / (-2 * dx_lo), da_hi / (-2 * dx_hi)
    deficit_lo = unit * t_lo / (1 + t_lo)
    deficit_hi = unit * t_hi / (1 + t_hi)
    return dict(length=length, denominator_power=bits, mass_squared=[lo, hi],
                magnetic_deficit=[deficit_lo.__floor__(), deficit_hi.__ceil__()])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("length", type=int, nargs="*", default=[2, 3, 4, 6, 8, 16, 32, 64])
    args = parser.parse_args()
    print(json.dumps([enclose(length) for length in args.length], indent=2))
