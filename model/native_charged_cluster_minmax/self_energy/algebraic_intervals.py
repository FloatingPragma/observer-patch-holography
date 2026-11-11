"""Rational enclosures of real algebraic coefficients, with exact root rounding.

Rational powers use the real nonnegative branch. Negative bases with
noninteger exponents are rejected rather than selecting a complex branch.
All arithmetic determining the bounds is integer or Fraction arithmetic.
"""
from fractions import Fraction
import sympy as sp


def nth_root_bounds(value, degree, digits=60):
    """Enclose a nonnegative rational's positive degree-th root.

    Let S=10**digits and k=floor((value*S**degree)**(1/degree)).
    Then k/S <= root(value) <= (k+1)/S; the upper endpoint equals the
    lower endpoint when its power is exactly value. The floor integer root
    is unchanged if value*S**degree is first rounded down to an integer.
    """
    value = Fraction(value)
    if not isinstance(degree, int) or degree < 1:
        raise ValueError('root degree must be a positive integer')
    if not isinstance(digits, int) or digits < 1:
        raise ValueError('decimal precision must be a positive integer')
    if value < 0:
        raise ValueError('nonnegative real root required')
    if value == 0:
        return Fraction(0), Fraction(0)
    scale = 10**digits
    scaled_numerator = value.numerator * scale**degree
    scaled_floor = scaled_numerator // value.denominator
    lower_integer, _ = sp.integer_nthroot(scaled_floor, degree)
    lower = Fraction(int(lower_integer), scale)
    exact = int(lower_integer)**degree * value.denominator == scaled_numerator
    upper = lower if exact else Fraction(int(lower_integer) + 1, scale)
    assert lower**degree <= value <= upper**degree
    return lower, upper


def integer_power_bounds(lower, upper, power):
    """Range of x**power on a closed rational interval, excluding poles."""
    if lower > upper:
        raise ValueError('reversed interval')
    if power == 0:
        return Fraction(1), Fraction(1)
    if power < 0:
        if lower <= 0 <= upper:
            raise ValueError('negative power interval contains zero')
        lower, upper = Fraction(1, upper), Fraction(1, lower)
        power = -power
    if power % 2:
        return lower**power, upper**power
    endpoints = (lower**power, upper**power)
    return (Fraction(0) if lower <= 0 <= upper else min(endpoints), max(endpoints))


def interval_expr(expr, digits=60):
    """Enclose a real expression of rationals, sums, products and powers."""
    expr = sp.sympify(expr)
    if expr.is_Rational:
        value = Fraction(int(expr.p), int(expr.q))
        return value, value
    if expr.is_Add:
        parts = [interval_expr(part, digits) for part in expr.args]
        return sum((part[0] for part in parts), Fraction(0)), sum((part[1] for part in parts), Fraction(0))
    if expr.is_Mul:
        lower = upper = Fraction(1)
        for part in expr.args:
            a, b = interval_expr(part, digits)
            products = (lower*a, lower*b, upper*a, upper*b)
            lower, upper = min(products), max(products)
        return lower, upper
    if expr.is_Pow:
        base, exponent = expr.args
        lower, upper = interval_expr(base, digits)
        if not exponent.is_Rational:
            raise ValueError('only rational powers are admitted')
        numerator, denominator = int(exponent.p), int(exponent.q)
        if denominator == 1:
            return integer_power_bounds(lower, upper, numerator)
        if lower < 0:
            raise ValueError('noninteger power requires a nonnegative base interval')
        if numerator < 0 and lower == 0:
            raise ValueError('negative rational power interval contains zero')
        lower = nth_root_bounds(lower, denominator, digits)[0]
        upper = nth_root_bounds(upper, denominator, digits)[1]
        return integer_power_bounds(lower, upper, numerator)
    raise ValueError(f'unsupported real algebraic expression: {expr}')
