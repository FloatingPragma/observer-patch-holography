"""Exact enclosure regression tests for source field coefficients."""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import pytest
import sympy as s
spec = importlib.util.spec_from_file_location('source_algebraic_intervals', Path(__file__).with_name('algebraic_intervals.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
interval_expr, nth_root_bounds = module.interval_expr, module.nth_root_bounds


def test_previously_unhandled_three_quarters_power():
    lower, upper = interval_expr(2**s.Rational(3,4))
    assert lower**4 <= 8 <= upper**4
    assert 0 < upper-lower < F(1,10**58)
    inverse_lower, inverse_upper = interval_expr(2**s.Rational(-3,4))
    assert inverse_lower**4 <= F(1,8) <= inverse_upper**4


@pytest.mark.parametrize('degree', range(1,10))
def test_root_bounds_straddle_exact_rationals(degree):
    for value in (F(0), F(1,19), F(2), F(1001,17), F(7**degree,11**degree)):
        lower, upper = nth_root_bounds(value, degree, digits=24)
        assert lower**degree <= value <= upper**degree
        assert upper-lower <= F(1,10**24)


@pytest.mark.parametrize('power',[s.Rational(1,2),s.Rational(-1,2),s.Rational(1,4),s.Rational(-1,4),s.Rational(5,6),s.Rational(-7,3)])
def test_power_endpoint_contract(power):
    value=s.Rational(37,23);lower,upper=interval_expr(value**power)
    exact=F(37,23)**int(power.p)
    assert lower**int(power.q) <= exact <= upper**int(power.q)


def test_perfect_roots_and_signed_integer_powers_are_exact():
    assert nth_root_bounds(F(81,16),4)==(F(3,2),F(3,2))
    assert module.integer_power_bounds(F(-3),F(2),2)==(F(0),F(9))
    assert module.integer_power_bounds(F(-3),F(-2),-3)==(F(-1,8),F(-1,27))
    assert interval_expr(s.Rational(9,4)**s.Rational(3,2))==(F(27,8),F(27,8))


def test_actual_nested_field_variance_encloses_exact_identity():
    lower,upper=interval_expr(s.Rational(1,16)+1/(12*s.sqrt(6)))
    # Transform back to sqrt(6) using a decreasing positive map.
    assert lower>F(1,16)
    root_lower=1/(12*(upper-F(1,16)));root_upper=1/(12*(lower-F(1,16)))
    assert root_lower**2 <= 6 <= root_upper**2


def test_invalid_branches_and_poles_rejected():
    with pytest.raises(ValueError):nth_root_bounds(F(-1),3)
    with pytest.raises(ValueError):module.integer_power_bounds(F(-1),F(1),-1)
    with pytest.raises(ValueError):interval_expr(s.Pow(-2,s.Rational(1,3),evaluate=False))
    with pytest.raises(ValueError):interval_expr(s.Pow(2,s.pi,evaluate=False))
