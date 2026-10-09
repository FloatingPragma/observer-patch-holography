"""Original-input scalar controls for lost vacuum-production signals.

The second-order oracle uses exact rational Chebyshev polynomials.  The
fourth-order oracle executes the elementary signed gates on a covariance at
independently chosen high precision; it imports no implementation helpers.
"""

from decimal import Decimal
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import sys

import mpmath as mp
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from m1_vacuum_fidelity import check, model


AXIS3 = ((F(1), 1), (F(28), 6), (F(55), 12), (F(82), 8))
NEAR_REVIVAL = F(
    24533837163709007127477488471341329157372681512876460105061276634582969781980378411110383555900361877,
    250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000,
)
# A rational decimal near trace(S4(z))/2 = cos(pi/32).  It is supplied as
# written; the reference never substitutes the exact transcendental revival.
FOURTH_ORDER_NEAR_REVIVAL = F(
    "0.098175374644468212088036360922066774532021049431377967701067713855249152709985645386084491613436792731"
)
EVALUATORS = (model.quantum_case, check.quantum)


def scalar(ctx, value):
    return ctx.mpf(value.numerator)/value.denominator


@lru_cache(None)
def original_observables(spectrum, tick, order):
    # Resolve the small trace excess before subtracting the original vacuum.
    # Squared steps are formed from the supplied rationals, never an exported
    # or implementation-generated matrix.  This generous precision is only an
    # independent test oracle, not the implementation's numerical policy.
    squares = [tick*tick*value for value, _ in spectrum]
    smallness = max(len(str(s.denominator))-len(str(s.numerator)) for s in squares)
    ctx = mp.mp.clone()
    ctx.dps = max(180, 6*smallness+100, 2*len(str(tick.denominator))+100)
    rows = [[ctx.mpf(0) for _ in range(3)] for _ in range(32)]
    for (value, multiplicity), square in zip(spectrum, squares):
        omega = ctx.sqrt(scalar(ctx, value))
        if order == 2:
            a = 1-square/2
            previous, current = F(0), F(1)
            populations = []
            for _ in range(32):
                populations.append(scalar(ctx, square**3*current**2/64))
                previous, current = current, 2*a*current-previous
        else:
            z = scalar(ctx, tick)*omega
            root = ctx.root(2, 3)
            weight = 1/(2-root)
            # Original ground-state covariance is I/2.  Store twice that
            # covariance and evolve through K,D,K independently of mode_step.
            covariance = ctx.eye(2)
            populations = []
            for _ in range(32):
                for coefficient in (weight, -root*weight, weight):
                    kick = ctx.matrix([[1, 0], [-coefficient*z/2, 1]])
                    drift = ctx.matrix([[1, coefficient*z], [0, 1]])
                    for gate in (kick, drift, kick):
                        covariance = gate*covariance*gate.T
                populations.append((covariance[0, 0]+covariance[1, 1]-2)/4)
        for row, number in zip(rows, populations):
            row[0] += multiplicity*number
            row[1] += multiplicity*omega*number
            row[2] += multiplicity*ctx.log1p(number)/2
    return ctx, rows


def assert_observables(actual, expected, ctx):
    assert actual["times"] == list(range(1, 33))
    assert len(actual["observations"]) == 32
    desired = expected+[[sum(row[i] for row in expected)/32 for i in range(3)]]
    reported = actual["observations"]+[actual["mean"]]
    for row_index, (reported_row, expected_row) in enumerate(zip(reported, desired), 1):
        assert len(reported_row) == 3
        for column, (observed, correct) in enumerate(zip(reported_row, expected_row)):
            assert isinstance(observed, str)
            value = ctx.mpf(observed)
            assert ctx.isfinite(value)
            assert correct > 0
            error = abs(value/correct-1)
            if error >= ctx.mpf("5e-23"):
                pytest.fail(
                    f"row={row_index}, observable={column}, actual={observed}, "
                    f"expected={ctx.nstr(correct, 40)}, relative_error={ctx.nstr(error, 10)}"
                )
            assert len(Decimal(observed).as_tuple().digits) == 23


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("spectrum,tick,order", [
    pytest.param(AXIS3, F(1, 10**40), 2, id="second-order-small-tick"),
    pytest.param(AXIS3, F(1, 10**20), 4, id="fourth-order-small-tick"),
    pytest.param(((F(10**484), 1),), F(1, 10**282), 2, id="ordinary-energy"),
    pytest.param(((F(1), 10**242),), F(1, 10**40), 2, id="ordinary-total"),
    pytest.param(((F(1), 1),), NEAR_REVIVAL, 2, id="rational-near-revival"),
])
def test_original_input_observables(evaluate, spectrum, tick, order):
    ctx, expected = original_observables(spectrum, tick, order)
    assert_observables(evaluate(spectrum, tick, order), expected, ctx)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("order,exponent", [
    (2, 10), (2, 20), (2, 30), (2, 33), (2, 34), (2, 100),
    (4, 6), (4, 10), (4, 12), (4, 14), (4, 16), (4, 18), (4, 40), (4, 100),
])
def test_small_signals_cross_lost_intermediate_precision(evaluate, order, exponent):
    tick = F(1, 10**exponent)
    ctx, expected = original_observables(AXIS3, tick, order)
    assert_observables(evaluate(AXIS3, tick, order), expected, ctx)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("spectrum,tick,order", [
    (((F(10**404), 1),), F(1, 10**222), 4),
    (((F(1), 10**202),), F(1, 10**20), 4),
    (((F(1, 10**484), 1),), F(10**202), 2),
    (((F(1, 10**404), 1),), F(10**182), 4),
    (((F(1, 4), 3), (F(1), 10**240), (F(82), 11)), F(1, 10**40), 2),
    (((F(1, 4), 3), (F(1), 10**200), (F(82), 11)), F(1, 10**20), 4),
])
def test_rescaling_and_mixed_modes_preserve_each_observable(evaluate, spectrum, tick, order):
    ctx, expected = original_observables(spectrum, tick, order)
    assert_observables(evaluate(spectrum, tick, order), expected, ctx)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("order,center", [(2, NEAR_REVIVAL), (4, FOURTH_ORDER_NEAR_REVIVAL)])
@pytest.mark.parametrize("side", [-1, 0, 1])
def test_near_revivals_resolve_the_supplied_rational(evaluate, order, center, side):
    spectrum = ((F(1), 1),)
    tick = center+side*F(1, 10**90)
    ctx, expected = original_observables(spectrum, tick, order)
    assert 0 < expected[31][0] < ctx.mpf("1e-180")
    assert_observables(evaluate(spectrum, tick, order), expected, ctx)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("order,exponent", [(2, 40), (4, 20)])
@pytest.mark.parametrize("ambient_dps", [7, 90, 300])
def test_scalar_results_isolate_and_restore_arithmetic_context(evaluate, order, exponent, ambient_dps):
    tick = F(1, 10**exponent)
    ctx, expected = original_observables(AXIS3, tick, order)
    with mp.workdps(ambient_dps):
        previous_iv_precision = mp.iv.prec
        mp.iv.prec = 37 if ambient_dps == 7 else 137
        try:
            before = (mp.mp.prec, mp.mp.trap_complex, mp.mp.pretty, mp.iv.prec, mp.iv.pretty)
            actual = evaluate(AXIS3, tick, order)
            assert (mp.mp.prec, mp.mp.trap_complex, mp.mp.pretty, mp.iv.prec, mp.iv.pretty) == before
        finally:
            mp.iv.prec = previous_iv_precision
    assert_observables(actual, expected, ctx)


def test_committed_catalog_agrees_with_original_input_control():
    packet = json.loads((Path(__file__).with_name("receipt.json")).read_text(encoding="utf-8"))
    for graph in packet["evidence"]["graphs"].values():
        spectrum = tuple((F(value), multiplicity) for value, multiplicity in graph["spectrum"])
        for key, case in graph["quantum"].items():
            order, tick = key.split(":")
            ctx, expected = original_observables(spectrum, F(tick), int(order))
            assert_observables(case, expected, ctx)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("order", [2, 4])
@pytest.mark.parametrize("side", [-1, 0, 1])
def test_original_exact_tick_domain(evaluate, order, side):
    tick = F(1, 10)+side*F(1, 10**100)
    if side > 0:
        with pytest.raises(ValueError, match="domain"):
            evaluate(((F(1), 1),), tick, order)
    else:
        ctx, expected = original_observables(((F(1), 1),), tick, order)
        assert_observables(evaluate(((F(1), 1),), tick, order), expected, ctx)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("order", [2, 4])
@pytest.mark.parametrize("side", [-1, 0, 1])
def test_original_exact_eigenvalue_domain(evaluate, order, side):
    spectrum = ((F(1)+side*F(1, 10**100), 1),)
    if side > 0:
        with pytest.raises(ValueError, match="domain"):
            evaluate(spectrum, F(1, 10), order)
    else:
        ctx, expected = original_observables(spectrum, F(1, 10), order)
        assert_observables(evaluate(spectrum, F(1, 10), order), expected, ctx)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("spectrum,tick,order", [
    (((1, 3), (F(1), 7), (F(4), 2)), F(1, 256), 2),
    (((F(4), 2), (F(1), 7), (1, 3)), F(1, 256), 4),
    (((F(1, 100), 1),), 1, 2),
    (((F(1, 100), 1),), 1, 4),
])
def test_integer_scalars_and_repeated_modes_are_valid(evaluate, spectrum, tick, order):
    ctx, expected = original_observables(spectrum, F(tick), order)
    assert_observables(evaluate(spectrum, tick, order), expected, ctx)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("bad", [True, False, 0.01, Decimal("0.01"), "1/100", None, 1j])
@pytest.mark.parametrize("field", ["tick", "eigenvalue"])
def test_inexact_or_non_numeric_scalars_are_refused(evaluate, field, bad):
    spectrum, tick = ((F(1), 1),), F(1, 256)
    if field == "tick":
        tick = bad
    else:
        spectrum = ((bad, 1),)
    with pytest.raises((TypeError, ValueError)):
        evaluate(spectrum, tick, 2)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("bad", [0, -1, F(0), F(-1, 10**100)])
@pytest.mark.parametrize("field", ["tick", "eigenvalue"])
def test_nonpositive_scalars_are_refused(evaluate, field, bad):
    spectrum, tick = ((F(1), 1),), F(1, 256)
    if field == "tick":
        tick = bad
    else:
        spectrum = ((bad, 1),)
    with pytest.raises(ValueError):
        evaluate(spectrum, tick, 2)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("multiplicity", [0, -1, True, False, F(1), 1.0, Decimal(1), "1", None])
def test_multiplicity_requires_positive_python_integer(evaluate, multiplicity):
    with pytest.raises((TypeError, ValueError)):
        evaluate(((F(1), multiplicity),), F(1, 256), 2)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("order", [0, 1, 3, -2, True, F(2), 2.0, "2", None])
def test_order_requires_supported_python_integer(evaluate, order):
    with pytest.raises((TypeError, ValueError)):
        evaluate(((F(1), 1),), F(1, 256), order)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("spectrum", [(), [], ((F(1),),), ((F(1), 1, 2),), (None,)])
def test_empty_or_malformed_spectrum_is_refused(evaluate, spectrum):
    with pytest.raises((TypeError, ValueError)):
        evaluate(spectrum, F(1, 256), 2)


@pytest.mark.parametrize("evaluate", EVALUATORS)
def test_domain_refusal_preserves_arithmetic_context(evaluate):
    with mp.workdps(11):
        previous_iv_precision = mp.iv.prec
        mp.iv.prec = 37
        try:
            before = (mp.mp.prec, mp.mp.trap_complex, mp.mp.pretty, mp.iv.prec, mp.iv.pretty)
            with pytest.raises(ValueError):
                evaluate(((F(1), 1),), F(1, 10)+F(1, 10**100), 2)
            assert (mp.mp.prec, mp.mp.trap_complex, mp.mp.pretty, mp.iv.prec, mp.iv.pretty) == before
        finally:
            mp.iv.prec = previous_iv_precision


def exact_interval_bounds(value):
    """Read binary endpoints as rationals without an ambient mp conversion."""
    result = []
    for sign, mantissa, exponent, bit_count in value._mpi_:
        assert bit_count >= 0  # The oracle tests only finite intervals.
        value = F(-mantissa if sign else mantissa)
        result.append(value*2**exponent if exponent >= 0 else value/F(2**(-exponent)))
    return result


@pytest.mark.parametrize("value", [F(1, 3), F(-1, 3), F(1, 10**1000), F(10**1000)])
def test_interval_rational_conversion_encloses_original_fraction(value):
    from mpmath.ctx_iv import MPIntervalContext
    from m1_vacuum_fidelity import numerics

    ctx = MPIntervalContext()
    ctx.dps = 40
    lower, upper = exact_interval_bounds(numerics.rational(ctx, value))
    assert lower <= value <= upper
    assert upper-lower < abs(value)*F(1, 10**35)


@pytest.mark.parametrize("value", [
    F(1, 10), F(1, 2**1074), F(1, 10**400), F(1, 10**1000),
    F(1, 8)-F(1, 10**100), F(1, 8), F(1, 8)+F(1, 10**100), F(10**1000),
])
def test_interval_logarithm_encloses_original_positive_input(value):
    from mpmath.ctx_iv import MPIntervalContext
    from m1_vacuum_fidelity import numerics

    ctx = MPIntervalContext()
    ctx.dps = 50
    lower, upper = exact_interval_bounds(numerics.log1p_positive(ctx, numerics.rational(ctx, value)))
    # Independent log1p receives the original exact fraction.  Computing
    # log(1+x) at this interval precision would erase the three small inputs.
    oracle = mp.mp.clone()
    oracle.dps = 250
    expected = oracle.log1p(scalar(oracle, value))
    assert 0 < scalar(oracle, lower) <= expected <= scalar(oracle, upper)
    assert scalar(oracle, upper-lower) < expected*oracle.mpf("1e-40")


def test_interval_logarithm_keeps_exact_zero_and_refuses_invalid_inputs():
    from mpmath.ctx_iv import MPIntervalContext
    from m1_vacuum_fidelity import numerics

    ctx = MPIntervalContext()
    ctx.dps = 50
    assert exact_interval_bounds(numerics.log1p_positive(ctx, ctx.mpf(0))) == [F(0), F(0)]
    for value in (ctx.mpf(-1), ctx.mpf([-1, 1]), ctx.mpf("inf"), ctx.mpf("nan")):
        with pytest.raises(ValueError, match="nonnegative"):
            numerics.log1p_positive(ctx, value)


def test_interval_logarithm_encloses_both_ends_of_a_positive_interval():
    from mpmath.ctx_iv import MPIntervalContext
    from m1_vacuum_fidelity import numerics

    ctx = MPIntervalContext()
    ctx.dps = 50
    original_lower, original_upper = F(1, 10**400), F(2, 10**400)
    value = ctx.mpf([numerics.rational(ctx, original_lower).a,
                    numerics.rational(ctx, original_upper).b])
    lower, upper = exact_interval_bounds(numerics.log1p_positive(ctx, value))
    oracle = mp.mp.clone()
    oracle.dps = 250
    assert 0 < scalar(oracle, lower) <= oracle.log1p(scalar(oracle, original_lower))
    assert oracle.log1p(scalar(oracle, original_upper)) <= scalar(oracle, upper)


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("order,near_revival", [(2, NEAR_REVIVAL), (4, FOURTH_ORDER_NEAR_REVIVAL)])
def test_precision_budget_refuses_unresolved_revival_and_preserves_context(
        evaluate, order, near_revival, monkeypatch):
    from m1_vacuum_fidelity import numerics

    monkeypatch.setattr(numerics, "PRECISION_STEPS", (90,))
    spectrum = ((F(1), 1),)
    # A budget refusal must leave an ordinary valid control available.
    ctx, expected = original_observables(spectrum, F(1, 256), order)
    assert_observables(evaluate(spectrum, F(1, 256), order), expected, ctx)
    with mp.workdps(11):
        previous_iv_precision = mp.iv.prec
        mp.iv.prec = 37
        try:
            before = (mp.mp.prec, mp.mp.trap_complex, mp.mp.pretty, mp.iv.prec, mp.iv.pretty)
            with pytest.raises(ValueError, match="precision|budget"):
                evaluate(spectrum, near_revival, order)
            assert (mp.mp.prec, mp.mp.trap_complex, mp.mp.pretty, mp.iv.prec, mp.iv.pretty) == before
        finally:
            mp.iv.prec = previous_iv_precision


def test_interval_formatter_requires_all_reported_digits_to_be_resolved():
    from mpmath.ctx_iv import MPIntervalContext
    from mpmath.ctx_mp import MPContext
    from m1_vacuum_fidelity import numerics

    ctx, formatter = MPIntervalContext(), MPContext()
    ctx.dps = 90
    # Positive endpoints on opposite sides of the final reported digit.
    crossing = ctx.mpf(["1.234567890123456789012349", "1.234567890123456789012351"])
    assert numerics._formatted_interval(crossing, formatter) is None
    resolved = ctx.mpf(["1.234567890123456789012341", "1.234567890123456789012342"])
    assert numerics._formatted_interval(resolved, formatter) == "1.2345678901234567890123"


@pytest.mark.parametrize("evaluate", EVALUATORS)
@pytest.mark.parametrize("order,center", [(2, NEAR_REVIVAL), (4, FOURTH_ORDER_NEAR_REVIVAL)])
def test_precision_budget_refuses_positive_but_unresolved_observations(
        evaluate, order, center, monkeypatch):
    from mpmath.ctx_iv import MPIntervalContext
    from mpmath.ctx_mp import MPContext
    from m1_vacuum_fidelity import numerics

    spectrum, tick = ((F(1), 1),), center+F(1, 10**70)
    ctx, expected = original_observables(spectrum, tick, order)
    assert_observables(evaluate(spectrum, tick, order), expected, ctx)
    # Unlike the closer revival control, this budget resolves positivity.
    # The formatter, rather than the zero-containing-interval guard, must
    # refuse it: the last population has fewer than 23 resolved digits.
    interval, formatter = MPIntervalContext(), MPContext()
    interval.dps = 90
    kernel = model._mode_numbers if evaluate is model.quantum_case else check._matrix_mode_numbers
    numbers = list(kernel(interval, spectrum, tick, order))[0][2]
    assert all(number.a > 0 for number in numbers)
    endpoints = [formatter.nstr(formatter.make_mpf(endpoint), 23, strip_zeros=False)
                 for endpoint in numbers[-1]._mpi_]
    assert endpoints[0] != endpoints[1]
    monkeypatch.setattr(numerics, "PRECISION_STEPS", (90,))
    with pytest.raises(ValueError, match="precision|budget"):
        evaluate(spectrum, tick, order)
