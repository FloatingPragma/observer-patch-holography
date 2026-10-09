"""Original-input scalar controls for lost vacuum-production signals.

The second-order oracle uses exact rational Chebyshev polynomials.  The
fourth-order oracle executes the elementary signed gates on a covariance at
independently chosen high precision; it imports no implementation helpers.
"""

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
            assert error < ctx.mpf("5e-23"), (
                f"row={row_index}, observable={column}, actual={observed}, "
                f"expected={ctx.nstr(correct, 40)}, relative_error={ctx.nstr(error, 10)}"
            )


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


def test_committed_catalog_agrees_with_original_input_control():
    packet = json.loads((Path(__file__).with_name("receipt.json")).read_text(encoding="utf-8"))
    for graph in packet["evidence"]["graphs"].values():
        spectrum = tuple((F(value), multiplicity) for value, multiplicity in graph["spectrum"])
        for key, case in graph["quantum"].items():
            order, tick = key.split(":")
            ctx, expected = original_observables(spectrum, F(tick), int(order))
            assert_observables(case, expected, ctx)


@pytest.mark.parametrize("side", [-1, 0, 1])
def test_checker_compares_original_exact_domain(side):
    tick = F(1, 10)+side*F(1, 10**100)
    if side > 0:
        with pytest.raises(ValueError, match="domain"):
            check.quantum(((F(1), 1),), tick, 2)
    else:
        ctx, expected = original_observables(((F(1), 1),), tick, 2)
        assert_observables(check.quantum(((F(1), 1),), tick, 2), expected, ctx)
