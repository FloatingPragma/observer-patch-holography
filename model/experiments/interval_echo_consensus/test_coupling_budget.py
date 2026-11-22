import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from . import alias_window, coupling_budget, persistence_forecast


def test_fractional_degree_and_constant_odds():
    t = sp.Symbol("t", positive=True)
    # Unnormalized masses from three independent words of a binary source.
    masses = (t**3, t**2*(1-t), (1-t)**3)
    assert coupling_budget.mixture_exponent(masses, t) == sp.Rational(1, 2)
    assert coupling_budget.mixture_exponent((2*masses[0], *masses[1:]), t) is None
    assert coupling_budget.mixture_exponent((t, t, t), t) is None
    assert coupling_budget.finite_resolution(9)[0] == [8]
    assert coupling_budget.finite_resolution(17)[0] == [8, 16]


def test_native_standard_and_other_material():
    result = coupling_budget.run()
    assert result["selected"] == 8
    assert min(result["alternate_material_second_pulse_tv"]) > .1
    assert result["word_sources"]["8_fair"]["exponent"] is None
    assert result["word_sources"]["16_common"]["degree"] == 17


def test_tail_and_held_counts():
    tail = alias_window.uniform_tail(*alias_window.setup())
    assert list(map(Fraction, tail["squared_norm_bounds"])) == [
        Fraction(17,4), Fraction(107,4), Fraction(5,16), Fraction(7003,16), Fraction(118043,256)]
    assert tail["lower_bound"] > .94
    assert tail["launch_lower"] > .945
    counts = json.loads(Path(__file__).with_name("persistence_counts.json").read_text())
    replay = persistence_forecast.run()
    assert set(counts) == set(replay)
    for arm, rows in counts.items():
        assert replay[arm]["counts"] == rows
        decision = persistence_forecast.infer_floor(rows[:2], rows[2:4], rows[4])
        assert decision["resolved"] == (arm in ("installed", "reversed", "other_alias"))
        assert replay[arm]["forecast"]["resolved"] == decision["resolved"]
        if decision["resolved"]:
            assert decision["home_probability_lower"] > .94036


def test_inference_rejects_broken_acquisitions():
    rows = json.loads(Path(__file__).with_name("persistence_counts.json").read_text())["installed"]
    quarter, irrational, full = rows[:2], rows[2:4], rows[4]
    infer = persistence_forecast.infer_floor
    assert not infer(quarter, irrational, [2**24-1, 1, 0])["resolved"]
    assert not infer(quarter, irrational, [2**24-1, 0, 1])["resolved"]
    assert not infer([[2**24, 0, 0]]*2, irrational, full)["resolved"]
    assert not infer(quarter, quarter, full)["resolved"]
    with pytest.raises(ValueError):
        infer(quarter, irrational, [2**24-1, 0, 0])
    with pytest.raises(ValueError):
        infer(quarter, irrational, np.array(full, float))
