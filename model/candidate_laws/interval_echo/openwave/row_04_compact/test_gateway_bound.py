from fractions import Fraction
import json

import numpy as np
import pytest
from scipy.linalg import expm

from model.executable_model import cell
from .gateway_bound import (
    ATTEMPTS, ROUNDS, SCALE, SCORE, board_charts, certify, dark_projector, internal_pulse, success_intervals,
)
from .gateway_release import HERE, experiment, verify


def test_full_gateway_receipt_against_native_acquisition():
    stored = json.loads((HERE/"gateway_results.json").read_text(encoding="utf-8"))
    actual = experiment()
    verify(actual, stored)
    stored["layouts"][0]["pulses"]["1.0"]["joint_home_floor_ppt"] += 1
    with pytest.raises(ValueError):
        verify(actual, stored)


def test_every_native_chart_and_all_three_cadences():
    boards = board_charts()
    assert {n: sum(rows.values()) for n, rows in boards.items()} == {20: 20, 80: 80, 320: 320}
    charts = set().union(*boards.values())
    assert len(charts) == 10
    for duration in (Fraction(1, 2), Fraction(1), Fraction(2)):
        for tips in charts:
            result = certify(tips, duration)
            assert result["joint_home_floor_ppt"] > 94*SCORE//100
            assert 0 < result["exhaustion_ceiling_ppt"] < SCORE//1000


def test_one_step_enclosure_and_conditional_weights_from_a_second_solver():
    tips = (1, 3, 10)
    adjacency = cell().adjacency().toarray()
    for duration in (Fraction(1, 2), Fraction(1), Fraction(2)):
        re, im, units = internal_pulse(duration)
        expected = expm(.5j*float(duration)*adjacency)
        computed = (re.astype(float) + 1j*im.astype(float))/SCALE
        # The independent dense double-precision exponential has its own rounding.
        assert np.linalg.norm(computed-expected, 2) < 3e-14
        assert 0 < units/SCALE < 1e-20
        expected[list(tips)] = 0
        columns = np.linalg.matrix_power(expected, ROUNDS)
        success = (abs(columns)**2).sum(axis=0)
        intervals = success_intervals(tips, duration)
        for value, (lo, hi) in zip(success, intervals, strict=True):
            assert float(lo)-4e-12 <= value <= float(hi)+4e-12
        diagonal = np.array(dark_projector(tips).diagonal(), dtype=float).ravel()
        floor = 0.
        reach = 1.
        for attempt in range(ATTEMPTS):
            port = attempt % 12
            floor += reach*diagonal[port]
            reach *= 1-success[port]
        enclosure = certify(tips, duration)
        assert -4e-12 < floor-enclosure["joint_home_floor_ppt"]/SCORE < 1e-10
        assert reach <= enclosure["exhaustion_ceiling_ppt"]/SCORE+4e-12


@pytest.mark.parametrize("chart", [(1, 1, 3), (3, 1, 10), (1, 3, 12), (True, 3, 10), ()])
def test_reject_malformed_gateway_instruments(chart):
    with pytest.raises(ValueError):
        success_intervals(chart, Fraction(1))
