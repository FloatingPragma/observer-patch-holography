"""Tests for the charged Z3 phase non-selection theorem."""

from __future__ import annotations

import math

import derive_charged_z3_phase_holonomy_no_go as lane


def test_balance_invariants_do_not_select_phase():
    for delta in (0.0, 2.0 / 9.0, 0.4):
        trace_one, trace_two = lane.trace_invariants(delta)
        assert math.isclose(trace_one, 3.0, rel_tol=0.0, abs_tol=1.0e-14)
        assert math.isclose(trace_two, 6.0, rel_tol=0.0, abs_tol=1.0e-14)


def test_two_distinct_positive_phase_countermodels_exist():
    zero = lane.roots(0.0)
    stage5 = lane.roots(2.0 / 9.0)
    assert all(value > 0.0 for value in zero + stage5)
    assert zero != stage5


def test_phase_remains_unpromoted():
    artifact = lane.build_artifact()
    assert artifact["public_phase_promotion_allowed"] is False
    assert artifact["status"] == "CLOSED_NO_GO_CURRENT_OPH_DATA_LEAVE_PHASE_CONTINUOUS"
