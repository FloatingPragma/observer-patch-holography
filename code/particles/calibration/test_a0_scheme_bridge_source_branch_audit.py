"""Focused tests for the a0 scheme-bridge source-branch audit."""

from __future__ import annotations

import math

import pytest

import derive_a0_scheme_bridge_source_branch_audit as lane


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    out = tmp_path_factory.mktemp("bridge") / "a0_scheme_bridge_source_branch_audit.json"
    return lane.build(out)


def test_one_loop_reproduction(result):
    rep = result["one_loop_reproduction"]
    assert abs(rep["residual"]) < 1e-6
    assert abs(rep["independent_alpha2_check_residual"]) < 1e-3


def test_two_loop_shift_has_right_sign_and_overshoots(result):
    lo, hi = result["anchor_shift_interval_over_yukawa_conventions"]
    gap_lo, gap_hi = result["membership_test"]["certified_gap_interval"]
    assert lo > 0.0, "two-loop shift must move the anchor toward the physical value"
    assert lo > gap_hi, "every Yukawa convention overshoots the certified gap"
    assert result["membership_test"]["any_convention_lands_inside"] is False


def test_yukawa_scan_monotone_decreasing(result):
    shifts = [row["anchor_shift_inv_alpha"] for row in result["two_loop_scan"]]
    assert shifts == sorted(shifts, reverse=True)


def test_balance_requirement_is_negative_and_threshold_sub_10_tev(result):
    balance = result["balance_requirement"]
    assert balance["required_threshold_shift_inv_alpha"] < 0.0
    assert 91.1876 < balance["equivalent_single_effective_threshold_gev"] < 1e4


def test_verdict_fail_closed(result):
    verdict = result["verdict"]
    assert verdict["status"] == "source_branch_not_closable_under_declared_conventions"
    assert verdict["refined_missing_objects"] == ["scheme_lock", "threshold_map"]
    guards = result["guards"]
    assert guards["source_only_bridge_emitted"] is False
    assert guards["public_promotion_allowed"] is False
    assert guards["measured_alpha_in_two_loop_computation"] is False
