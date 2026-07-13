"""Tests for the up-type register exponent scan."""

from __future__ import annotations

import functools

import derive_up_type_register_exponent_scan as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_scan_is_negative_and_fail_closed():
    artifact = _artifact()
    assert artifact["status"] == "SCAN_NEGATIVE_FAMILY_REMOVED"
    assert artifact["promotion_allowed"] is False
    assert not any(row["certified"] for row in artifact["scan_rows"].values())


def test_all_frozen_bases_report_both_residuals():
    rows = _artifact()["scan_rows"]
    assert set(rows) == {
        "sqrt_alpha_U",
        "exp_minus_three_halves",
        "alpha_U_times_P",
        "exp_minus_pi_over_two",
    }
    assert all(len(row["integer_residuals"]) == 2 for row in rows.values())


def test_cross_register_observations_are_labeled_post_exposure():
    cross = _artifact()["cross_register_observations"]
    assert cross["classification"] == "post_exposure_observation_not_certified"
    assert "trial" in cross["trials_ledger"]


def test_compare_only_inputs_are_declared():
    artifact = _artifact()
    assert artifact["measured_compare_only_inputs"]["mc_mc_gev"] == 1.27
    assert (
        artifact["row_class"]
        == "compare_only_candidate_scan_never_a_prediction_ancestor"
    )
