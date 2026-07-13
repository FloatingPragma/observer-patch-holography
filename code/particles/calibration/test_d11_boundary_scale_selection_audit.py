"""Tests for the boundary-scale selection audit."""

from __future__ import annotations

import functools

import derive_d11_boundary_scale_selection_audit as lane


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build()


def test_flow_internal_no_go_is_certified_with_margin():
    no_go = _artifact()["flow_internal_selection_no_go"]
    assert no_go["no_root_certified"] is True
    assert no_go["minimum_flow_derivative"] > 1.0e-5
    assert all(
        row["flow_dbeta_lambda_dt"] > 0 for row in no_go["sampled_rows"]
    )


def test_registry_is_frozen_with_classifications():
    registry = _artifact()["candidate_registry_frozen"]
    assert set(registry) == {
        "mu_U_gauge_unification",
        "log_midpoint_half_turn",
        "E_cell_pixel_energy",
        "E_star",
    }
    assert (
        registry["log_midpoint_half_turn"]["classification"]
        == "prospective_post_exposure_candidate"
    )
    assert (
        registry["mu_U_gauge_unification"]["classification"]
        == "archived_declared_choice"
    )


def test_log_midpoint_two_loop_values_are_recorded():
    entry = _artifact()["candidate_registry_frozen"]["log_midpoint_half_turn"]
    assert abs(entry["two_loop"]["mt_pole_gev"] - 172.63) < 0.05
    assert abs(entry["two_loop"]["mh_tree_gev"] - 125.77) < 0.05


def test_selection_stays_open_and_fail_closed():
    artifact = _artifact()
    assert artifact["status"] == "SELECTION_OPEN_CANDIDATES_FROZEN_PROSPECTIVELY"
    assert artifact["promotion_allowed"] is False
    assert (
        artifact["registered_discriminating_test"]["object"]
        == "FROZEN_THREE_LOOP_RG_MATCHING_PACKET"
    )
    assert (
        artifact["adopted_working_conditional_branch"]["name"]
        == "log_midpoint_half_turn"
    )
