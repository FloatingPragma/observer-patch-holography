"""Tests for the Ward-normalized charged determinant-line lane."""

from __future__ import annotations

import functools

import mpmath as mp
import pytest

import derive_charged_ward_determinant_line as lane


EXPECTED_OPEN_PARENT_IDS = [
    "R_CHARGED_CENTERED_SHAPE",
    "R_WARD_SPACELIKE_ENDPOINT_PAIR",
    "R_WARD_NONLEPTONIC_SUBTRACTION",
    "R_WARD_HIGHER_ORDER_MONOTONICITY",
    "R_ELECTRON_ATOMIC_MASS_TRANSPORT",
    "R_CLEAN_P_ALPHA_ROOT",
    "R_NO_TARGET_CHARGED_DAG",
]


@pytest.fixture(autouse=True)
def _scoped_precision():
    """Set a high mpmath precision for each test and restore it on teardown.

    The tight-tolerance comparisons in this module require extended precision.
    Scoping keeps the global mpmath precision unchanged for other test modules.
    """

    with mp.workdps(120):
        yield


@functools.lru_cache(maxsize=1)
def _artifact() -> dict:
    return lane.build_artifact()


def test_kernel_closed_form_matches_direct_quadrature_to_1e40():
    block = _artifact()["kernel_identity"]
    assert block["all_within_tolerance"] is True
    assert [row["z"] for row in block["rows"]] == ["1e-8", "0.1", "1", "10", "1e6"]
    for row in block["rows"]:
        assert row["within_tolerance"] is True
        assert mp.mpf(row["absolute_residual"]) < mp.mpf("1e-40")
    assert block["series_matches_closed_form"] is True


def test_ward_response_is_strictly_decreasing_with_declared_limits():
    block = _artifact()["monotonicity"]
    assert block["derivative_strictly_negative"] is True
    assert all(mp.mpf(value) < 0 for value in block["derivatives"])
    assert block["derivative_consistent_with_response"] is True
    assert block["upper_limit_is_zero"] is True
    assert block["divergence_lower_bound_holds"] is True
    assert block["high_energy_slope_matches_minus_2_over_pi"] is True
    assert mp.mpf(block["high_energy_slope_residual_at_mu_minus_12"]) < mp.mpf("1e-8")


def test_synthetic_round_trip_recovers_mu_to_1e30():
    block = _artifact()["synthetic_round_trip"]
    assert block["packet_label"] == "synthetic"
    assert block["recovered_within_1e-30"] is True
    assert mp.mpf(block["absolute_recovery_residual"]) < mp.mpf("1e-30")
    assert block["seed_within_asymptotic_window"] is True
    assert block["det_line_recovered"] is True


def test_direct_inversion_functions_recover_the_synthetic_coordinate():
    q = "100"
    shape = ("-2", "0.25", "1.75")
    mu_true = mp.mpf("-3.25")
    response = lane.ward_response(mu_true, q, shape)
    assert response > 0
    recovered = lane.solve_mu(response, q, shape)
    assert abs(recovered - mu_true) < mp.mpf("1e-30")


def test_interval_enclosure_brackets_the_synthetic_root():
    block = _artifact()["interval_enclosure"]
    assert block["certified"] is True
    assert block["encloses_synthetic_mu"] is True
    assert block["width_below_1e-5"] is True
    mu_lo, mu_hi = (mp.mpf(value) for value in block["mu_interval"])
    assert mu_lo <= mp.mpf("-3.25") <= mu_hi


def test_open_parent_gate_lists_the_exact_parents():
    gate = _artifact()["physical_source_packet"]
    assert gate["status"] == "ABSENT"
    assert gate["gate"] == "fail_closed"
    assert [row["id"] for row in gate["open_parents"]] == EXPECTED_OPEN_PARENT_IDS
    assert all(row["reason"] for row in gate["open_parents"])


def test_promotion_flags_are_false_and_status_is_declared():
    artifact = _artifact()
    assert artifact["status"] == "MATHEMATICAL_INVERSION_CLOSED_PHYSICAL_PARENTS_OPEN"
    assert artifact["public_charged_mass_promotion_allowed"] is False
    assert artifact["charged_reference_data_consumed"] is False
    assert (
        artifact["lane_role"]["x2_2_route"]
        == "X2.2 electron-absolute-ratio route B producer"
    )
    assert artifact["checks_pass"] is True
    registry = artifact["theorem_registry"]
    assert "CHARGED_WARD_DETERMINANT_LINE_TO_ELECTRON_SOURCE_RATIO_THEOREM" in registry
    assert "WARD_STRICT_MONOTONICITY_ORBIT_BREAKING_THEOREM" in registry
    assert "WARD_HIGHER_ORDER_STABILITY_CRITERION" in registry
    assert "AUDIT_CORRECTION_COMMON_SHIFT_NO_GO_SCOPE" in registry
