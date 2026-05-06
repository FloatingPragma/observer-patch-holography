#!/usr/bin/env python3
"""Smoke tests for the constructive Thomson endpoint contract."""

from __future__ import annotations

from thomson_endpoint_contract import build_contract


def test_contract_requires_constructive_worker_outputs() -> None:
    payload = build_contract()

    assert payload["artifact"] == "oph_ward_projected_thomson_endpoint_contract"
    assert payload["github_issue"] == 235
    assert payload["closed_blocker_isolation_issue"] == 223
    assert payload["promotion_allowed"] is False
    assert payload["computed_package"] == "code/P_derivation/runtime/thomson_endpoint_package_current.json"
    assert payload["worker_result_policy"]["obstruction_only_result_allowed"] is False
    object_ids = {entry["id"] for entry in payload["constructive_objects"]}
    assert "rho_had_spectral_measure" in object_ids
    assert "delta_qcd_screening_and_endpoint_remainder" in object_ids
    assert "full_endpoint_interval_certificate" in object_ids
    assert "measured_alpha_0" in payload["forbidden_solver_inputs"]
