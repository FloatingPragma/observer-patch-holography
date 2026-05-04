#!/usr/bin/env python3
"""Smoke tests for the constructive Thomson endpoint contract."""

from __future__ import annotations

from thomson_endpoint_contract import build_contract


def test_contract_requires_constructive_worker_outputs() -> None:
    payload = build_contract()

    assert payload["artifact"] == "oph_ward_projected_thomson_endpoint_contract"
    assert payload["promotion_allowed"] is False
    assert payload["worker_result_policy"]["obstruction_only_result_allowed"] is False
    object_ids = {entry["id"] for entry in payload["constructive_objects"]}
    assert "rho_had_spectral_measure" in object_ids
    assert "full_endpoint_interval_certificate" in object_ids
    assert "measured_alpha_0" in payload["forbidden_solver_inputs"]
