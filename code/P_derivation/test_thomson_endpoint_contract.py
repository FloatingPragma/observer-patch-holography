#!/usr/bin/env python3
"""Smoke tests for the constructive Thomson endpoint contract."""

from __future__ import annotations

from thomson_endpoint_contract import build_contract


def test_contract_records_scientific_boundary_without_project_state() -> None:
    payload = build_contract()

    assert payload["artifact"] == "oph_ward_projected_thomson_endpoint_contract"
    assert payload["promotion_allowed"] is False
    assert payload["status"] == "source_spectral_reduction_closed_measure_payload_absent"
    assert payload["promotion_boundary"]["minimal_required_theorem"] == (
        "WardProjectedHadronicSpectralEmission_Q"
    )
    assert set(payload) == {
        "artifact",
        "blocking_artifacts",
        "computed_package",
        "constructive_objects",
        "endpoint_formula",
        "forbidden_solver_inputs",
        "no_go_results",
        "promotion_allowed",
        "promotion_boundary",
        "status",
    }
    assert payload["computed_package"] == "code/P_derivation/runtime/thomson_endpoint_package_current.json"
    assert payload["blocking_artifacts"]["screening_invariant_no_go"].endswith(
        "screening_invariant_no_go_current.json"
    )
    assert payload["blocking_artifacts"]["source_spectral_theorem"].endswith(
        "source_spectral_theorem_current.json"
    )
    assert payload["no_go_results"]["detuning_only_bypass"] == "closed_no_go"
    assert "hadronic measure" in payload["promotion_boundary"]["reason"]
    object_ids = {entry["id"] for entry in payload["constructive_objects"]}
    assert "rho_had_spectral_measure" in object_ids
    assert "source_spectral_reduction_theorem" in object_ids
    assert "screening_invariant_no_go" in object_ids
    assert "delta_qcd_screening_and_endpoint_remainder" in object_ids
    assert "full_endpoint_interval_certificate" in object_ids
    assert "measured_alpha_0" in payload["forbidden_solver_inputs"]
    assert "c_Q_target" in payload["forbidden_solver_inputs"]
