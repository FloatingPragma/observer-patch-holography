#!/usr/bin/env python3
"""Smoke tests for the conditional Thomson endpoint package."""

from __future__ import annotations

from decimal import Decimal

from thomson_endpoint_package import build_endpoint_package


def test_endpoint_package_is_conditional_and_source_guarded() -> None:
    payload = build_endpoint_package(
        {
            "mode": "synthetic",
            "precision": 10,
            "alpha_inv": "10",
            "source_anchor_alpha_inv": "8",
            "p": "1.2",
            "phi": "1",
            "sqrt_pi": "2",
            "structured_running": {"total_delta_alpha_inv": "2"},
        },
        compare_alpha_inv=Decimal("11"),
        compare_alpha_inv_uncertainty=Decimal("0.1"),
        precision=10,
        su2_cutoff=6,
        su3_cutoff=4,
    )

    assert payload["artifact"] == "oph_ward_projected_thomson_endpoint_package"
    assert payload["github_issue"] == 223
    assert payload["successor_github_issue"] == 235
    assert payload["claim_status"] == "endpoint_package_computed_blocker_isolated"
    assert payload["promotion_allowed"] is False
    assert payload["source_only_guard"]["codata_enters_solver"] is False
    assert payload["source_only_guard"]["measured_endpoint_allowed_as_transport_input"] is False
    assert payload["implemented_fixed_point_gap"]["missing_transport_delta_alpha_inv"] == "1"
    assert (
        payload["first_non_internalized_object"]["id"]
        == "ward_projected_qcd_screening_and_endpoint_remainder"
    )
    assert payload["first_non_internalized_object"]["successor_github_issue"] == 235
    assert payload["issue_223_acceptance"]["theorem_grade_object_defined"] is True
    assert payload["issue_223_acceptance"]["closable_as_measured_alpha_derivation"] is False
    assert payload["issue_223_acceptance"]["closable_as_blocker_isolation_package"] is True
    assert payload["issue_223_acceptance"]["successor_issue_for_source_residual_map"] == 235
