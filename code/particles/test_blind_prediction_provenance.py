#!/usr/bin/env python3
"""Smoke tests for the quantitative particle provenance ledger."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from build_blind_prediction_provenance import build_payload  # noqa: E402


def test_blind_prediction_provenance_records_target_use_and_declared_sensitivity_taxonomy() -> None:
    payload = build_payload()

    assert payload["artifact"] == "oph_blind_prediction_provenance_ledger"
    assert payload["promotion_allowed"] is False
    assert payload["scientific_boundary"]["public_rows_classified"] is True
    assert payload["scientific_boundary"]["numeric_sensitivity_intervals_supplied"] is False
    assert "pipeline_closure_status" not in payload["source_surfaces"]
    assert "github_issue" not in payload
    assert "closure_gate" not in payload
    assert "finalization_gates" not in payload
    assert payload["convention_sensitivity"]["classification"] == (
        "declared_taxonomy__numeric_sweep_not_performed"
    )
    row_map = {row["particle_id"]: row for row in payload["rows"]}
    carrier_map = {row["carrier_id"]: row for row in payload["carrier_mode_rows"]}
    withheld_map = {row["particle_id"]: row for row in payload["withheld_rows"]}
    assert "photon" not in row_map
    assert carrier_map["photon"]["blind_status"] == "not_a_quantum_particle_prediction"
    assert carrier_map["photon"]["particle_promotion_allowed"] is False
    assert carrier_map["gluon"]["quantum_particle_gate"] == "not_passed"
    assert "w_boson" not in row_map
    assert "z_boson" not in row_map
    assert "electron" not in row_map
    assert withheld_map["electron"]["target_use"] == "target_values_or_target_derived_datum_used"
    assert withheld_map["electron"]["blind_status"] == "withheld_not_blind"
    assert row_map["higgs"]["blind_status"] == "conditionally_blind_on_declared_surface"
    assert "top_quark" not in row_map
    assert withheld_map["top_quark"]["target_use"] == "target_values_or_target_derived_datum_used"
    assert "electron_neutrino" not in row_map
    assert withheld_map["electron_neutrino"]["target_use"] == "target_ranked_selector_development_and_correlated_profile_rejection"
    assert withheld_map["electron_neutrino"]["blind_status"] == "withheld_not_blind_rejected_candidate"
    workflows = {workflow["id"]: workflow for workflow in payload["prospective_comparison_workflows"]}
    assert workflows["new_quantity_pre_reference_provenance"]["classification"] == (
        "protocol_defined__not_exercised"
    )
    assert workflows["convention_sensitivity_sweep"]["classification"] == (
        "taxonomy_declared__certified_intervals_required"
    )
