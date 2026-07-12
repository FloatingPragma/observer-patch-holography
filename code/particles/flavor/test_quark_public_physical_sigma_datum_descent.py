#!/usr/bin/env python3
"""Tests for the direct public sigma-datum descent theorem."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_public_physical_sigma_datum_descent import build_artifact


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "particles" / "runs" / "flavor"


def _load(name: str) -> dict:
    return json.loads((RUNS / name).read_text(encoding="utf-8"))


def test_quark_public_physical_sigma_datum_descent_closes_selected_public_class() -> None:
    payload = build_artifact(
        _load("quark_current_family_transport_frame_sector_attached_lift.json"),
        _load("quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"),
        _load("overlap_edge_line_lift.json"),
    )

    assert payload["artifact"] == "oph_quark_public_physical_sigma_datum_descent"
    assert payload["proof_status"] == "blocked_target_derived_sigma_datum_descent"
    assert payload["theorem_id"] == "selected_bridge_fiber_sigma_descent_not_source_selection"
    assert payload["claim_tier"] == "selected_class_conditional_on_source_sigma"
    assert payload["public_promotion_allowed"] is False
    assert payload["source_only_sigma_emitted"] is False
    assert payload["downstream_algebra_closed"] is True
    assert payload["sigma_descent_is_sigma_selection"] is False
    assert payload["display_allowed_as_selected_class_witness"] is True
    assert payload["non_circularity_status"]["target_derived_sigma_datum_used"] is True
    assert payload["non_circularity_status"]["source_sigma_selector_closed"] is False
    assert payload["non_circularity_status"]["missing_source_object"] == (
        "quark_sigma_source_datum_no_target_leak_required"
    )
    assert "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM" in payload["missing_for_promotion"]
    assert payload["induces_global_contract"]["id"] == "strengthened_quark_physical_sigma_ud_lift"
    assert payload["realized_transport_frame_section_uniqueness"]["common_refinement_level"] == 1
    assert payload["realized_transport_frame_section_uniqueness"]["refinement_functoriality_closed"] is True
    assert payload["realized_transport_frame_section_uniqueness"]["common_refinement_invariance_closed_on_current_family"] is True
    assert payload["declared_bridge_fiber"]["name"] == "R_decl(f_P)"
    assert payload["declared_bridge_fiber"]["realized_section_is_singleton_mod_diagonal_phase"] is True
    assert payload["declared_bridge_fiber_invariance_theorem"]["id"] == "declared_selected_public_bridge_fiber_sigma_constancy"
    assert "For any r,r' in R_decl(f_P)" in payload["declared_bridge_fiber_invariance_theorem"]["statement"]
    sigma = payload["descended_physical_sigma_datum"]
    assert abs(float(sigma["sigma_u"]) - 5.579692209267639) < 1.0e-12
    assert abs(float(sigma["sigma_d"]) - 3.300314452061615) < 1.0e-12
    assert payload["selected_public_physical_frame_class"]["selected_by"] == "P"
