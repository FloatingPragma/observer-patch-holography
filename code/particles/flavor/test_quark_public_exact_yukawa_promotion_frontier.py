#!/usr/bin/env python3
"""Tests for the final public exact-Yukawa promotion frontier."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_public_physical_sigma_datum_descent import build_artifact as build_public_sigma_descent
from derive_quark_public_exact_yukawa_end_to_end_theorem import build_artifact as build_public_exact_yukawa
from derive_quark_public_exact_yukawa_promotion_frontier import build_artifact


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "particles" / "runs" / "flavor"


def _load(name: str) -> dict:
    return json.loads((RUNS / name).read_text(encoding="utf-8"))


def test_quark_public_exact_yukawa_promotion_frontier_records_both_public_burdens() -> None:
    public_sigma_theorem = build_public_sigma_descent(
        _load("quark_current_family_transport_frame_sector_attached_lift.json"),
        _load("quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"),
        _load("overlap_edge_line_lift.json"),
    )
    public_exact_yukawa_theorem = build_public_exact_yukawa(
        public_sigma_theorem,
        _load("quark_exact_pdg_end_to_end_theorem.json"),
        _load("quark_exact_yukawa_end_to_end_theorem.json"),
    )
    payload = build_artifact(
        _load("quark_public_strengthened_physical_sigma_lift_frontier.json"),
        public_sigma_theorem,
        _load("quark_exact_yukawa_end_to_end_theorem.json"),
        public_exact_yukawa_theorem,
    )

    assert payload["artifact"] == "oph_quark_public_exact_yukawa_promotion_frontier"
    assert payload["proof_status"] == (
        "blocked_mixed_scheme_dimensionful_mass_textures_not_physical_yukawas"
    )
    assert payload["target_name"] == "selected_class_mixed_scheme_mass_texture_audit"
    assert payload["public_promotion_allowed"] is False
    assert payload["non_circularity_status"]["missing_source_object"] == (
        "quark_sigma_source_datum_no_target_leak_required"
    )
    assert payload["resolved_by_theorem_artifact"] == "oph_quark_public_exact_yukawa_end_to_end_theorem"
    assert payload["final_public_theorem_candidate"]["id"] == "selected_bridge_fiber_sigma_descent_not_source_selection"
    assert payload["alternate_upstream_route"]["id"] == "oph_generation_bundle_branch_generator_splitting"
    assert payload["closed_local_endpoint"]["artifact"] == "oph_quark_exact_yukawa_end_to_end_theorem"
    assert payload["closed_local_endpoint"]["certified_physical_yukawa_matrices"] is False
    assert "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM" in payload[
        "minimal_exact_blocker_set"
    ]
    assert "QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE" in payload[
        "minimal_exact_blocker_set"
    ]
    induced = payload["closed_public_endpoint"]["public_exact_outputs"]["forward_yukawa_artifact"]
    assert induced["artifact"] == "oph_quark_current_family_transport_frame_exact_forward_yukawas"
    assert induced["forward_certified"] is True
    assert induced["certification_status"] == "forward_matrix_certified"
