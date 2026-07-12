#!/usr/bin/env python3
"""Tests for the public exact Yukawa theorem wrapper."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_public_physical_sigma_datum_descent import build_artifact as build_public_sigma_descent
from derive_quark_public_exact_yukawa_end_to_end_theorem import build_artifact


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "particles" / "runs" / "flavor"


def _load(name: str) -> dict:
    return json.loads((RUNS / name).read_text(encoding="utf-8"))


def test_quark_public_exact_yukawa_end_to_end_theorem_fails_closed_on_surface_mismatch() -> None:
    public_sigma_theorem = build_public_sigma_descent(
        _load("quark_current_family_transport_frame_sector_attached_lift.json"),
        _load("quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"),
        _load("overlap_edge_line_lift.json"),
    )
    exact_pdg = _load("quark_exact_pdg_end_to_end_theorem.json")
    exact_pdg["exact_running_values_gev"]["s"] *= 1.01
    exact_pdg["exact_running_values_gev"]["t"] *= 1.01
    payload = build_artifact(
        public_sigma_theorem,
        exact_pdg,
        _load("quark_exact_yukawa_end_to_end_theorem.json"),
        _load("quark_running_mass_scheme_convention_obstruction.json"),
    )

    assert payload["artifact"] == "oph_quark_public_exact_yukawa_end_to_end_theorem"
    assert payload["proof_status"] == "blocked_inconsistent_exact_mass_yukawa_surfaces"
    assert payload["target_name"] == "inconsistent_selected_class_mass_yukawa_support_surfaces"
    assert payload["claim_tier"] == "inconsistent_exact_support_surfaces"
    assert payload["public_promotion_allowed"] is False
    assert payload["source_only_sigma_emitted"] is False
    assert payload["downstream_algebra_closed"] is False
    assert payload["mass_texture_algebra_closed"] is False
    assert payload["physical_yukawa_construction_closed"] is False
    assert payload["display_allowed_as_selected_class_exact_witness"] is False
    assert "QUARK_SIGMA_SOURCE_SELECTOR" in payload["minimal_exact_blocker_set"]
    assert "NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE" in payload["minimal_exact_blocker_set"]
    assert "QUARK_EXACT_MASS_YUKAWA_SURFACE_CONSISTENCY" in payload["minimal_exact_blocker_set"]
    assert "QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE" in payload["minimal_exact_blocker_set"]
    assert payload["non_circularity_status"]["target_derived_sigma_datum_used"] is True
    consistency = payload["mass_yukawa_consistency"]
    assert consistency["status"] == "mismatched_mass_sextet_and_forward_yukawa_singular_values"
    assert consistency["consistent"] is False
    assert consistency["single_exact_sextet_matrix_pair_claim_allowed"] is False
    assert consistency["comparisons"]["u"]["consistent_within_tolerance"] is True
    assert consistency["comparisons"]["s"]["consistent_within_tolerance"] is False
    assert consistency["comparisons"]["t"]["consistent_within_tolerance"] is False
    assert payload["public_exact_outputs"]["single_exact_mass_texture_pair_claim_allowed"] is False
    assert payload["public_exact_outputs"]["single_exact_physical_yukawa_pair_claim_allowed"] is False
    outputs = payload["public_exact_outputs"]["forward_yukawa_artifact"]
    assert outputs["artifact"] == "oph_quark_current_family_transport_frame_exact_forward_yukawas"
    assert outputs["forward_certified"] is True
    assert payload["public_exact_outputs"]["exact_running_values_gev"]["d"] == 0.004699999999999999


def test_quark_public_exact_yukawa_consistency_accepts_one_common_sextet() -> None:
    public_sigma_theorem = build_public_sigma_descent(
        _load("quark_current_family_transport_frame_sector_attached_lift.json"),
        _load("quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"),
        _load("overlap_edge_line_lift.json"),
    )
    exact_pdg = _load("quark_exact_pdg_end_to_end_theorem.json")
    exact_yukawa = _load("quark_exact_yukawa_end_to_end_theorem.json")
    forward = exact_yukawa["forward_yukawa_artifact"]
    exact_pdg["exact_running_values_gev"] = {
        "u": forward["singular_values_u"][0],
        "c": forward["singular_values_u"][1],
        "t": forward["singular_values_u"][2],
        "d": forward["singular_values_d"][0],
        "s": forward["singular_values_d"][1],
        "b": forward["singular_values_d"][2],
    }

    payload = build_artifact(
        public_sigma_theorem,
        exact_pdg,
        exact_yukawa,
        _load("quark_running_mass_scheme_convention_obstruction.json"),
    )

    assert payload["mass_yukawa_consistency"]["consistent"] is True
    assert payload["public_exact_outputs"]["single_exact_mass_texture_pair_claim_allowed"] is True
    assert payload["public_exact_outputs"]["single_exact_physical_yukawa_pair_claim_allowed"] is False
    assert payload["proof_status"] == (
        "blocked_mixed_scheme_dimensionful_mass_textures_not_physical_yukawas"
    )
    assert payload["claim_tier"] == "mixed_scheme_dimensionful_mass_texture_audit"
    assert payload["mass_texture_algebra_closed"] is True
    assert payload["physical_yukawa_construction_closed"] is False
    assert "QUARK_EXACT_MASS_YUKAWA_SURFACE_CONSISTENCY" not in payload["minimal_exact_blocker_set"]
    assert "QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE" in payload["minimal_exact_blocker_set"]
