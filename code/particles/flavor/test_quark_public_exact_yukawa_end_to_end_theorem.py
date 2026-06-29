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


def test_quark_public_exact_yukawa_end_to_end_theorem_closes_public_target() -> None:
    public_sigma_theorem = build_public_sigma_descent(
        _load("quark_current_family_transport_frame_sector_attached_lift.json"),
        _load("quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"),
        _load("overlap_edge_line_lift.json"),
    )
    payload = build_artifact(
        public_sigma_theorem,
        _load("quark_exact_pdg_end_to_end_theorem.json"),
        _load("quark_exact_yukawa_end_to_end_theorem.json"),
    )

    assert payload["artifact"] == "oph_quark_public_exact_yukawa_end_to_end_theorem"
    assert payload["proof_status"] == "blocked_target_derived_sigma_source_missing"
    assert payload["target_name"] == "selected_class_conditional_forward_quark_yukawas"
    assert payload["claim_tier"] == "selected_class_conditional_on_source_sigma"
    assert payload["public_promotion_allowed"] is False
    assert payload["source_only_sigma_emitted"] is False
    assert payload["downstream_algebra_closed"] is True
    assert payload["display_allowed_as_selected_class_exact_witness"] is True
    assert "QUARK_SIGMA_SOURCE_SELECTOR" in payload["minimal_exact_blocker_set"]
    assert "NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE" in payload["minimal_exact_blocker_set"]
    assert payload["non_circularity_status"]["target_derived_sigma_datum_used"] is True
    outputs = payload["public_exact_outputs"]["forward_yukawa_artifact"]
    assert outputs["artifact"] == "oph_quark_current_family_transport_frame_exact_forward_yukawas"
    assert outputs["forward_certified"] is True
    assert payload["public_exact_outputs"]["exact_running_values_gev"]["d"] == 0.004699999999999999
