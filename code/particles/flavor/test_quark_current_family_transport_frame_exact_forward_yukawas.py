#!/usr/bin/env python3
"""Tests for legacy-named target-audit mass textures on the transport frame."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_current_family_transport_frame_exact_forward_yukawas import build_artifact


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "particles" / "runs" / "flavor"


def _load(name: str) -> dict:
    return json.loads((RUNS / name).read_text(encoding="utf-8"))


def test_current_family_transport_frame_exact_forward_yukawas_close() -> None:
    payload = build_artifact(
        _load("quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"),
        _load("quark_current_family_transport_frame_exact_pdg_completion.json"),
    )

    assert payload["proof_status"] == "closed_target_anchored_mixed_convention_mass_texture_audit"
    assert payload["scope"] == "current_family_common_refinement_transport_frame_only"
    assert payload["forward_certified"] is True
    assert payload["certification_status"] == "forward_matrix_certified"
    assert payload["matrix_classification"] == "mixed_scheme_GeV_mass_texture_matrices"
    assert payload["physical_yukawa_certified"] is False
    assert payload["single_common_scale_running_sextet"] is False
    assert payload["promotion_blockers"] == ["QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE"]
    assert abs(float(payload["singular_values_u"][0]) - 0.00216) < 1.0e-8
    assert abs(float(payload["singular_values_d"][0]) - 0.0047) < 1.0e-8
