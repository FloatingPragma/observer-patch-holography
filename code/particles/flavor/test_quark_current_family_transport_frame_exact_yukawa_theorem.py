#!/usr/bin/env python3
"""Tests for the exact current-family transport-frame Yukawa theorem artifact."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_current_family_transport_frame_exact_yukawa_theorem import build_artifact


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "particles" / "runs" / "flavor"


def _load(name: str) -> dict:
    return json.loads((RUNS / name).read_text(encoding="utf-8"))


def test_current_family_transport_frame_exact_yukawa_theorem_closes() -> None:
    payload = build_artifact(
        _load("quark_current_family_transport_frame_exact_forward_yukawas.json"),
        _load("quark_current_family_end_to_end_exact_pdg_derivation_chain.json"),
    )

    assert payload["proof_status"] == "closed_target_anchored_mixed_convention_mass_texture_audit"
    assert payload["target_name"] == "mixed_convention_quark_mass_textures_on_declared_current_family_transport_frame"
    assert payload["forward_yukawa_artifact"]["artifact"] == "oph_quark_current_family_transport_frame_exact_forward_yukawas"
    assert payload["forward_yukawa_artifact"]["forward_certified"] is True
    assert payload["forward_yukawa_artifact"]["physical_yukawa_certified"] is False
    assert payload["minimal_exact_blocker_set"] == ["QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE"]
