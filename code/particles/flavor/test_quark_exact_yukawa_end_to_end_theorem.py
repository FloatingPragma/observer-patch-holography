#!/usr/bin/env python3
"""Tests for the end-to-end exact Yukawa theorem artifact."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_exact_yukawa_end_to_end_theorem import build_artifact


ROOT = Path(__file__).resolve().parents[2]
YUKAWA_THEOREM_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_yukawa_theorem.json"


def test_quark_exact_yukawa_end_to_end_theorem_closes() -> None:
    theorem = json.loads(YUKAWA_THEOREM_JSON.read_text(encoding="utf-8"))
    payload = build_artifact(theorem)

    assert payload["proof_status"] == "closed_target_anchored_mixed_convention_mass_texture_audit"
    assert payload["target_name"] == "mixed_convention_quark_mass_textures_on_declared_current_family_transport_frame"
    assert payload["supporting_theorem_artifact"] == "oph_quark_current_family_transport_frame_exact_yukawa_theorem"
    assert payload["forward_yukawa_artifact"]["forward_certified"] is True
    assert payload["forward_yukawa_artifact"]["physical_yukawa_certified"] is False
