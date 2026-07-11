#!/usr/bin/env python3
"""Tests for the strengthened restricted-carrier physical sigma-lift theorem."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem import build_artifact


ROOT = Path(__file__).resolve().parents[2]
PHYSICAL_SIGMA_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_physical_sigma_lift_theorem.json"
)
SIGMA_TARGET_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_sigma_target.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_strengthened_physical_sigma_lift_merges_sigma_target() -> None:
    payload = build_artifact(_load(PHYSICAL_SIGMA_THEOREM_JSON), _load(SIGMA_TARGET_JSON))

    assert payload["proof_status"] == "closed_current_family_transport_frame_strengthened_physical_sigma_lift_theorem"
    assert payload["theorem_scope"] == "current_family_common_refinement_transport_frame_only"
    assert payload["compressed_global_contract"]["id"] == "strengthened_quark_physical_sigma_ud_lift"
    assert payload["restricted_sigma_ud_phys_element"]["sigma_id"] == "sigma_phys_transport_frame_current_family"
    sigma = payload["theorem_grade_physical_sigma_datum"]
    assert abs(float(sigma["sigma_u"]) - 5.579692209267639) < 1.0e-12
    assert abs(float(sigma["sigma_d"]) - 3.300314452061615) < 1.0e-12
