#!/usr/bin/env python3
"""Tests for the merged current-family transport-frame exact completion artifact."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_current_family_transport_frame_exact_pdg_completion import build_artifact


ROOT = Path(__file__).resolve().parents[2]
STRENGTHENED_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"
)
ABS_COLLAPSE_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_absolute_readout_algebraic_collapse.json"
PDG_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_pdg_theorem.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_transport_frame_exact_pdg_completion_closes_declared_surface() -> None:
    payload = build_artifact(_load(STRENGTHENED_THEOREM_JSON), _load(ABS_COLLAPSE_JSON), _load(PDG_JSON))

    assert payload["proof_status"] == "closed_target_anchored_mixed_convention_coordinate_audit"
    assert payload["target_anchored"] is True
    assert payload["source_only_prediction"] is False
    assert payload["single_running_quark_sextet_claim_allowed"] is False
    assert payload["theorem_scope"] == "current_family_common_refinement_transport_frame_only"
    assert payload["strengthened_physical_sigma_lift"]["artifact"] == "oph_quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem"
    masses = payload["exact_running_values_gev"]
    assert abs(float(masses["u"]) - 0.00216) < 1.0e-12
    assert abs(float(masses["d"]) - 0.0047) < 1.0e-12
    assert abs(float(masses["s"]) - 0.0929) < 1.0e-12
    assert abs(float(masses["c"]) - 1.2729) < 1.0e-12
    assert abs(float(masses["b"]) - 4.186) < 1.0e-12
