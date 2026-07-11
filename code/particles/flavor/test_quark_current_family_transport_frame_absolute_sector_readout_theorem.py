#!/usr/bin/env python3
"""Tests for the restricted-carrier absolute sector-readout theorem artifact."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_current_family_transport_frame_absolute_sector_readout_theorem import build_artifact


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "particles" / "runs" / "flavor"


def test_current_family_transport_frame_absolute_sector_readout_theorem_closes() -> None:
    physical_sigma = json.loads(
        (RUNS / "quark_current_family_transport_frame_physical_sigma_lift_theorem.json").read_text(encoding="utf-8")
    )
    sigma_target = json.loads((RUNS / "quark_current_family_exact_sigma_target.json").read_text(encoding="utf-8"))
    absolute_collapse = json.loads((RUNS / "quark_absolute_readout_algebraic_collapse.json").read_text(encoding="utf-8"))
    payload = build_artifact(physical_sigma, sigma_target, absolute_collapse)

    assert payload["artifact"] == "oph_quark_current_family_transport_frame_absolute_sector_readout_theorem"
    assert payload["proof_status"] == "closed_current_family_transport_frame_absolute_sector_readout_theorem"
    assert payload["corresponds_to_global_contract"]["id"] == "quark_absolute_sector_readout_theorem"
    assert payload["theorem_scope"] == "current_family_common_refinement_transport_frame_only"
    scales = payload["emitted_absolute_sector_scales"]
    assert abs(float(scales["g_u"]) - 0.7792490879740623) < 1.0e-12
    assert abs(float(scales["g_d"]) - 0.12226559489803177) < 1.0e-12
