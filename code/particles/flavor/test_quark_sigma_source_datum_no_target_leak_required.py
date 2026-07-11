#!/usr/bin/env python3
"""Tests for the missing quark sigma source theorem gate."""

from __future__ import annotations

import json
from pathlib import Path

from derive_quark_sigma_source_datum_no_target_leak_required import build_artifact


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "particles" / "runs" / "flavor"


def _load(name: str) -> dict:
    return json.loads((RUNS / name).read_text(encoding="utf-8"))


def test_quark_sigma_source_gate_records_target_free_missing_theorem() -> None:
    payload = build_artifact(
        _load("quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"),
        _load("quark_edge_statistics_spread_candidate.json"),
    )

    assert payload["artifact"] == "oph_quark_sigma_source_datum_no_target_leak_required"
    assert payload["status"] == "missing_theorem"
    assert payload["claim_tier"] == "selected_class_conditional_on_source_sigma"
    assert payload["source_only_sigma_emitted"] is False
    assert payload["downstream_algebra_closed"] is True
    assert payload["selected_class"] == "f_P"
    assert payload["promotion_allowed"] is False
    assert "QUARK_SIGMA_SOURCE_SELECTOR" in payload["missing_for_promotion"]
    assert "NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE" in payload["missing_for_promotion"]
    assert "quark_current_family_exact_sigma_target" in payload["forbidden_ancestors"]
    assert "PDG" in payload["forbidden_ancestors"]

    target = payload["target_values_for_future_source_theorem"]
    assert abs(target["sigma_u"] - 5.579692209267639) < 1.0e-12
    assert abs(target["sigma_d"] - 3.300314452061615) < 1.0e-12
    assert abs(target["sigma_seed_ud"] - 4.440003330664627) < 1.0e-12
    assert abs(target["eta_ud"] - 1.139688878603012) < 1.0e-12

    candidate = payload["strongest_current_source_candidate"]
    assert abs(candidate["sigma_u_edge"] - 5.578418804072826) < 1.0e-12
    assert abs(candidate["sigma_d_edge"] - 3.4210589139721543) < 1.0e-12
    assert abs(candidate["required_R_u"] - (-0.004490377677282886)) < 1.0e-12
    assert abs(candidate["required_R_d"] - (-0.1247947151634663)) < 1.0e-12
    assert abs(candidate["required_R_seed"] - (-0.06464254642037481)) < 1.0e-12
    assert abs(candidate["required_R_eta"] - 0.060152168743091705) < 1.0e-12
