#!/usr/bin/env python3
"""Tests for the compatibility projection of the quark-spread obstruction."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from derive_quark_sigma_source_datum_no_target_leak_required import build_artifact


ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "particles" / "runs" / "flavor"


def _load(name: str) -> dict[str, Any]:
    return json.loads((RUNS / name).read_text(encoding="utf-8"))


def _walk_keys(value: Any) -> set[str]:
    if isinstance(value, dict):
        keys = set(value)
        for child in value.values():
            keys.update(_walk_keys(child))
        return keys
    if isinstance(value, list):
        keys: set[str] = set()
        for child in value:
            keys.update(_walk_keys(child))
        return keys
    return set()


def test_quark_sigma_source_gate_projects_the_target_free_obstruction() -> None:
    payload = build_artifact(_load("quark_sigma_source_nonidentifiability_obstruction.json"))

    assert payload["artifact"] == "oph_quark_sigma_source_datum_no_target_leak_required"
    assert payload["status"] == "closed_current_corpus_nonidentifiability_obstruction"
    assert payload["closure_kind"] == "theorem_grade_sharper_obstruction"
    assert payload["claim_tier"] == "selected_class_conditional_support_source_sigma_nonidentifiable"
    assert payload["source_only_sigma_emitted"] is False
    assert payload["downstream_algebra_closed"] is True
    assert payload["promotion_allowed"] is False
    assert payload["resolved_github_issues"] == [377, 379, 380]
    assert payload["compatible_spread_fiber"] == "(R_{>0})^2"
    assert payload["compatible_spread_fiber_dimension"] == 2
    assert payload["independent_unselected_coordinates"] == ["sigma_u", "sigma_d"]
    assert payload["affine_downstream_injective"] is True
    assert "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM" in payload["missing_for_promotion"]
    assert "NO_TARGET_LEAK_DAG_QUARK_SOURCE_SPREAD" in payload["missing_for_promotion"]
    assert payload["dependency_audit"]["no_target_leak"] is True
    assert payload["dependency_audit"]["loads_running_quark_reference_rows"] is False

    keys = _walk_keys(payload)
    assert "target_values_for_future_source_theorem" not in keys
    assert "strongest_current_source_candidate" not in keys
    assert "required_R_u" not in keys
    assert "required_R_d" not in keys
