#!/usr/bin/env python3
"""Validate the direct-top codomain closure certificate."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_direct_top_bridge_contract.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "direct_top_bridge_contract.json"
REFERENCE = ROOT / "particles" / "data" / "particle_reference_values.json"


def test_direct_top_bridge_contract_records_auxiliary_row_no_go() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    references = json.loads(REFERENCE.read_text(encoding="utf-8"))["entries"]
    expected_delta = (
        float(references["top_quark_direct_aux"]["value_gev"])
        - float(references["top_quark"]["value_gev"])
    )
    assert payload["artifact"] == "oph_direct_top_bridge_contract"
    assert payload["status"] == "hard_no_go_current_corpus_compare_only_direct_top_codomain"
    assert payload["promotion_allowed"] is False
    assert set(payload) == {
        "artifact",
        "auxiliary_direct_top_coordinate",
        "comparison_only_readout",
        "constructive_objects",
        "current_target_audit_coordinate",
        "forbidden_solver_inputs",
        "formal_nonidentifiability_witness",
        "promotion_allowed",
        "promotion_boundary",
        "scientific_result",
        "status",
    }
    assert payload["scientific_result"]["result_kind"] == "hard_no_go_current_corpus"
    assert payload["scientific_result"]["auxiliary_row_policy"] == "compare_only_not_promotable"
    assert payload["current_target_audit_coordinate"]["pdg_summary_id"] == "Q007TP4"
    assert payload["auxiliary_direct_top_coordinate"]["pdg_summary_id"] == "Q007TP"
    assert payload["current_target_audit_coordinate"]["value_gev"] == pytest.approx(
        references["top_quark"]["value_gev"]
    )
    assert payload["auxiliary_direct_top_coordinate"]["value_gev"] == pytest.approx(
        references["top_quark_direct_aux"]["value_gev"]
    )
    assert payload["comparison_only_readout"]["direct_minus_current_coordinate_gev"] == pytest.approx(
        expected_delta,
        abs=1.0e-12,
    )
    assert payload["comparison_only_readout"]["within_combined_one_sigma"] is True
    assert payload["promotion_boundary"]["promotion_allowed"] is False
    assert payload["formal_nonidentifiability_witness"]["lambda_matching_auxiliary_central_value_gev"] == pytest.approx(
        expected_delta,
        abs=1.0e-12,
    )
    object_ids = {item["id"] for item in payload["constructive_objects"]}
    assert "cross_section_to_direct_top_response_kernel" in object_ids
    assert "direct_top_uncertainty_propagation" in object_ids
    assert "Q007TP_direct_top_central_value_as_calibration_input" in payload["forbidden_solver_inputs"]
