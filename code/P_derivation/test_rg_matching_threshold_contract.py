#!/usr/bin/env python3
"""Validate the RG matching and threshold constructive contract."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parent
SCRIPT = ROOT / "rg_matching_threshold_contract.py"
OUTPUT = ROOT / "runtime" / "rg_matching_threshold_contract_current.json"


def test_rg_matching_threshold_contract_is_partial_and_nonpromoting() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_rg_matching_threshold_contract"
    assert payload["status"] == "open_source_rg_frontier_partial"
    assert payload["promotion_allowed"] is False
    assert set(payload) == {
        "artifact",
        "constructive_objects",
        "forbidden_promotions",
        "promotion_allowed",
        "promotion_boundary",
        "scientific_boundary",
        "source_frontier",
        "status",
    }
    assert payload["source_frontier"]["status"] == (
        "PARTIAL_EXACT_REPRESENTATION_INDICES__SOURCE_MATCHING_OPEN"
    )
    assert payload["promotion_boundary"]["promotion_allowed"] is False
    object_ids = {item["id"] for item in payload["constructive_objects"]}
    assert object_ids == {
        "representation_index_frontier",
        "scheme_lock",
        "threshold_map",
        "beta_provenance_table",
        "matching_interval_composition_certificate",
    }
    statuses = {
        item["id"]: item["current_status"] for item in payload["constructive_objects"]
    }
    assert statuses["representation_index_frontier"] == (
        "complete_at_finite_representation_scope"
    )
    assert statuses["beta_provenance_table"] == (
        "partial_parametric_gauge_one_loop_only"
    )
    assert statuses["threshold_map"] == "not_emitted"
    assert "zero-vertex decoupling" in payload["promotion_boundary"]["reason"]
    assert "representation indices alone do not determine" in payload[
        "scientific_boundary"
    ]["corpus_limited_no_go"]
    assert "using_threshold_choices_as_hidden_fit_parameters" in payload["forbidden_promotions"]
    assert (
        "reusing_an_external_validation_packet_as_an_OPH_source"
        in payload["forbidden_promotions"]
    )
