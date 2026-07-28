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


def test_rg_matching_threshold_contract_is_open_partial_and_nonpromoting() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_rg_matching_threshold_contract"
    assert payload["github_issue"] == 32
    assert payload["status"] == "open_source_rg_frontier_partial"
    assert payload["promotion_allowed"] is False
    assert payload["github_issue_state"] == "open"
    assert payload["github_dependencies"] == [569, 630, 631, 632, 634]
    assert payload["source_frontier"]["status"] == (
        "PARTIAL_EXACT_REPRESENTATION_INDICES__SOURCE_MATCHING_OPEN"
    )
    assert payload["worker_result_policy"]["partial_frontier_allowed"] is True
    assert payload["worker_result_policy"]["closure_requires_complete_source_packet"] is True
    assert payload["closure_gate"]["closable_now"] is False
    assert payload["closure_gate"]["closed_as"] is None
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
    assert "zero-vertex decoupling" in payload["closure_gate"]["reason"]
    assert "using_threshold_choices_as_hidden_fit_parameters" in payload["forbidden_promotions"]
    assert (
        "reusing_the_issue_593_external_validation_packet_as_an_OPH_source"
        in payload["forbidden_promotions"]
    )
