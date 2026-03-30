#!/usr/bin/env python3
"""Guard the neutrino absolute-attachment scaffold."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_absolute_attachment_scaffold.json"


def test_neutrino_absolute_attachment_scaffold_contract() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_neutrino_absolute_attachment_scaffold"
    assert payload["status"] == "minimal_constructive_extension"
    assert payload["exact_missing_object"] == "neutrino_weighted_cycle_absolute_attachment"
    assert payload["equivalent_scalar"]["name"] == "lambda_nu"
    assert payload["current_no_go"]["current_candidate_interface_artifact"] == "oph_majorana_overlap_defect_scalar_evaluator"
    assert payload["current_no_go"]["exact_next_theorem_object"] == "oph_same_label_overlap_defect_weight_normalizer"
    assert payload["current_no_go"]["strictly_smaller_missing_clause"] == "selector_overlap_phase_coboundary_trivializes_same_label_edge_transport"
    assert payload["extension_contract"]["must_emit"].startswith("lambda_nu")
    stack = payload["extension_contract"]["current_theorem_stack"]
    assert stack[0]["id"] == "selector_overlap_phase_coboundary_trivializes_same_label_edge_transport"
    assert stack[1]["id"] == "oph_same_label_overlap_defect_weight_normalizer"
