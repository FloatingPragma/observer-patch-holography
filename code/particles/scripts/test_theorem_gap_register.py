"""Tests for the theorem gap register."""

from __future__ import annotations

import build_theorem_gap_register as lane


def test_register_counts_are_consistent():
    register = lane.build()
    counts = register["counts"]
    assert counts["closed"] == len(register["closed_this_program"])
    assert sum(
        v for k, v in counts.items() if k != "closed"
    ) <= len(register["open_register"])


def test_closures_carry_artifacts():
    register = lane.build()
    assert all("artifact" in e for e in register["closed_this_program"])
    ids = {e["id"] for e in register["closed_this_program"]}
    assert "FLOW_INTERNAL_SELECTION_NO_GO" in ids
    assert "ANCHOR_RECONCILIATION_MIDPOINT_IMPLICATION" in ids
    assert "PHOTON_EXACT_MASSLESSNESS_RECEIPT" in ids


def test_open_register_names_the_faithfulness_gate_and_the_test():
    register = lane.build()
    ids = {e["id"] for e in register["open_register"]}
    assert "CARRIER_MODEL_FAITHFULNESS" in ids
    assert "CRITICALITY_BOUNDARY_SCALE_SELECTION" in ids
    closed_ids = {e["id"] for e in register["closed_this_program"]}
    assert "AR_PREMISE_REDUCTION" in closed_ids
    assert "CF1_CF2_MODEL_LEVEL_CENSUS" in closed_ids
    assert register["promotion_allowed"] is False


def test_markdown_renders_both_sections():
    register = lane.build()
    text = lane.render_markdown(register)
    assert "## Closed" in text
    assert "## Open" in text
