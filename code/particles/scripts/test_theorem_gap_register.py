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
    assert "MAXWELL_CLASSICAL_ZERO_HARD_MASS_MODE_RECEIPT" in ids
    assert "YANG_MILLS_CLASSICAL_ZERO_HARD_MASS_MODE_RECEIPT" in ids
    assert "EINSTEIN_CLASSICAL_ZERO_HARD_MASS_MODE_RECEIPT" in ids
    assert "PHOTON_EXACT_MASSLESSNESS_RECEIPT" not in ids
    carrier_rows = [
        entry
        for entry in register["closed_this_program"]
        if entry["family"] == "conditional classical carrier modes"
    ]
    assert len(carrier_rows) == 3
    assert all(
        entry["artifact"] == "runs/status/carrier_mode_acceptance.json"
        for entry in carrier_rows
    )
    carrier_by_id = {entry["id"]: entry for entry in carrier_rows}
    assert "quantum-particle" in carrier_by_id[
        "MAXWELL_CLASSICAL_ZERO_HARD_MASS_MODE_RECEIPT"
    ]["statement"]
    assert "asymptotic gluon" in carrier_by_id[
        "YANG_MILLS_CLASSICAL_ZERO_HARD_MASS_MODE_RECEIPT"
    ]["statement"]
    assert "graviton Hilbert space" in carrier_by_id[
        "EINSTEIN_CLASSICAL_ZERO_HARD_MASS_MODE_RECEIPT"
    ]["statement"]


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
