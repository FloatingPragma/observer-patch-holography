"""Drift and tamper rejection tests for the V3 observation ledger surface."""

from __future__ import annotations

import json

import pytest

import build_observation_ledger as tool


def _ledger() -> dict:
    return json.loads(tool.LEDGER_PATH.read_text(encoding="utf-8"))


def test_committed_surface_matches_rebuild() -> None:
    rows = tool.validate(_ledger())
    assert tool.SURFACE_PATH.read_bytes() == tool.render(rows).encode("utf-8")
    assert tool.main(["--check"]) == 0


def test_ol_a1_failed_control_and_unfrozen_followup_are_exact() -> None:
    rows = {row["id"]: row for row in _ledger()["rows"]}
    row = rows["OL-A1"]
    assert row["status"] == "owed"
    assert row["rung"] == "emergent"
    assert "claims/emergent_instrument_register.json" in row["evidence"]
    assert "docs/OL_A1_FACTORIAL_FOLLOWUP_DESIGN.md" not in row["evidence"]
    assert "threshold form on two of five" in row["notes"]
    assert "robust reference on one of five" in row["notes"]
    assert "does not prove carrier count is the cause" in row["notes"]
    assert "carrier-count-by-absolute-support-size" in row["notes"]
    assert "not frozen or authorized to run" in row["notes"]


def test_governance_row_and_lane_738_are_absent() -> None:
    rows = _ledger()["rows"]
    assert "OL-I1" not in {row["id"] for row in rows}
    assert 738 not in {row["lane_issue"] for row in rows}
    i2 = next(row for row in rows if row["id"] == "OL-I2")
    assert "owning physics lane" in i2["notes"]


def test_predictive_rows_name_only_registered_frozen_targets() -> None:
    rows = [row for row in _ledger()["rows"] if row["rung"] == "predictive"]
    assert rows
    assert all("frozen_targets" in row for row in rows)
    assert all(len(row["frozen_targets"]) == len(set(row["frozen_targets"])) for row in rows)
    known = set(tool.load_frozen_target_rows())
    assert all(target in known for row in rows for target in row["frozen_targets"])


def test_unknown_frozen_target_fails() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["rung"] == "predictive")
    row["frozen_targets"] = ["FZ-99"]
    with pytest.raises(SystemExit, match="is not on the register"):
        tool.validate(ledger)


def test_duplicate_frozen_target_fails() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["id"] == "OL-F4")
    row["frozen_targets"] = ["FZ-11", "FZ-11"]
    with pytest.raises(SystemExit, match="duplicate-free"):
        tool.validate(ledger)


def test_attained_predictive_row_requires_locked_target() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["id"] == "OL-B4")
    row["status"] = "attained"
    row["open_premises"] = []
    assert row["frozen_targets"] == ["FZ-06"]
    with pytest.raises(SystemExit, match="all targets in its fixed contract"):
        tool.validate(ledger)


def test_attained_predictive_row_accepts_locked_target(monkeypatch) -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["id"] == "OL-F4")
    row["status"] = "attained"
    row["open_premises"] = []
    premise_rows, _ = tool.load_premise_register()
    premise_rows = json.loads(json.dumps(premise_rows))
    pr15 = next(item for item in premise_rows if item["id"] == "PR-15")
    pr15["consuming_lanes"].remove(733)
    monkeypatch.setattr(
        tool,
        "load_premise_register",
        lambda: (premise_rows, {item["id"]: item for item in premise_rows}),
    )
    tool.validate(ledger)


def test_unrelated_frozen_target_cannot_false_green_predictive_row() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["id"] == "OL-H7")
    row["frozen_targets"] = ["FZ-11"]
    row["status"] = "attained"
    row["open_premises"] = []
    with pytest.raises(SystemExit, match="fixed scientific contract"):
        tool.validate(ledger)


def test_partial_predictive_target_set_cannot_false_green_row() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["id"] == "OL-F4")
    row["frozen_targets"] = ["FZ-11"]
    row["status"] = "attained"
    row["open_premises"] = []
    with pytest.raises(SystemExit, match="fixed scientific contract"):
        tool.validate(ledger)


def test_attained_predictive_row_requires_every_target_locked() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["id"] == "OL-F4")
    row["status"] = "attained"
    row["open_premises"] = []
    frozen = tool.load_frozen_target_rows()
    frozen["FZ-12"]["status"] = "registered_pending_freeze"
    original = tool.load_frozen_target_rows
    try:
        tool.load_frozen_target_rows = lambda: frozen
        with pytest.raises(SystemExit, match="all targets in its fixed contract"):
            tool.validate(ledger)
    finally:
        tool.load_frozen_target_rows = original


def test_attained_row_cannot_retain_open_premise() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["id"] == "OL-C1")
    row["status"] = "attained"
    row["open_premises"] = ["PR-04"]
    with pytest.raises(SystemExit, match="cannot retain open premises"):
        tool.validate(ledger)


def test_nonpredictive_row_cannot_carry_frozen_targets() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["rung"] != "predictive")
    row["frozen_targets"] = []
    with pytest.raises(SystemExit, match="keys mismatch"):
        tool.validate(ledger)


def test_check_catches_mutated_status(tmp_path, monkeypatch) -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["status"] == "attained")
    row["status"] = "partial"
    mutated = tmp_path / "observation_ledger.json"
    mutated.write_text(json.dumps(ledger), encoding="utf-8")
    monkeypatch.setattr(tool, "LEDGER_PATH", mutated)
    assert tool.main(["--check"]) == 1


def test_missing_evidence_path_fails() -> None:
    ledger = _ledger()
    ledger["rows"][0]["evidence"] = ["docs/NO_SUCH_RECEIPT.md"]
    with pytest.raises(SystemExit, match="evidence path missing"):
        tool.validate(ledger)


def test_invalid_rung_fails() -> None:
    ledger = _ledger()
    ledger["rows"][0]["rung"] = "cosmic"
    with pytest.raises(SystemExit, match="rung"):
        tool.validate(ledger)


def test_owed_row_may_name_registered_open_premises() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["id"] == "OL-A4")
    assert row["status"] == "owed"
    assert row["premises"] == []
    assert row["open_premises"] == ["PR-16", "PR-52"]
    tool.validate(ledger)


def test_reverse_premise_lane_mismatch_fails() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["id"] == "OL-A2")
    row["open_premises"] = ["PR-43"]
    with pytest.raises(SystemExit, match="does not declare consuming lane #728"):
        tool.validate(ledger)


def test_dropped_observation_row_fails() -> None:
    ledger = _ledger()
    ledger["rows"] = [row for row in ledger["rows"] if row["id"] != "OL-E2"]
    with pytest.raises(SystemExit, match="fixed ordered observation inventory"):
        tool.validate(ledger)


def test_new_reverse_consumer_gap_fails() -> None:
    ledger = _ledger()
    for row in ledger["rows"]:
        if row["id"] in {"OL-K2", "OL-K3", "OL-K5", "OL-K6"}:
            row["open_premises"] = []
    with pytest.raises(SystemExit, match="reverse-map exceptions drifted"):
        tool.validate(ledger)


def test_ol_e1_cannot_hide_the_calibration_premise_in_another_lane_row() -> None:
    ledger = _ledger()
    row = next(row for row in ledger["rows"] if row["id"] == "OL-E1")
    row["premises"].remove("PR-15")
    with pytest.raises(SystemExit, match="fixed row-level contract"):
        tool.validate(ledger)


def test_duplicate_json_key_fails_at_load(tmp_path) -> None:
    path = tmp_path / "duplicate.json"
    path.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key"):
        tool.load_ledger(path)
