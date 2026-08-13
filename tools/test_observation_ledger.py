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


def test_unknown_architecture_version_fails() -> None:
    ledger = _ledger()
    ledger["rows"][0]["architecture_version"] = "AV-999"
    with pytest.raises(SystemExit, match="architecture_version"):
        tool.validate(ledger)


def test_dropped_observation_row_fails() -> None:
    ledger = _ledger()
    ledger["rows"] = [row for row in ledger["rows"] if row["id"] != "OL-E2"]
    with pytest.raises(SystemExit, match="fixed ordered observation inventory"):
        tool.validate(ledger)


def test_new_reverse_consumer_gap_fails() -> None:
    ledger = _ledger()
    for row in ledger["rows"]:
        if row["id"] in {"OL-K2", "OL-K3"}:
            row["open_premises"] = []
    with pytest.raises(SystemExit, match="reverse-map exceptions drifted"):
        tool.validate(ledger)


def test_duplicate_json_key_fails_at_load(tmp_path) -> None:
    path = tmp_path / "duplicate.json"
    path.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key"):
        tool.load_ledger(path)
