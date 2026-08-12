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
