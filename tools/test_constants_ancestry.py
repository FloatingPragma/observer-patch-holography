"""Drift and tamper rejection tests for the V3 constants ancestry surface."""

from __future__ import annotations

import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_constants_ancestry as ancestry_tool


def _register() -> dict:
    return copy.deepcopy(ancestry_tool.load_json(ancestry_tool.REGISTER_PATH))


def _row(register: dict, row_id: str) -> dict:
    for row in register["rows"]:
        if row["id"] == row_id:
            return row
    raise AssertionError(f"row {row_id} absent")


def test_committed_register_validates() -> None:
    rows = ancestry_tool.validate(_register())
    assert [row["id"] for row in rows] == [
        expected[0] for expected in ancestry_tool.EXPECTED_ROWS
    ]


def test_rebuild_parity_with_committed_surface() -> None:
    rows = ancestry_tool.validate(_register())
    rendered = ancestry_tool.render(rows).encode("utf-8")
    committed = ancestry_tool.SURFACE_PATH.read_bytes()
    assert rendered == committed


def test_check_mode_passes_on_committed_artifacts() -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parent / "build_constants_ancestry.py"),
            "--check",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_status_tamper_rejected() -> None:
    register = _register()
    _row(register, "CA-01")["status"] = "conditional_determination"
    with pytest.raises(SystemExit, match="status must equal"):
        ancestry_tool.validate(register)


def test_unknown_status_rejected() -> None:
    register = _register()
    _row(register, "CA-01")["status"] = "prediction"
    with pytest.raises(SystemExit, match="not in the enum"):
        ancestry_tool.validate(register)


def test_unknown_class_rejected() -> None:
    register = _register()
    _row(register, "CA-01")["ancestry"][0]["class"] = "vibes"
    with pytest.raises(SystemExit, match="not in the enum"):
        ancestry_tool.validate(register)


def test_missing_evidence_path_rejected() -> None:
    register = _register()
    _row(register, "CA-01")["evidence"].append("docs/NO_SUCH_EVIDENCE_FILE.md")
    with pytest.raises(SystemExit, match="evidence path missing"):
        ancestry_tool.validate(register)


def test_missing_falsifier_artifact_rejected() -> None:
    register = _register()
    _row(register, "CA-02")["falsifier"]["artifact"] = "docs/NO_SUCH_RULE.md"
    with pytest.raises(SystemExit, match="artifact path missing"):
        ancestry_tool.validate(register)


def test_unknown_register_row_rejected() -> None:
    register = _register()
    _row(register, "CA-06")["ancestry"][0]["register_row"] = "PR-99"
    with pytest.raises(SystemExit, match="neither an existing"):
        ancestry_tool.validate(register)


def test_proposed_entry_citing_pr_id_rejected() -> None:
    register = _register()
    row = _row(register, "CA-01")
    for entry in row["ancestry"]:
        if entry.get("register_row", "").startswith("PR-3"):
            entry["register_row"] = "proposed"
            entry["input"] = "declared premise cited as PR-20"
            break
    else:
        raise AssertionError("no register_import entry available to repurpose")
    with pytest.raises(SystemExit, match="cites no PR id"):
        ancestry_tool.validate(register)


def test_pn_entry_without_pr17_rejected() -> None:
    register = _register()
    _row(register, "CA-01")["ancestry"][0]["register_row"] = "PR-14"
    with pytest.raises(SystemExit, match="must cite register row PR-17"):
        ancestry_tool.validate(register)


def test_owed_row_with_ancestry_rejected() -> None:
    register = _register()
    _row(register, "CA-07")["ancestry"].append(
        {"input": "a smuggled input", "class": "derived"}
    )
    with pytest.raises(SystemExit, match="consumes nothing"):
        ancestry_tool.validate(register)


def test_postdiction_without_measured_target_rejected() -> None:
    register = _register()
    row = _row(register, "CA-04")
    row["ancestry"] = [
        entry
        for entry in row["ancestry"]
        if entry["class"] != "measured_comparison_target"
    ]
    with pytest.raises(SystemExit, match="plain sight"):
        ancestry_tool.validate(register)


def test_measured_target_on_conditional_determination_rejected() -> None:
    register = _register()
    _row(register, "CA-08")["ancestry"].append(
        {
            "input": "a measured witness on a determination row",
            "class": "measured_comparison_target",
        }
    )
    with pytest.raises(SystemExit, match="postdiction ledger"):
        ancestry_tool.validate(register)


def test_ledger_premise_mismatch_rejected() -> None:
    register = _register()
    row = _row(register, "CA-03")
    row["ancestry"] = [
        entry for entry in row["ancestry"] if entry.get("register_row") != "PR-13"
    ]
    with pytest.raises(SystemExit, match="observation-ledger premises"):
        ancestry_tool.validate(register)


def test_diagnostic_import_must_be_empirical_import() -> None:
    register = _register()
    _row(register, "CA-06")["diagnostic_import_register_rows"] = ["PR-11"]
    with pytest.raises(SystemExit, match="registered empirical_import"):
        ancestry_tool.validate(register)


def test_diagnostic_import_forbidden_on_postdiction_row() -> None:
    register = _register()
    _row(register, "CA-03")["diagnostic_import_register_rows"] = ["PR-14"]
    with pytest.raises(SystemExit, match="only on a diagnostic row"):
        ancestry_tool.validate(register)


def test_duplicate_json_key_rejected(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"first","schema":"second"}\n', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key 'schema'"):
        ancestry_tool.load_json(duplicate)


def test_nested_duplicate_json_key_rejected(tmp_path: Path) -> None:
    duplicate = tmp_path / "nested-duplicate.json"
    duplicate.write_text(
        json.dumps({"row": {"id": "CA-01"}}).replace(
            '"id": "CA-01"', '"id": "CA-01", "id": "CA-02"'
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(SystemExit, match="duplicate JSON key 'id'"):
        ancestry_tool.load_json(duplicate)


def test_dropped_row_rejected() -> None:
    register = _register()
    register["rows"].pop()
    with pytest.raises(SystemExit, match="exactly 8 entries"):
        ancestry_tool.validate(register)


def test_stale_surface_fails_check(tmp_path, monkeypatch) -> None:
    stale = tmp_path / "CONSTANTS_ANCESTRY_V3.md"
    stale.write_bytes(b"stale\n")
    monkeypatch.setattr(ancestry_tool, "SURFACE_PATH", stale)
    monkeypatch.setattr(sys, "argv", ["build_constants_ancestry.py", "--check"])
    assert ancestry_tool.main() == 1
