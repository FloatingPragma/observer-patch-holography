"""Drift and tamper rejection tests for the V3 premise register surface."""

from __future__ import annotations

import copy
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_premise_register as register_tool


def _register() -> dict:
    return copy.deepcopy(register_tool.load_json(register_tool.REGISTER_PATH))


def test_committed_register_validates() -> None:
    rows = register_tool.validate(_register())
    assert [row["id"] for row in rows] == [
        expected[0] for expected in register_tool.EXPECTED_ROWS
    ]


def test_rebuild_parity_with_committed_surface() -> None:
    rows = register_tool.validate(_register())
    rendered = register_tool.render(rows).encode("utf-8")
    committed = register_tool.SURFACE_PATH.read_bytes()
    assert rendered == committed


def test_discharge_link_names_durable_register_custody() -> None:
    rendered = register_tool.render(register_tool.validate(_register()))
    assert "durable architecture/audit-register custody" in rendered
    assert "issue #741 decision/version custody" not in rendered


def test_check_mode_passes_on_committed_artifacts() -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).parent / "build_premise_register.py"),
            "--check",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_disposition_tamper_rejected() -> None:
    register = _register()
    assert register["rows"][2]["id"] == "PR-03"
    register["rows"][2]["disposition"] = "axiomatize"
    with pytest.raises(SystemExit, match="disposition must equal"):
        register_tool.validate(register)


def test_unknown_type_rejected() -> None:
    register = _register()
    register["rows"][0]["type"] = "aesthetic_preference"
    with pytest.raises(SystemExit, match="not in the enum"):
        register_tool.validate(register)


def test_unknown_disposition_rejected() -> None:
    register = _register()
    register["rows"][0]["disposition"] = "postpone"
    with pytest.raises(SystemExit, match="not in the enum"):
        register_tool.validate(register)


def test_missing_evidence_path_rejected() -> None:
    register = _register()
    register["rows"][0]["evidence"].append("docs/NO_SUCH_EVIDENCE_FILE.md")
    with pytest.raises(SystemExit, match="evidence path missing"):
        register_tool.validate(register)


def test_lane_out_of_range_rejected() -> None:
    register = _register()
    register["rows"][0]["consuming_lanes"] = [727]
    with pytest.raises(SystemExit, match="outside"):
        register_tool.validate(register)


def test_renamed_row_rejected() -> None:
    register = _register()
    register["rows"][4]["name"] = "path reference"
    with pytest.raises(SystemExit, match="name must equal"):
        register_tool.validate(register)


def test_dropped_row_rejected() -> None:
    register = _register()
    register["rows"].pop()
    with pytest.raises(SystemExit, match="exactly 57 entries"):
        register_tool.validate(register)


def test_stale_surface_fails_check(tmp_path, monkeypatch) -> None:
    stale = tmp_path / "PREMISE_REGISTER_V3.md"
    stale.write_bytes(b"stale\n")
    monkeypatch.setattr(register_tool, "SURFACE_PATH", stale)
    monkeypatch.setattr(sys, "argv", ["build_premise_register.py", "--check"])
    assert register_tool.main() == 1


def test_duplicate_json_key_rejected(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key"):
        register_tool.load_json(duplicate)
