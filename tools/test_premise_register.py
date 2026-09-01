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


def test_intro_keeps_register_scientific_and_explains_evidence_roles() -> None:
    rendered = register_tool.render(register_tool.validate(_register()))
    assert "Established under [issue #727]" in rendered
    assert "maintained as a scientific register" in rendered
    assert "implicit reverse-consumer edge" in rendered
    assert "premise-discharge queue" not in rendered


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
    with pytest.raises(SystemExit, match="not a current scientific lane"):
        register_tool.validate(register)


def test_retired_lane_738_rejected() -> None:
    register = _register()
    register["rows"][0]["consuming_lanes"] = [738]
    with pytest.raises(SystemExit, match="not a current scientific lane"):
        register_tool.validate(register)


def test_evidence_role_path_parity_is_exact() -> None:
    register = _register()
    row = register["rows"][0]
    row["evidence_roles"].pop(row["evidence"][0])
    with pytest.raises(SystemExit, match="exact evidence-path parity"):
        register_tool.validate(register)


def test_unknown_evidence_role_rejected() -> None:
    register = _register()
    row = register["rows"][0]
    row["evidence_roles"][row["evidence"][0]] = "related"
    with pytest.raises(SystemExit, match="must be one of"):
        register_tool.validate(register)


def test_pr53_has_no_retired_custody_lane() -> None:
    register = _register()
    row = next(row for row in register["rows"] if row["id"] == "PR-53")
    assert row["consuming_lanes"] == [733, 736, 740]
    assert 738 not in row["consuming_lanes"]


def test_renamed_row_rejected() -> None:
    register = _register()
    register["rows"][4]["name"] = "path reference"
    with pytest.raises(SystemExit, match="name must equal"):
        register_tool.validate(register)


def test_dropped_row_rejected() -> None:
    register = _register()
    register["rows"].pop()
    with pytest.raises(SystemExit, match="exactly 82 entries"):
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
