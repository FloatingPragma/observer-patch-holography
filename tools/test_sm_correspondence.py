"""Drift and tamper tests for the SM Lagrangian correspondence surface."""

from __future__ import annotations

import json

import pytest

import build_sm_correspondence as tool


def _data() -> dict:
    return json.loads(tool.TABLE_PATH.read_text(encoding="utf-8"))


def test_committed_table_validates() -> None:
    rows = tool.validate(_data())
    assert [row["id"] for row in rows] == [
        f"SML-{index:02d}" for index in range(1, len(rows) + 1)
    ]


def test_rendered_surface_matches_committed_doc() -> None:
    data = _data()
    rows = tool.validate(data)
    assert tool.SURFACE_PATH.read_bytes() == tool.render(data, rows).encode("utf-8")


def test_check_mode_passes_on_committed_artifacts(capsys) -> None:
    assert tool.main(["--check"]) == 0
    assert "surface is current" in capsys.readouterr().out


def test_classification_outside_enum_is_rejected() -> None:
    data = _data()
    data["rows"][0]["classification"] = "derived"
    with pytest.raises(SystemExit, match="classification"):
        tool.validate(data)


def test_upgrading_an_absent_row_without_premises_is_rejected() -> None:
    data = _data()
    row = next(r for r in data["rows"] if r["classification"] == "absent")
    row["classification"] = "derived_conditional"
    with pytest.raises(SystemExit, match="premise"):
        tool.validate(data)


def test_absent_row_naming_premises_is_rejected() -> None:
    data = _data()
    row = next(r for r in data["rows"] if r["classification"] == "absent")
    row["premises"] = ["PR-11"]
    with pytest.raises(SystemExit, match="absent"):
        tool.validate(data)


def test_classification_tamper_breaks_surface_parity() -> None:
    data = _data()
    row = next(r for r in data["rows"] if r["classification"] == "partial")
    row["classification"] = "registered_premise"
    rows = tool.validate(data)
    assert tool.render(data, rows).encode("utf-8") != tool.SURFACE_PATH.read_bytes()


def test_missing_evidence_path_is_rejected() -> None:
    data = _data()
    data["rows"][0]["evidence"] = ["docs/DOES_NOT_EXIST_SML.md"]
    with pytest.raises(SystemExit, match="evidence path missing"):
        tool.validate(data)


def test_unknown_premise_id_is_rejected() -> None:
    data = _data()
    data["rows"][0]["premises"] = ["PR-99"]
    with pytest.raises(SystemExit, match="unknown premise"):
        tool.validate(data)


def test_premise_register_source_tamper_is_rejected() -> None:
    data = _data()
    data["premise_register_source"] = "tracking/another_register.json"
    with pytest.raises(SystemExit, match="premise_register_source"):
        tool.validate(data)


def test_dropped_term_is_rejected() -> None:
    data = _data()
    data["rows"].pop()
    with pytest.raises(SystemExit, match="exactly 19 entries"):
        tool.validate(data)


def test_renamed_term_is_rejected() -> None:
    data = _data()
    data["rows"][0]["term"] = "Color term: renamed"
    with pytest.raises(SystemExit, match="canonical inventory entry"):
        tool.validate(data)


def test_unknown_open_premise_id_is_rejected() -> None:
    data = _data()
    data["rows"][0]["open_premises"] = ["PR-99"]
    with pytest.raises(SystemExit, match="unknown open premise"):
        tool.validate(data)


def test_banned_wording_is_rejected() -> None:
    data = _data()
    banned = "hon" + "est"
    data["rows"][0]["boundary"] = f"This shape is {banned} and derived."
    with pytest.raises(SystemExit, match="banned word"):
        tool.validate(data)


def test_duplicate_json_key_is_rejected(tmp_path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key"):
        tool.load_json(duplicate)
