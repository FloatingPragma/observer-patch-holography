from __future__ import annotations

from datetime import datetime, timezone

import pytest

from reproducible_build_env import build_environment


def test_release_date_fixes_pdf_clock_and_timezone(tmp_path, monkeypatch) -> None:
    paper = tmp_path / "paper"
    paper.mkdir()
    (paper / "release_info.tex").write_text(
        r"\newcommand{\OPHPaperReleaseDate}{July 23, 2026}" + "\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "1")
    monkeypatch.setenv("TZ", "Europe/Vienna")

    env = build_environment(tmp_path)

    expected = datetime(2026, 7, 23, tzinfo=timezone.utc)
    assert env["SOURCE_DATE_EPOCH"] == str(int(expected.timestamp()))
    assert env["TZ"] == "UTC"


def test_missing_release_date_fails_closed(tmp_path) -> None:
    paper = tmp_path / "paper"
    paper.mkdir()
    (paper / "release_info.tex").write_text("% missing\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="OPHPaperReleaseDate"):
        build_environment(tmp_path)
