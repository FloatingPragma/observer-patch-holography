#!/usr/bin/env python3
"""Regression tests for exact warning-allowlist generator membership."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
import gen_warning_allowlist as gen  # noqa: E402


def test_registry_covers_exact_all_papers_plus_book() -> None:
    registered = gen.registered_logs()
    paper_logs = {
        tex_path.with_suffix(".log").resolve()
        for tex_path in gen.build_tex_papers.ALL_PAPERS.values()
    }
    assert set(registered) == paper_logs | {gen.BOOK_LOG.resolve()}


def _fixture_registry(tmp_path: Path) -> dict[Path, str]:
    paths = {
        tmp_path / "paper-a.log": "paper:a",
        tmp_path / "paper-b.log": "paper:b",
        tmp_path / "book.log": "book:reader-facing",
    }
    for path in paths:
        path.write_text("fixture log\n", encoding="utf-8")
    return {path.resolve(): label for path, label in paths.items()}


def test_exact_registered_log_set_is_accepted(tmp_path: Path) -> None:
    expected = _fixture_registry(tmp_path)
    supplied = list(reversed(expected))
    assert gen.validate_exact_coverage(supplied, expected) == sorted(expected)


def test_missing_registered_log_is_rejected(tmp_path: Path) -> None:
    expected = _fixture_registry(tmp_path)
    omitted = next(path for path, label in expected.items() if label == "paper:b")
    with pytest.raises(ValueError, match="missing registered logs") as exc:
        gen.validate_exact_coverage(
            [path for path in expected if path != omitted], expected
        )
    assert str(omitted) in str(exc.value)


def test_duplicate_registered_log_is_rejected(tmp_path: Path) -> None:
    expected = _fixture_registry(tmp_path)
    duplicate = next(iter(expected))
    with pytest.raises(ValueError, match="duplicate logs") as exc:
        gen.validate_exact_coverage([*expected, duplicate], expected)
    assert str(duplicate) in str(exc.value)


def test_unregistered_log_is_rejected(tmp_path: Path) -> None:
    expected = _fixture_registry(tmp_path)
    extra = tmp_path / "old-unregistered.log"
    extra.write_text("stale log\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unregistered logs") as exc:
        gen.validate_exact_coverage([*expected, extra], expected)
    assert str(extra.resolve()) in str(exc.value)


def test_registered_but_unbuilt_log_is_rejected(tmp_path: Path) -> None:
    expected = _fixture_registry(tmp_path)
    absent = next(iter(expected))
    absent.unlink()
    with pytest.raises(ValueError, match="not found on disk"):
        gen.validate_exact_coverage(list(expected), expected)


def test_help_exits_cleanly() -> None:
    script = TOOLS / "gen_warning_allowlist.py"
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        cwd=gen.REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "usage:" in result.stdout.lower()
    assert "exactly match the registry" in result.stdout
