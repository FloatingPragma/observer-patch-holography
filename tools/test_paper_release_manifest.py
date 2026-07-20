#!/usr/bin/env python3
"""Regression tests for tools/validate_paper_release_manifest.py (issue #514).

Run with: python -m pytest tools/test_paper_release_manifest.py
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import validate_paper_release_manifest as validator

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = REPO_ROOT / "paper" / "paper_release_manifest.json"


def _base() -> dict:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def _write(tmp_path: Path, manifest: dict) -> Path:
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return path


def test_committed_manifest_matches_source_set(tmp_path: Path) -> None:
    # The manifest membership is derived from build_tex_papers.py, so the
    # committed manifest must validate with no problems.
    assert validator.validate(MANIFEST) == []


def test_expected_sets_are_derived_not_fixed() -> None:
    sections = validator.expected_sections()
    assert sections["papers"] == set(validator.source.RELEASE_TRACKED)
    assert sections["extra_papers"] == set(validator.source.EXTRA_PAPERS)
    assert sections["supplemental_papers"] == set(validator.source.PAPERS) - set(validator.source.RELEASE_TRACKED)


def test_rejects_missing_paper(tmp_path: Path) -> None:
    manifest = _base()
    removed = next(iter(manifest["papers"]))
    manifest["papers"].pop(removed)
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("missing" in p and removed in p for p in problems)


def test_rejects_unexpected_paper(tmp_path: Path) -> None:
    manifest = _base()
    manifest["papers"]["not_a_source_paper"] = {"pdf_path": "paper/x.pdf", "sha256": "x", "size_bytes": 1}
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("unexpected" in p and "not_a_source_paper" in p for p in problems)


def test_rejects_absent_artifact(tmp_path: Path) -> None:
    manifest = _base()
    paper_id = next(iter(manifest["papers"]))
    manifest["papers"][paper_id] = {"pdf_path": "paper/DOES_NOT_EXIST.pdf", "sha256": "x", "size_bytes": 1}
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("missing on disk" in p for p in problems)


def test_rejects_sha256_mismatch(tmp_path: Path) -> None:
    # A listed PDF whose content no longer matches its declared digest (silent rebuild /
    # swap / tamper) must be rejected, not silently accepted because the path still exists.
    manifest = _base()
    paper_id = next(iter(manifest["papers"]))
    manifest["papers"][paper_id]["sha256"] = "0" * 64
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("sha256 mismatch" in p and paper_id in p for p in problems)


def test_rejects_size_mismatch(tmp_path: Path) -> None:
    # A truncated / regrown artifact whose byte count no longer matches the manifest is rejected.
    manifest = _base()
    paper_id = next(iter(manifest["papers"]))
    manifest["papers"][paper_id]["size_bytes"] = "1"
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("size_bytes mismatch" in p and paper_id in p for p in problems)
