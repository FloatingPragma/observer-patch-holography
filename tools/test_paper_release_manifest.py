#!/usr/bin/env python3
"""Regression tests for tools/validate_paper_release_manifest.py (issue #514).

Run with: python -m pytest tools/test_paper_release_manifest.py
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

import generate_paper_release_manifest as generator
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
    assert set(sections["papers"]) == set(validator.source.RELEASE_TRACKED)
    assert set(sections["extra_papers"]) == set(validator.source.EXTRA_PAPERS)
    assert set(sections["supplemental_papers"]) == set(validator.source.PAPERS) - set(validator.source.RELEASE_TRACKED)


def test_manifest_generation_timestamp_is_release_derived() -> None:
    assert (
        generator.deterministic_generated_at("July 25, 2026")
        == "2026-07-25T00:00:00Z"
    )


def test_generator_rejects_changed_pdf_behind_existing_tag() -> None:
    tagged = _base()
    current = copy.deepcopy(tagged)
    paper_id = next(iter(current["papers"]))
    current["papers"][paper_id]["sha256"] = "0" * 64

    with pytest.raises(SystemExit, match="immutable Git tag.*--allow-same-release"):
        generator.enforce_tag_immutability(tagged, current)


def test_generator_accepts_unchanged_pdf_behind_existing_tag() -> None:
    tagged = _base()
    generator.enforce_tag_immutability(tagged, copy.deepcopy(tagged))


def test_generator_accepts_new_untagged_release() -> None:
    current = _base()
    current["release_id"] = "r-next"
    generator.enforce_tag_immutability(_base(), current)


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


def test_rejects_cross_paper_artifact_mapping(tmp_path: Path) -> None:
    """A valid digest for the wrong paper cannot satisfy source-derived membership."""
    manifest = _base()
    left, right = list(manifest["papers"])[:2]
    manifest["papers"][left] = copy.deepcopy(manifest["papers"][right])
    problems = validator.validate(_write(tmp_path, manifest))
    assert any(left in p and "source derives" in p for p in problems)


def test_rejects_sha256_mismatch(tmp_path: Path) -> None:
    # A listed PDF whose content differs from its declared digest (silent rebuild /
    # swap / tamper) must be rejected even when the path exists.
    manifest = _base()
    paper_id = next(iter(manifest["papers"]))
    manifest["papers"][paper_id]["sha256"] = "0" * 64
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("sha256 mismatch" in p and paper_id in p for p in problems)


def test_rejects_size_mismatch(tmp_path: Path) -> None:
    # A truncated / regrown artifact whose byte count differs from the manifest is rejected.
    manifest = _base()
    paper_id = next(iter(manifest["papers"]))
    manifest["papers"][paper_id]["size_bytes"] = "1"
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("size_bytes mismatch" in p and paper_id in p for p in problems)


def test_generator_rejects_pdf_not_implied_by_source(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Exercise the actual release-surface stray-artifact boundary."""
    paper_dir = tmp_path / "paper"
    extra_dir = tmp_path / "extra"
    paper_dir.mkdir()
    extra_dir.mkdir()
    expected_pdf = paper_dir / "expected.pdf"
    expected_pdf.write_bytes(b"%PDF expected")
    manifest = {
        "papers": {
            "expected": {
                "pdf_path": "paper/expected.pdf",
                "sha256": "unused-by-this-check",
                "size_bytes": expected_pdf.stat().st_size,
            }
        },
        "supplemental_papers": {},
        "extra_papers": {},
    }
    monkeypatch.setattr(generator, "NON_TEX_SOURCE_PDFS", {})

    generator.verify_no_stray_pdfs(tmp_path, manifest)
    (extra_dir / "not_implied_by_any_source.pdf").write_bytes(b"%PDF stray")

    with pytest.raises(SystemExit, match="stray PDFs.*not_implied_by_any_source"):
        generator.verify_no_stray_pdfs(tmp_path, manifest)


def test_generator_rejects_missing_registered_non_tex_source(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "paper").mkdir()
    (tmp_path / "extra").mkdir()
    monkeypatch.setattr(
        generator,
        "NON_TEX_SOURCE_PDFS",
        {
            Path("extra/derived.pdf"): Path("extra/derived/build.sh"),
        },
    )

    with pytest.raises(SystemExit, match="registered source.*missing"):
        generator.verify_no_stray_pdfs(
            tmp_path,
            {
                "papers": {},
                "supplemental_papers": {},
                "extra_papers": {},
            },
        )


def test_generator_rejects_missing_registered_non_tex_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "paper").mkdir()
    source = tmp_path / "extra" / "derived" / "build.sh"
    source.parent.mkdir(parents=True)
    source.write_text("#!/bin/sh\n", encoding="utf-8")
    monkeypatch.setattr(
        generator,
        "NON_TEX_SOURCE_PDFS",
        {
            Path("extra/derived.pdf"): Path("extra/derived/build.sh"),
        },
    )

    with pytest.raises(SystemExit, match="registered non-TeX output is missing"):
        generator.verify_no_stray_pdfs(
            tmp_path,
            {
                "papers": {},
                "supplemental_papers": {},
                "extra_papers": {},
            },
        )


def test_validator_rejects_missing_registered_non_tex_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "paper").mkdir()
    source_path = tmp_path / "extra" / "derived" / "build.sh"
    source_path.parent.mkdir(parents=True)
    source_path.write_text("#!/bin/sh\n", encoding="utf-8")
    monkeypatch.setattr(validator, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(
        validator.source,
        "NON_TEX_SOURCE_PDFS",
        {
            Path("extra/derived.pdf"): Path("extra/derived/build.sh"),
        },
    )
    problems: list[str] = []
    validator.check_release_surface(
        {
            "papers": {},
            "supplemental_papers": {},
            "extra_papers": {},
        },
        problems,
    )
    assert any("registered non-TeX output is missing" in p for p in problems)


def test_publication_ci_rejects_committed_artifact_drift_after_rebuild() -> None:
    workflow = (
        REPO_ROOT / ".github" / "workflows" / "publication-build.yml"
    ).read_text(encoding="utf-8")
    assert "git diff --exit-code --" in workflow
    assert "book/reverse-engineering-reality-book.pdf" in workflow
    assert "python tools/refresh_paper_release.py" in workflow
    assert "python tools/build_book_pdf.py" in workflow
