#!/usr/bin/env python3
"""Regression tests for the canonical book PDF asset boundary."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

import book_pdf_assets as assets


def fixture_root(tmp_path: Path) -> Path:
    root = tmp_path / "repo"
    (root / "assets" / "book_diagrams").mkdir(parents=True)
    sources = [
        Path("assets/book-cover.svg"),
        Path("assets/pixel-constant.svg"),
        Path("assets/OPH_Unification_Diagram.svg"),
        Path("assets/book_diagrams/example.svg"),
    ]
    records = []
    for index, source_relative in enumerate(sources):
        source = root / source_relative
        source.write_text(f"<svg><text>{index}</text></svg>\\n", encoding="utf-8")
        rendering_relative = Path(
            assets.rendering_relative_for(
                assets.PurePosixPath(source_relative.as_posix())
            )
        )
        rendering = root / rendering_relative
        rendering.parent.mkdir(parents=True, exist_ok=True)
        rendering.write_bytes(f"%PDF-1.5\\nasset {index}\\n".encode())
        records.append(
            {
                "source": source_relative.as_posix(),
                "source_sha256": assets.sha256_file(source),
                "rendering": rendering_relative.as_posix(),
                "rendering_sha256": assets.sha256_file(rendering),
            }
        )
    manifest = root / Path(assets.MANIFEST_RELATIVE)
    manifest.write_text(
        json.dumps({"schema_version": 1, "assets": records}),
        encoding="utf-8",
    )
    return root


def test_valid_inventory_stages_exact_rendering_bytes(tmp_path: Path) -> None:
    root = fixture_root(tmp_path)
    validated = assets.validate_book_pdf_assets(root)
    assert len(validated) == 4

    staged_root = tmp_path / "staged"
    mapping = assets.stage_book_pdf_assets(staged_root, root)
    assert mapping["../assets/book-cover.svg"] == staged_root / "book-cover.pdf"
    assert (staged_root / "book-cover.pdf").read_bytes() == (
        root / "assets/book_pdf_renderings/book-cover.pdf"
    ).read_bytes()


def test_changed_svg_requires_regeneration(tmp_path: Path) -> None:
    root = fixture_root(tmp_path)
    (root / "assets/book-cover.svg").write_text("<svg/>\\n", encoding="utf-8")
    with pytest.raises(assets.BookAssetValidationError, match="source SHA-256 mismatch"):
        assets.validate_book_pdf_assets(root)


def test_changed_rendering_is_rejected(tmp_path: Path) -> None:
    root = fixture_root(tmp_path)
    (root / "assets/book_pdf_renderings/book-cover.pdf").write_bytes(
        b"%PDF-1.5\\nsubstitute\\n"
    )
    with pytest.raises(
        assets.BookAssetValidationError,
        match="rendering SHA-256 mismatch",
    ):
        assets.validate_book_pdf_assets(root)


def test_missing_and_extra_sources_are_rejected(tmp_path: Path) -> None:
    root = fixture_root(tmp_path)
    (root / "assets/book_diagrams/new.svg").write_text("<svg/>\\n", encoding="utf-8")
    with pytest.raises(
        assets.BookAssetValidationError,
        match="source inventory mismatch",
    ):
        assets.validate_book_pdf_assets(root)


def test_noncanonical_output_path_is_rejected(tmp_path: Path) -> None:
    root = fixture_root(tmp_path)
    manifest_path = root / Path(assets.MANIFEST_RELATIVE)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["assets"][0]["rendering"] = "assets/book_pdf_renderings/wrong.pdf"
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(assets.BookAssetValidationError, match="rendering must be"):
        assets.validate_book_pdf_assets(root)


def test_stray_pdf_is_rejected(tmp_path: Path) -> None:
    root = fixture_root(tmp_path)
    stray = root / "assets/book_pdf_renderings/stray.pdf"
    shutil.copyfile(root / "assets/book_pdf_renderings/book-cover.pdf", stray)
    with pytest.raises(
        assets.BookAssetValidationError,
        match="output inventory mismatch",
    ):
        assets.validate_book_pdf_assets(root)
