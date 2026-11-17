#!/usr/bin/env python3
"""Validate and stage the book's committed SVG-to-PDF renderings.

SVG text layout depends on the host's librsvg, Cairo, Pango, and font stack.
The book therefore embeds committed PDF renderings instead of converting SVG
files during every publication build.  This module binds each rendering to its
source SVG by path and SHA-256 digest.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
RENDERING_ROOT_RELATIVE = PurePosixPath("assets/book_pdf_renderings")
MANIFEST_RELATIVE = RENDERING_ROOT_RELATIVE / "manifest.json"
MANIFEST_SCHEMA_VERSION = 1
REQUIRED_TOP_LEVEL_SOURCES = (
    PurePosixPath("assets/book-cover.svg"),
    PurePosixPath("assets/pixel-constant.svg"),
    PurePosixPath("assets/OPH_Unification_Diagram.svg"),
)
ASSET_RECORD_KEYS = {
    "source",
    "source_sha256",
    "rendering",
    "rendering_sha256",
}


class BookAssetValidationError(RuntimeError):
    """Raised when the committed rendering inventory is incomplete or stale."""


@dataclass(frozen=True)
class ValidatedBookAsset:
    source_relative: PurePosixPath
    source_path: Path
    rendering_relative: PurePosixPath
    rendering_path: Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover_book_svg_sources(repo_root: Path = REPO_ROOT) -> tuple[PurePosixPath, ...]:
    """Return the exact SVG inventory consumed by the book builder."""
    paths = list(REQUIRED_TOP_LEVEL_SOURCES)
    diagram_root = repo_root / "assets" / "book_diagrams"
    paths.extend(
        PurePosixPath(path.relative_to(repo_root).as_posix())
        for path in sorted(diagram_root.glob("*.svg"))
    )
    return tuple(paths)


def rendering_relative_for(source_relative: PurePosixPath) -> PurePosixPath:
    if not source_relative.is_relative_to(PurePosixPath("assets")):
        raise ValueError(f"book SVG source is outside assets/: {source_relative}")
    inside_assets = source_relative.relative_to(PurePosixPath("assets"))
    return RENDERING_ROOT_RELATIVE / inside_assets.with_suffix(".pdf")


def _load_manifest(path: Path) -> dict[str, Any]:
    if path.is_symlink():
        raise BookAssetValidationError(
            f"canonical book rendering manifest must not be a symlink: {path}"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise BookAssetValidationError(
            f"missing canonical book rendering manifest: {path}"
        ) from exc
    except (OSError, json.JSONDecodeError) as exc:
        raise BookAssetValidationError(
            f"cannot read canonical book rendering manifest {path}: {exc}"
        ) from exc
    if not isinstance(payload, dict):
        raise BookAssetValidationError("book rendering manifest root must be an object")
    if set(payload) != {"schema_version", "assets"}:
        raise BookAssetValidationError(
            "book rendering manifest must contain exactly schema_version and assets"
        )
    if payload["schema_version"] != MANIFEST_SCHEMA_VERSION:
        raise BookAssetValidationError(
            "unsupported book rendering manifest schema_version: "
            f"{payload['schema_version']!r}"
        )
    if not isinstance(payload["assets"], list):
        raise BookAssetValidationError("book rendering manifest assets must be an array")
    return payload


def validate_book_pdf_assets(
    repo_root: Path = REPO_ROOT,
) -> tuple[ValidatedBookAsset, ...]:
    """Validate exact source/output coverage and every recorded digest."""
    repo_root = repo_root.resolve()
    manifest_path = repo_root / Path(MANIFEST_RELATIVE)
    payload = _load_manifest(manifest_path)
    expected_sources = discover_book_svg_sources(repo_root)
    expected_source_names = {path.as_posix() for path in expected_sources}

    records_by_source: dict[str, dict[str, Any]] = {}
    for index, raw_record in enumerate(payload["assets"]):
        if not isinstance(raw_record, dict):
            raise BookAssetValidationError(
                f"book rendering manifest asset {index} must be an object"
            )
        if set(raw_record) != ASSET_RECORD_KEYS:
            raise BookAssetValidationError(
                f"book rendering manifest asset {index} must contain exactly "
                f"{sorted(ASSET_RECORD_KEYS)}"
            )
        if not all(isinstance(raw_record[key], str) for key in ASSET_RECORD_KEYS):
            raise BookAssetValidationError(
                f"book rendering manifest asset {index} fields must be strings"
            )
        source_name = raw_record["source"]
        if source_name in records_by_source:
            raise BookAssetValidationError(
                f"duplicate book SVG source in rendering manifest: {source_name}"
            )
        records_by_source[source_name] = raw_record

    actual_source_names = set(records_by_source)
    if actual_source_names != expected_source_names:
        missing = sorted(expected_source_names - actual_source_names)
        extra = sorted(actual_source_names - expected_source_names)
        raise BookAssetValidationError(
            f"book rendering source inventory mismatch; missing={missing}, extra={extra}"
        )

    validated: list[ValidatedBookAsset] = []
    expected_renderings: set[str] = set()
    for source_relative in expected_sources:
        source_name = source_relative.as_posix()
        record = records_by_source[source_name]
        expected_rendering = rendering_relative_for(source_relative)
        expected_rendering_name = expected_rendering.as_posix()
        expected_renderings.add(expected_rendering_name)
        if record["rendering"] != expected_rendering_name:
            raise BookAssetValidationError(
                f"{source_name}: rendering must be {expected_rendering_name}, "
                f"not {record['rendering']!r}"
            )

        source_path = repo_root / Path(source_relative)
        rendering_path = repo_root / Path(expected_rendering)
        if not source_path.is_file() or source_path.is_symlink():
            raise BookAssetValidationError(
                f"{source_name}: source must be a regular tracked file"
            )
        if not rendering_path.is_file() or rendering_path.is_symlink():
            raise BookAssetValidationError(
                f"{expected_rendering_name}: rendering must be a regular tracked file"
            )
        source_digest = sha256_file(source_path)
        if record["source_sha256"] != source_digest:
            raise BookAssetValidationError(
                f"{source_name}: source SHA-256 mismatch; regenerate canonical renderings"
            )
        rendering_digest = sha256_file(rendering_path)
        if record["rendering_sha256"] != rendering_digest:
            raise BookAssetValidationError(
                f"{expected_rendering_name}: rendering SHA-256 mismatch"
            )
        if not rendering_path.read_bytes().startswith(b"%PDF-"):
            raise BookAssetValidationError(
                f"{expected_rendering_name}: canonical rendering is not a PDF"
            )

        validated.append(
            ValidatedBookAsset(
                source_relative=source_relative,
                source_path=source_path,
                rendering_relative=expected_rendering,
                rendering_path=rendering_path,
            )
        )

    rendering_root = repo_root / Path(RENDERING_ROOT_RELATIVE)
    actual_renderings = {
        PurePosixPath(path.relative_to(repo_root).as_posix()).as_posix()
        for path in rendering_root.rglob("*.pdf")
    }
    if actual_renderings != expected_renderings:
        missing = sorted(expected_renderings - actual_renderings)
        extra = sorted(actual_renderings - expected_renderings)
        raise BookAssetValidationError(
            f"book rendering output inventory mismatch; missing={missing}, extra={extra}"
        )

    return tuple(validated)


def stage_book_pdf_assets(
    out_dir: Path,
    repo_root: Path = REPO_ROOT,
) -> dict[str, Path]:
    """Copy validated canonical bytes into the book build tree."""
    converted: dict[str, Path] = {}
    for asset in validate_book_pdf_assets(repo_root):
        inside_assets = asset.source_relative.relative_to(
            PurePosixPath("assets")
        ).with_suffix(".pdf")
        output_path = out_dir / Path(inside_assets)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(asset.rendering_path, output_path)
        original_ref = f"../{asset.source_relative.as_posix()}"
        converted[original_ref] = output_path
    return converted


def main() -> int:
    try:
        assets = validate_book_pdf_assets()
    except BookAssetValidationError as exc:
        print(f"book PDF asset validation failed: {exc}")
        return 1
    print(
        "book PDF asset validation passed: "
        f"{len(assets)} source-bound canonical renderings"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
