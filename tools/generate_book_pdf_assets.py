#!/usr/bin/env python3
"""Regenerate committed canonical PDF renderings for the book's SVG assets."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

from book_pdf_assets import (
    MANIFEST_RELATIVE,
    MANIFEST_SCHEMA_VERSION,
    RENDERING_ROOT_RELATIVE,
    REPO_ROOT,
    discover_book_svg_sources,
    rendering_relative_for,
    sha256_file,
    validate_book_pdf_assets,
)
from reproducible_build_env import build_environment


BUILD_DIR = REPO_ROOT.parent / "temp" / "book_pdf_asset_generation"


def ensure_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"required tool not found in PATH: {name}")


def run(command: list[str]) -> None:
    subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=True,
        env=build_environment(REPO_ROOT),
    )


def reusable_renderings() -> dict[str, Path]:
    """Return existing outputs whose paths and source/output hashes still bind."""
    manifest_path = REPO_ROOT / Path(MANIFEST_RELATIVE)
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(payload, dict) or not isinstance(payload.get("assets"), list):
        return {}

    sources_by_name = {
        source.as_posix(): source for source in discover_book_svg_sources()
    }
    reusable: dict[str, Path] = {}
    for record in payload["assets"]:
        if not isinstance(record, dict):
            continue
        source_name = record.get("source")
        if not isinstance(source_name, str):
            continue
        source_relative = sources_by_name.get(source_name)
        if source_relative is None:
            continue
        expected_rendering = rendering_relative_for(source_relative)
        if record.get("rendering") != expected_rendering.as_posix():
            continue
        source_path = REPO_ROOT / Path(source_relative)
        rendering_path = REPO_ROOT / Path(expected_rendering)
        if not source_path.is_file() or not rendering_path.is_file():
            continue
        if record.get("source_sha256") != sha256_file(source_path):
            continue
        if record.get("rendering_sha256") != sha256_file(rendering_path):
            continue
        reusable[source_name] = rendering_path
    return reusable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all",
        action="store_true",
        help="rerender every SVG instead of preserving unchanged hash-bound outputs",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensure_tool("rsvg-convert")
    ensure_tool("gs")
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    sources = discover_book_svg_sources()
    reusable = {} if args.all else reusable_renderings()
    generated_count = 0

    with tempfile.TemporaryDirectory(
        prefix="renderings-",
        dir=BUILD_DIR,
    ) as temporary:
        staging_root = Path(temporary)
        staged: dict[str, Path] = {}
        for source_relative in sources:
            source_path = REPO_ROOT / Path(source_relative)
            source_name = source_relative.as_posix()
            if source_name in reusable:
                staged[source_name] = reusable[source_name]
                continue
            relative_output = rendering_relative_for(source_relative)
            staged_output = staging_root / Path(relative_output)
            staged_output.parent.mkdir(parents=True, exist_ok=True)
            raw_output = staged_output.with_name(f"{staged_output.stem}.raw.pdf")
            run(
                [
                    "rsvg-convert",
                    "-f",
                    "pdf",
                    "-o",
                    str(raw_output),
                    str(source_path),
                ]
            )
            run(
                [
                    "gs",
                    "-q",
                    "-dSAFER",
                    "-dBATCH",
                    "-dNOPAUSE",
                    "-sDEVICE=pdfwrite",
                    "-dCompatibilityLevel=1.5",
                    f"-sOutputFile={staged_output}",
                    str(raw_output),
                ]
            )
            raw_output.unlink()
            staged[source_name] = staged_output
            generated_count += 1

        rendering_root = REPO_ROOT / Path(RENDERING_ROOT_RELATIVE)
        rendering_root.mkdir(parents=True, exist_ok=True)
        expected_outputs: set[Path] = set()
        manifest_assets: list[dict[str, str]] = []
        for source_relative in sources:
            source_name = source_relative.as_posix()
            rendering_relative = rendering_relative_for(source_relative)
            rendering_path = REPO_ROOT / Path(rendering_relative)
            rendering_path.parent.mkdir(parents=True, exist_ok=True)
            if staged[source_name] != rendering_path:
                shutil.copyfile(staged[source_name], rendering_path)
            expected_outputs.add(rendering_path)
            manifest_assets.append(
                {
                    "source": source_name,
                    "source_sha256": sha256_file(REPO_ROOT / Path(source_relative)),
                    "rendering": rendering_relative.as_posix(),
                    "rendering_sha256": sha256_file(rendering_path),
                }
            )

        for existing in rendering_root.rglob("*.pdf"):
            if existing not in expected_outputs:
                existing.unlink()

        manifest_path = REPO_ROOT / Path(MANIFEST_RELATIVE)
        manifest_path.write_text(
            json.dumps(
                {
                    "schema_version": MANIFEST_SCHEMA_VERSION,
                    "assets": manifest_assets,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    validated = validate_book_pdf_assets()
    print(
        "Regenerated canonical book PDF assets: "
        f"{len(validated)} source-bound renderings "
        f"({generated_count} rendered, {len(validated) - generated_count} preserved)"
    )


if __name__ == "__main__":
    main()
