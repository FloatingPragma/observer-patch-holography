#!/usr/bin/env python3
"""Generate the post-integration manifest for the hardened package."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


SOURCE_ARCHIVE_SHA256 = "fb871c5ac810f8703fd49b9fdcf621096c1bd60bc0c426d34ba42dc2fd3bae8c"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_manifest(root: Path, output: Path) -> dict[str, object]:
    files: dict[str, dict[str, int | str]] = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if (
            not path.is_file()
            or path == output
            or "__pycache__" in relative.parts
            or path.suffix == ".pyc"
        ):
            continue
        files[relative.as_posix()] = {
            "bytes": path.stat().st_size,
            "sha256": digest(path),
        }
    return {
        "schema": "physical_wz_strict_one_loop_rer_port_manifest_v1",
        "package_id": "physical_wz_strict_one_loop_pole_map_rer_port_2026-07-20",
        "created": "2026-07-20",
        "source_archive_sha256": SOURCE_ARCHIVE_SHA256,
        "upstream_manifest": "UPSTREAM_PACKAGE_MANIFEST.json",
        "upstream_manifest_validation": "PASS_12_OF_12_BEFORE_INTEGRATION",
        "port_location": "code/particles/calibration/strict_one_loop_pole_map",
        "scientific_classification": "CONDITIONAL_STRICT_1L_POLE_MAP_NOT_OPH_NATIVE_PHYSICAL",
        "physical_promotion_allowed": False,
        "file_count": len(files),
        "files": files,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = root / "INTEGRATION_MANIFEST.json"
    manifest = build_manifest(root, output)
    if args.check:
        stored = json.loads(output.read_text(encoding="utf-8"))
        if stored != manifest:
            raise SystemExit("FAIL: integration manifest drift")
        print(f"PASS: integration manifest {manifest['file_count']}/{manifest['file_count']}")
        return 0
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
