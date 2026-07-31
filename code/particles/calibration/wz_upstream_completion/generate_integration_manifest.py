#!/usr/bin/env python3
"""Generate the post-audit manifest without self-referential runtime logs."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "INTEGRATION_MANIFEST.json"
EXCLUDED_PARTS = {".pytest_cache", "__pycache__"}
EXCLUDED_RELATIVE = {
    "INTEGRATION_MANIFEST.json",
    "outputs/latest_validation.log",
}


def included_files() -> list[Path]:
    files: list[Path] = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT).as_posix()
        if rel in EXCLUDED_RELATIVE:
            continue
        if any(part in EXCLUDED_PARTS for part in path.relative_to(ROOT).parts):
            continue
        if path.suffix in {".pyc", ".pyo"}:
            continue
        files.append(path)
    return sorted(files, key=lambda p: p.relative_to(ROOT).as_posix())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build() -> dict:
    records = [
        {
            "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in included_files()
    ]
    return {
        "schema_version": "physical_wz_upstream_integration_manifest_v1",
        "date": "2026-07-20",
        "upstream_archive_sha256": "e6184ec1a535244c0f805c14e275b085a64c7ed078d780828c038d4f21c86a0d",
        "scientific_status": "DRAFT_SUFFICIENCY_STACK_DEFINED__SIMULATION_RECEIPTS_OPEN__NO_OPH_NATIVE_POLE_PROMOTION",
        "promotion_allowed": False,
        "baseline_specification_schema_documents": 9,
        "diagnostic_schema_documents": 1,
        "certified_contour_schema_documents": 2,
        "total_schema_documents": 12,
        "baseline_specification_receipt_instances": 10,
        "wz_boundary_diagnostic_receipts": 6,
        "baseline_exact_symbolic_checks": 8,
        "wz_exact_finite_corrections": 3,
        "baseline_specification_tests": 38,
        "package_test_count_policy": (
            "the current executable count is reported by pytest collection; "
            "38 is the preserved baseline specification-suite count"
        ),
        "runtime_log_excluded": "outputs/latest_validation.log",
        "file_count_excluding_manifest_and_runtime_log": len(records),
        "total_bytes_excluding_manifest_and_runtime_log": sum(
            record["bytes"] for record in records
        ),
        "files": records,
    }


def main() -> None:
    MANIFEST.write_text(
        json.dumps(build(), indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    print(MANIFEST)


if __name__ == "__main__":
    main()
