#!/usr/bin/env python3
"""Smoke-test the JWST compact-object source-release receipt scaffold."""

from __future__ import annotations

import hashlib
import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "jwst" / "build_compact_object_source_release_receipts.py"


def test_build_compact_object_source_release_receipts(tmp_path: pathlib.Path) -> None:
    out_dir = tmp_path / "jwst"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(out_dir)],
        check=True,
        cwd=ROOT,
    )

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["milestone"] == "Q3_JWST_COMPACT_OBJECT_SOURCE_RELEASE_AUDIT"
    assert manifest["strongest_allowed_claim"] == "J0_DIAGNOSTIC_PROXY"
    assert manifest["first_blocked_gate"] == "CATALOG_INGESTION_RECEIPT"
    assert manifest["promotion_allowed"] is False
    assert manifest["missing_files"] == []

    for rel_path in manifest["required_files"]:
        assert (out_dir / rel_path).is_file(), rel_path
    assert "manifest.json" not in manifest["file_hashes"]
    for rel_path, expected in manifest["file_hashes"].items():
        actual = "sha256:" + hashlib.sha256((out_dir / rel_path).read_bytes()).hexdigest()
        assert actual == expected, rel_path

    source = json.loads((out_dir / "source_artifact.json").read_text(encoding="utf-8"))
    assert source["readiness_gates"]["NO_TARGET_LEAKAGE_RECEIPT"] is True
    assert source["readiness_gates"]["OBJECT_SOURCE_LAW_RECEIPT"] is False

    ladder = json.loads((out_dir / "claim_ladder.json").read_text(encoding="utf-8"))
    assert ladder["strongest_allowed_claim"] == "J0_DIAGNOSTIC_PROXY"
    assert "luminosity is not stellar mass" in ladder["nonclaims"]

    claim = (out_dir / "claim.md").read_text(encoding="utf-8").strip()
    assert claim == "J0_DIAGNOSTIC_PROXY"


def test_rejects_jwst_catalog_as_source_input(tmp_path: pathlib.Path) -> None:
    config = tmp_path / "source_config.json"
    config.write_text(json.dumps({"source_inputs": ["jwst_catalog_counts"]}), encoding="utf-8")
    out_dir = tmp_path / "leaky"

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--output",
            str(out_dir),
            "--config",
            str(config),
        ],
        check=True,
        cwd=ROOT,
    )

    dag = json.loads((out_dir / "no_target_leak_dag.json").read_text(encoding="utf-8"))
    assert dag["NO_TARGET_LEAKAGE_RECEIPT"] is False
    assert dag["status"] == "FAIL_FORBIDDEN_SOURCE_INPUT"
    assert "jwst_catalog" in dag["target_leak_hits"]
    assert "catalog_counts" in dag["target_leak_hits"]
