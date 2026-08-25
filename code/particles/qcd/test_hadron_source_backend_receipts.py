#!/usr/bin/env python3
"""Smoke-test the hadronic source backend receipt scaffold."""

from __future__ import annotations

import hashlib
import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "qcd" / "build_hadron_source_backend_receipts.py"


def test_build_hadron_source_backend_receipts(tmp_path: pathlib.Path) -> None:
    out_dir = tmp_path / "hadron_source_backend"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--output",
            str(out_dir),
            "--claim",
            "SOURCE_PROTOTYPE_NOT_PROMOTED",
            "--tier",
            "H2",
        ],
        check=True,
        cwd=ROOT,
    )

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["milestone"] == "HVP_ALPHA_SOURCE_PROTOTYPE"
    assert manifest["claim"] == "SOURCE_PROTOTYPE_NOT_PROMOTED"
    assert manifest["claim_tier"] == "H2"
    assert manifest["promotion_allowed"] is False
    assert manifest["missing_files"] == []

    for rel_path in manifest["required_files"]:
        assert (out_dir / rel_path).is_file(), rel_path
    assert "manifest.json" not in manifest["file_hashes"]
    for rel_path, expected in manifest["file_hashes"].items():
        actual = "sha256:" + hashlib.sha256((out_dir / rel_path).read_bytes()).hexdigest()
        assert actual == expected, rel_path

    source_dag = json.loads((out_dir / "source_dag.json").read_text(encoding="utf-8"))
    assert source_dag["comparison_data_excluded"] is True
    assert source_dag["no_target_leak_status"] == "PASS_EMPTY_COMPARISON_DAG"
    assert set(source_dag["forbidden_targets"]) == set(manifest["forbidden_targets"])

    claim = (out_dir / "claim.md").read_text(encoding="utf-8").strip()
    assert claim == "SOURCE_PROTOTYPE_NOT_PROMOTED"

    j24 = json.loads((out_dir / "spectral" / "J24Q.json").read_text(encoding="utf-8"))
    assert j24["status"] == "MUST_DERIVE_FROM_MOMENTS_OR_LANCZOS"


def test_rejects_unknown_claim(tmp_path: pathlib.Path) -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--output",
            str(tmp_path / "bad"),
            "--claim",
            "SOURCE_ONLY_BUT_MAGIC",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    assert result.returncode != 0
