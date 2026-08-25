#!/usr/bin/env python3
"""Smoke-test the gamma morphology receipt scaffold."""

from __future__ import annotations

import hashlib
import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "gamma" / "build_gamma_morphology_receipts.py"


def test_build_gamma_morphology_receipts(tmp_path: pathlib.Path) -> None:
    out_dir = tmp_path / "gamma"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(out_dir)],
        check=True,
        cwd=ROOT,
    )

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["milestone"] == "GAMMA_MORPHOLOGY_AUDIT"
    assert manifest["strongest_allowed_claim"] == "DIAGNOSTIC_GAMMA_MAP"
    assert manifest["first_blocked_gate"] == "GAMMA_SOURCE_ARTIFACT_RECEIPT"
    assert manifest["promotion_allowed"] is False
    assert manifest["missing_files"] == []

    for rel_path in manifest["required_files"]:
        assert (out_dir / rel_path).is_file(), rel_path
    assert "manifest.json" not in manifest["file_hashes"]
    for rel_path, expected in manifest["file_hashes"].items():
        actual = "sha256:" + hashlib.sha256((out_dir / rel_path).read_bytes()).hexdigest()
        assert actual == expected, rel_path

    source_dag = json.loads((out_dir / "source_dag.json").read_text(encoding="utf-8"))
    assert source_dag["status"] == "PASS_EMPTY_COMPARISON_DAG"

    photon = json.loads((out_dir / "photon_response.json").read_text(encoding="utf-8"))
    assert photon["direct_anomaly_gamma_default"] == 0.0
    assert photon["nonzero_direct_requires"] == "ANOMALY_EM_CURRENT_RECEIPT"

    claim = (out_dir / "claim.md").read_text(encoding="utf-8").strip()
    assert claim == "DIAGNOSTIC_GAMMA_MAP"


def test_rejects_gamma_residual_as_source_input(tmp_path: pathlib.Path) -> None:
    config = tmp_path / "source_config.json"
    config.write_text(json.dumps({"source_inputs": ["gamma_residual_maps"]}), encoding="utf-8")
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

    dag = json.loads((out_dir / "source_dag.json").read_text(encoding="utf-8"))
    receipts = json.loads((out_dir / "receipts.json").read_text(encoding="utf-8"))
    assert dag["status"] == "FAIL_FORBIDDEN_SOURCE_INPUT"
    assert "gamma_residual_maps" in dag["target_leak_hits"]
    assert receipts["readiness_gates"]["GAMMA_NO_DATA_USE_RECEIPT"] is False
