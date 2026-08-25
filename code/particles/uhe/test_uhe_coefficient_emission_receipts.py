#!/usr/bin/env python3
"""Smoke-test the high-energy messenger coefficient-emission receipt scaffold."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "uhe" / "build_uhe_coefficient_emission_receipts.py"


def test_build_uhe_coefficient_emission_receipts(tmp_path: pathlib.Path) -> None:
    out_dir = tmp_path / "uhe"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(out_dir)],
        check=True,
        cwd=ROOT,
    )

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["milestone"] == "UHE_COEFFICIENT_EMISSION_AUDIT"
    assert manifest["strongest_allowed_claim"] == "SYNTHETIC_COEFFICIENT_SCAFFOLD_ONLY"
    assert manifest["first_blocked_gate"] == "SOURCE_ARTIFACT_AND_COMPUTED_READINESS_RECEIPTS"
    assert manifest["artifact_type"] == "UHE_COEFFICIENT_EMISSION_SCAFFOLD"
    assert manifest["promotion_allowed"] is False
    assert manifest["synthetic_fixture"] is True
    assert manifest["source_artifact_bound"] is False
    assert manifest["physical_claim"] is False
    assert manifest["missing_files"] == []
    assert "manifest.json" not in manifest["file_hashes"]

    for rel_path in manifest["required_files"]:
        assert (out_dir / rel_path).is_file(), rel_path

    coeffs = json.loads((out_dir / "emitted_coefficients.json").read_text(encoding="utf-8"))
    assert coeffs["common_source_lock"] is False
    assert coeffs["coefficient_status"] == "PLANTED_EXAMPLES_NOT_SOURCE_DERIVED"
    assert coeffs["coefficients"]["finite_maxent_eta"] == [0.25, 0.35, -0.20, 0.10, 0.05]

    dag = json.loads((out_dir / "source_dag.json").read_text(encoding="utf-8"))
    assert dag["status"] == "PASS_NO_FORBIDDEN_INPUTS_SCAFFOLD_ONLY"
    assert dag["readiness_gates"]["NO_UHE_DATA_USE"] is True
    assert all(
        value is False
        for name, value in dag["readiness_gates"].items()
        if name != "NO_UHE_DATA_USE"
    )


def test_rejects_uhe_event_data_as_source_input(tmp_path: pathlib.Path) -> None:
    config = tmp_path / "source_config.json"
    config.write_text(json.dumps({"inputs": ["event_coordinates", "likelihood_values"]}), encoding="utf-8")
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
    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert dag["status"] == "FAIL_FORBIDDEN_SOURCE_INPUT"
    assert dag["readiness_gates"]["NO_UHE_DATA_USE"] is False
    assert manifest["strongest_allowed_claim"] == "INVALIDATED_COEFFICIENT_DAG"
    assert manifest["promotion_allowed"] is False
    assert "event_coordinates" in dag["target_leak_hits"]
    assert "likelihood_values" in manifest["target_leak_hits"]
