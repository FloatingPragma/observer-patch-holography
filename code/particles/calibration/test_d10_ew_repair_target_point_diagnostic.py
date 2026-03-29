#!/usr/bin/env python3
"""Validate the D10 frozen-target-point diagnostic artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_repair_target_point_diagnostic.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_repair_target_point_diagnostic.json"


def test_d10_repair_target_point_is_unique_once_target_spec_is_frozen() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_d10_ew_repair_target_point_diagnostic"
    assert payload["status"] == "diagnostic_only_inverse_target_point"
    assert payload["object_id"] == "EWRepairTargetPointDiagnostic_D10"
    assert payload["spec_id"] == "official_pdg_current_surface_2026_03_28"
    assert payload["tau2_tree_exact_target"] != 0.0
    assert "MW_target^2 / (pi * v_inherited^2 * alpha2_mz) - 1" in payload["formulas"]["tau2_tree_exact_target"]

