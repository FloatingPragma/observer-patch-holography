#!/usr/bin/env python3
"""Ensure the D11 predictive artifact does not export reference-fit witnesses."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DECLARED_SURFACE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_declared_calibration_surface.py"
SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_critical_surface_readout.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_critical_surface_readout.json"


def test_d11_predictive_artifact_omits_reference_exact_witness() -> None:
    subprocess.run([sys.executable, str(DECLARED_SURFACE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert "current_reference_exact_witness" not in payload["readout_kernel"]
    assert "current_reference_exact_witness" not in payload
