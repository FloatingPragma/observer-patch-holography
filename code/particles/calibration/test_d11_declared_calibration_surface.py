#!/usr/bin/env python3
"""Guard the declared D11 calibration surface artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_declared_calibration_surface.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_declared_calibration_surface.json"


def test_d11_declared_calibration_surface_exports_core_and_jacobian() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_declared_calibration_surface"
    assert payload["proof_status"] == "declared_d10_d11_calibration_surface_fixed"
    assert payload["surface_kind"] == "declared_running_matching_threshold_surface"
    assert payload["core"]["y_t_core_mt"] == 0.92046435
    assert payload["jacobian"]["d_mH_d_lambda"] == 480.0
