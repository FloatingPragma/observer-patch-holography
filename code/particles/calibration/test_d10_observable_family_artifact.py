#!/usr/bin/env python3
"""Validate the D10 electroweak observable-family artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_observable_family.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"


def test_d10_running_family_artifact_is_scheme_clean() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d10_ew_observable_family"
    assert payload["observable_family_id"] == "d10_running_tree"
    assert payload["promotion_gate"]["mixed_scheme"] is False
    assert payload["mixed_reporting_surface"]["mixed_scheme"] is True
    assert payload["mixed_reporting_surface"]["family_purity_violation"] is True
    assert payload["mixed_reporting_surface"]["alpha_em_row_current_source"] == "d10_running_tree"
    assert payload["mixed_reporting_surface"]["exact_missing_transport_object_for_pole_family"] == "EWTransportKernel_D10"
    assert abs(payload["identity_residuals"]["mw_from_v_alpha2"]) < 1.0e-12
    assert abs(payload["identity_residuals"]["alpha_em_from_alpha1_alpha2"]) < 1.0e-12
    assert payload["coherence_witness"]["mixed_sources_detected"] is True
    assert abs(payload["coherence_witness"]["running_mass_ratio_residual"]) < 1.0e-12
    assert abs(payload["coherence_witness"]["stage3_mass_ratio_residual"]) > 1.0e-3
