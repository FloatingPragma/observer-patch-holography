#!/usr/bin/env python3
"""Guard the exact compare-only neutrino bridge-coordinate sidecar."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
ADAPTER_SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_two_parameter_exact_adapter.py"
CORRECTION_SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_bridge_correction_invariant_scaffold.py"
SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_exact_adapter_bridge_coordinate.py"
OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_exact_adapter_bridge_coordinate.json"


def test_neutrino_exact_adapter_bridge_coordinate() -> None:
    subprocess.run([sys.executable, str(ADAPTER_SCRIPT)], check=True, capture_output=True, text=True)
    subprocess.run([sys.executable, str(CORRECTION_SCRIPT)], check=True, capture_output=True, text=True)
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(OUTPUT)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "saved:" in completed.stdout
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_neutrino_exact_adapter_bridge_coordinate"
    assert payload["scope"] == "compare_only_two_parameter_segment_adapter"
    assert payload["promotable"] is False
    bridge = payload["bridge_coordinates"]
    # Every bridge coordinate is built from tau_nu, lambda_nu and p_nu of the
    # exact adapter, so it inherits that fit's resolution: the ratio it solves
    # is a difference of np.linalg.eigh eigenvalues, whose binary64 noise floor
    # puts tau_nu at about 1e-12 and these values a few times 1e-15 apart
    # between builds. Assert relative agreement at 1e-12 instead of an absolute
    # 1e-15 that the arithmetic cannot deliver.
    assert bridge["paper_facing_amplitude"]["value"] == pytest.approx(6.696759752761619, rel=1.0e-12)
    assert bridge["internal_proxy"]["value"] == pytest.approx(6.699912221878757, rel=1.0e-12)
    assert bridge["reduced_correction_invariant"]["value"] == pytest.approx(0.9995294760568888, rel=1.0e-12)
    checks = payload["consistency_checks"]
    assert checks["reconstruction_error_abs"] == pytest.approx(0.0, abs=1.0e-15)
    assert checks["within_current_compare_only_c_window"] is True
    assert checks["within_current_compare_only_b_window"] is True
