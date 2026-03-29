#!/usr/bin/env python3
"""Validate the charged eta source-readback primitive."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
VALUE_LAW_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_value_law.py"
MINIMAL_EXTENSION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_minimal_source_support_extension_emitter.py"
COMPLETION_LAW_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_completion_law.py"
ENDPOINT_RATIO_BREAKER_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_endpoint_ratio_breaker.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_eta_source_readback.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_eta_source_readback.json"


def test_eta_source_readback_reduces_charged_support_shell_to_sigma() -> None:
    subprocess.run([sys.executable, str(VALUE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MINIMAL_EXTENSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(COMPLETION_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ENDPOINT_RATIO_BREAKER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_charged_sector_local_support_extension_eta_source_readback"
    assert payload["status"] == "closed_smaller_primitive"
    assert payload["proof_status"] == "eta_source_readback_formula_closed_value_open"
    assert payload["eta_source_support_extension_log_per_side"] is None
    assert payload["next_single_residual_object_after_eta"] == "sigma_source_support_extension_total_log_per_side"
    assert payload["downstream_sigma_artifact"] == "oph_charged_sector_local_support_extension_endpoint_ratio_breaker"
