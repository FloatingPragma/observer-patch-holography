#!/usr/bin/env python3
"""Validate the charged support-extension endpoint-ratio breaker artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SUPPORT_EXTENSION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_minimal_source_support_extension_emitter.py"
COMPLETION_LAW_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_completion_law.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_endpoint_ratio_breaker.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_endpoint_ratio_breaker.json"


def test_endpoint_ratio_breaker_reduces_support_extension_to_sigma() -> None:
    subprocess.run([sys.executable, str(SUPPORT_EXTENSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(COMPLETION_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_charged_sector_local_support_extension_endpoint_ratio_breaker"
    assert payload["status"] == "closed_smaller_primitive"
    assert payload["precondition_residual_object"] == "eta_source_support_extension_log_per_side"
    assert payload["smallest_constructive_missing_object_within_primitive"] == "sigma_source_support_extension_total_log_per_side"
