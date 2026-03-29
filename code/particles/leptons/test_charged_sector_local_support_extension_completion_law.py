#!/usr/bin/env python3
"""Validate the charged support-extension completion-law shell."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
VALUE_LAW_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_value_law.py"
OBSTRUCTION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_current_support_obstruction_certificate.py"
MINIMAL_EXTENSION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_minimal_source_support_extension_emitter.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_completion_law.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_completion_law.json"


def test_charged_support_extension_completion_law_is_two_scalar_shell() -> None:
    subprocess.run([sys.executable, str(VALUE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(OBSTRUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MINIMAL_EXTENSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_charged_sector_local_support_extension_completion_law"
    assert payload["proof_status"] == "two_scalar_support_extension_completion_law_closed_source_scalars_open"
    assert payload["smallest_constructive_missing_object"] == "eta_source_support_extension_log_per_side"
    assert payload["next_residual_after_debug_eta_promotion"] == "sigma_source_support_extension_total_log_per_side"
    assert payload["sigma_source_support_extension_total_log_per_side"] is None
    assert payload["eta_source_support_extension_log_per_side"] is None
