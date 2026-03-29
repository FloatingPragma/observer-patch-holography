#!/usr/bin/env python3
"""Smoke-test the charged-lepton excitation-gap map artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_source_emission.py"
VALUE_LAW_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_value_law.py"
SUPPORT_EXTENSION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_minimal_source_support_extension_emitter.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_excitation_gap_map.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_excitation_gap_map.json"


def main() -> int:
    subprocess.run([sys.executable, str(SOURCE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(VALUE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SUPPORT_EXTENSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_charged_lepton_excitation_gap_map":
        print("wrong lepton excitation-gap map artifact id", file=sys.stderr)
        return 1
    if payload.get("a_e_log_coeff") is None or payload.get("b_e_log_coeff") is None:
        print("charged gap map should expose the current source-side coefficient readback", file=sys.stderr)
        return 1
    if payload.get("source_side_ordered_package_value_law_artifact") != "oph_charged_sector_local_ordered_package_value_law":
        print("charged excitation-gap map should reference the ordered-package value-law artifact", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "eta_source_support_extension_log_per_side":
        print("charged excitation-gap map should now point to the support-extension scalar", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
