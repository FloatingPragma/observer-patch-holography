#!/usr/bin/env python3
"""Smoke-test the minimal charged support-extension emitter artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
VALUE_LAW_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_value_law.py"
OBSTRUCTION_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_current_support_obstruction_certificate.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_minimal_source_support_extension_emitter.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_minimal_source_support_extension_emitter.json"


def main() -> int:
    subprocess.run([sys.executable, str(VALUE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(OBSTRUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_charged_sector_local_minimal_source_support_extension_emitter":
        print("wrong charged support-extension emitter artifact id", file=sys.stderr)
        return 1
    if payload.get("proof_status") != "minimal_support_extension_formula_closed_source_scalar_open":
        print("charged support-extension emitter should close the formula but keep the source scalar open", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "eta_source_support_extension_log_per_side":
        print("charged support-extension emitter should reduce to eta_source_support_extension_log_per_side", file=sys.stderr)
        return 1
    if payload.get("eta_source_support_extension_log_per_side") is not None or payload.get("kappa_ext") is not None:
        print("charged support-extension source scalar should remain unset", file=sys.stderr)
        return 1
    if payload.get("eta_source_support_extension_log_per_side_candidate") is None:
        print("charged support-extension emitter should expose the rigid ordered-ratio eta candidate", file=sys.stderr)
        return 1
    if payload.get("candidate_next_single_residual_object") != "sigma_source_support_extension_total_log_per_side":
        print("charged support-extension emitter should expose the next total-span residual after the eta candidate", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
