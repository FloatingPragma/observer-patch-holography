#!/usr/bin/env python3
"""Smoke-test the charged midpoint-defect emitter artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_source_emission.py"
VALUE_LAW_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_value_law.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_midpoint_defect_emitter.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_midpoint_defect_emitter.json"


def main() -> int:
    subprocess.run([sys.executable, str(SOURCE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(VALUE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_charged_sector_local_ordered_package_midpoint_defect_emitter":
        print("wrong charged midpoint-defect emitter artifact id", file=sys.stderr)
        return 1
    if payload.get("proof_status") != "current_support_midpoint_defect_emitter_closed":
        print("charged midpoint-defect emitter should now be closed on the current support", file=sys.stderr)
        return 1
    if abs(float(payload.get("delta_midpoint_log_per_side_emitted", 1.0))) > 1.0e-12:
        print("current midpoint-defect emitter should expose the presently collapsed zero defect", file=sys.stderr)
        return 1
    if payload.get("delta_midpoint_zero_on_current_support") is not True:
        print("charged midpoint-defect emitter should certify zero defect on the current support", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") is not None:
        print("charged midpoint-defect emitter should no longer identify itself as the remaining predictive object", file=sys.stderr)
        return 1
    if payload.get("next_single_residual_object") != "oph_charged_sector_local_current_support_obstruction_certificate":
        print("charged midpoint-defect emitter should point to the current-support obstruction certificate next", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
