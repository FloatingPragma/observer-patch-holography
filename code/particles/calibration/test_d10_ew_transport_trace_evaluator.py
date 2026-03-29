#!/usr/bin/env python3
"""Smoke-test the D10 transport-trace evaluator artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
READOUT_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_readout.py"
EVALUATOR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_transport_trace_evaluator.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_transport_trace_evaluator.json"


def main() -> int:
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EVALUATOR_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_d10_ew_transport_trace_evaluator":
        print("wrong D10 transport-trace evaluator artifact id", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "EWSinglePostTransportTreeIdentity_D10":
        print("D10 transport-trace evaluator should now point at the unsplit post-transport tree identity", file=sys.stderr)
        return 1
    if payload.get("insufficiency_certificate", {}).get("sigma_not_selected_by_current_family") is not True:
        print("D10 transport-trace evaluator should carry the current-family underdetermination certificate", file=sys.stderr)
        return 1
    if payload.get("diagnostic_only") is not True:
        print("D10 transport-trace evaluator should now be marked diagnostic-only", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
