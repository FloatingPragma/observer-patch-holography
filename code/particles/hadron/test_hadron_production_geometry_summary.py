#!/usr/bin/env python3
"""Validate the hadron production geometry summary."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
PAYLOAD_SCRIPT = ROOT / "particles" / "hadron" / "derive_stable_channel_cfg_source_measure_payload.py"
RECEIPT_SCRIPT = ROOT / "particles" / "hadron" / "derive_runtime_schedule_receipt_n_therm_and_n_sep.py"
SCRIPT = ROOT / "particles" / "hadron" / "derive_hadron_production_geometry_summary.py"
OUTPUT = ROOT / "particles" / "runs" / "hadron" / "production_geometry_summary.json"


def main() -> int:
    subprocess.run([sys.executable, str(PAYLOAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(RECEIPT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PAYLOAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_hadron_production_geometry_summary":
        print("unexpected production geometry summary artifact", file=sys.stderr)
        return 1
    totals = payload.get("totals") or {}
    if totals.get("n_ensembles") != 3 or totals.get("total_cfg") != 6:
        print("production geometry summary should expose the frozen 3-ensemble / 6-cfg family", file=sys.stderr)
        return 1
    if totals.get("total_raw_gauge_bytes_all_cfg_naive", 0) <= totals.get("total_correlator_bytes_float64_backend_dump", 0):
        print("gauge storage should dominate backend correlator dump size", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
