#!/usr/bin/env python3
"""Smoke-test the quark diagonal B-odd source scalar evaluator artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
MAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_map.py"
SOURCE_LAW_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_law.py"
SOURCE_READBACK_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_readback.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_B_odd_source_scalar_evaluator.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_B_odd_source_scalar_evaluator.json"


def main() -> int:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_quark_diagonal_B_odd_source_scalar_evaluator":
        print("wrong quark B-odd source scalar evaluator artifact id", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "source_readback_u_log_per_side_and_source_readback_d_log_per_side":
        print("quark B-odd evaluator should reduce to the emitted pure-B payload pair", file=sys.stderr)
        return 1
    if payload.get("J_B_on_B_ord") != 1.0 or payload.get("J_B_on_center_vector") != 0.0 or payload.get("J_B_on_Q_ord") != 0.0:
        print("quark B-odd evaluator should expose the projector normalization certificates", file=sys.stderr)
        return 1
    if payload.get("predictive_J_B_source_law_status") != "missing":
        print("quark B-odd evaluator should still mark the predictive J_B source law as missing", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
