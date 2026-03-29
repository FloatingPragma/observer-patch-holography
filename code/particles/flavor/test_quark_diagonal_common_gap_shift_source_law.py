#!/usr/bin/env python3
"""Smoke-test the quark diagonal common gap-shift source-law artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
MAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_map.py"
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_law.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_law.json"


def main() -> int:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_family_excitation_diagonal_common_gap_shift_source_law":
        print("wrong quark diagonal common gap-shift source-law artifact id", file=sys.stderr)
        return 1
    if payload.get("proof_status") != "source_law_closed_waiting_J_B_source_pair":
        print("quark diagonal common gap-shift source law should now wait only on the J_B source pair", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "J_B_source_u_and_J_B_source_d":
        print("quark diagonal common gap-shift source law should now point to the J_B source pair", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
