#!/usr/bin/env python3
"""Smoke-test the quark diagonal common gap-shift source-readback artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
MAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_map.py"
SOURCE_LAW_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_law.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_readback.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_readback.json"


def main() -> int:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_family_excitation_diagonal_common_gap_shift_source_readback":
        print("wrong quark diagonal source-readback artifact id", file=sys.stderr)
        return 1
    if payload.get("proof_status") != "source_readback_law_closed_waiting_pure_B_payload_pair":
        print("quark source-readback law should now be closed and wait on the pure-B payload pair", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "source_readback_u_log_per_side_and_source_readback_d_log_per_side":
        print("quark source-readback artifact should now reduce to the emitted pure-B payload pair", file=sys.stderr)
        return 1
    if payload.get("first_data_bearing_primitive_beneath_scalar_pair") != "source_readback_u_log_per_side_and_source_readback_d_log_per_side":
        print("quark source-readback artifact should expose the pure-B payload pair as the first data-bearing primitive", file=sys.stderr)
        return 1
    if payload.get("J_B_functional_kind") != "pure_B_odd_point_separating_projection":
        print("quark source-readback artifact should expose the pure-B odd projector", file=sys.stderr)
        return 1
    if payload.get("source_readback_u_log_per_side") is not None or payload.get("source_readback_d_log_per_side") is not None:
        print("quark source-readback arrays should remain unset until a predictive readback is closed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
