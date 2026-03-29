#!/usr/bin/env python3
"""Smoke-test the quark diagonal source-values artifact."""

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
SOURCE_EMISSION_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_emission.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_values.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_values.json"


def main() -> int:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_EMISSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_family_excitation_diagonal_common_gap_shift_source_values":
        print("wrong quark diagonal source-values artifact id", file=sys.stderr)
        return 1
    if payload.get("beta_u_diag_B_source") is not None or payload.get("beta_d_diag_B_source") is not None:
        print("quark diagonal source values should remain unset until the B-mode emitter is populated", file=sys.stderr)
        return 1
    if payload.get("source_emission_artifact") != "oph_family_excitation_diagonal_common_gap_shift_source_emission":
        print("quark diagonal source-values artifact should consume the source-emission layer", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "source_readback_u_log_per_side_and_source_readback_d_log_per_side":
        print("quark diagonal source-values artifact should point to the emitted pure-B payload pair", file=sys.stderr)
        return 1
    if payload.get("source_readback_artifact") != "oph_family_excitation_diagonal_common_gap_shift_source_readback":
        print("quark diagonal source-values artifact should reference the source-readback layer", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
