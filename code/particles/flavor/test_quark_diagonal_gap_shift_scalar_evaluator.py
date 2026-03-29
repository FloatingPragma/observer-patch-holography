#!/usr/bin/env python3
"""Smoke-test the quark diagonal gap-shift scalar-evaluator artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
MAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_map.py"
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
SOURCE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_law.py"
SOURCE_READBACK_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_readback.py"
SOURCE_EMISSION_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_emission.py"
SOURCE_VALUES_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_values.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_scalar_evaluator.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_scalar_evaluator.json"


def main() -> int:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_EMISSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_VALUES_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_family_excitation_diagonal_gap_shift_scalar_evaluator":
        print("wrong quark diagonal scalar-evaluator artifact id", file=sys.stderr)
        return 1
    if payload.get("tau_u_log_per_side") is not None or payload.get("tau_d_log_per_side") is not None:
        print("quark scalar evaluator should remain unset until source values are emitted", file=sys.stderr)
        return 1
    if payload.get("source_values_artifact") != "oph_family_excitation_diagonal_common_gap_shift_source_values":
        print("quark scalar evaluator should reference the diagonal common gap-shift source values", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "tau_u_log_per_side_and_tau_d_log_per_side":
        print("quark scalar evaluator should point to the tau-pair as the next predictive object on the active builder path", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
