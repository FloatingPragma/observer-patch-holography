#!/usr/bin/env python3
"""Smoke-test the quark diagonal gap-shift tau-map artifact."""

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
SCALAR_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_scalar_evaluator.py"
TAU_MAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_tau_map.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_tau_map.json"


def main() -> int:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_EMISSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_VALUES_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCALAR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(TAU_MAP_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_family_excitation_diagonal_gap_shift_tau_map":
        print("wrong quark diagonal gap-shift tau-map artifact id", file=sys.stderr)
        return 1
    if payload.get("tau_u_log_per_side") is not None or payload.get("tau_d_log_per_side") is not None:
        print("tau-map coefficients should remain unset until emitted from OPH inputs", file=sys.stderr)
        return 1
    if payload.get("scalar_evaluator_artifact") != "oph_family_excitation_diagonal_gap_shift_scalar_evaluator":
        print("tau-map should reference the scalar evaluator artifact", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "source_readback_u_log_per_side_and_source_readback_d_log_per_side":
        print("tau-map should now identify the emitted pure-B payload pair as the next predictive object on the active builder path", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
