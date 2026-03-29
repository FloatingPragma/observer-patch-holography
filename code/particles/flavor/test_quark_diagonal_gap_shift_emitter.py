#!/usr/bin/env python3
"""Smoke-test the quark diagonal gap-shift emitter artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
MAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_map.py"
SOURCE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_law.py"
SOURCE_READBACK_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_readback.py"
SOURCE_EMISSION_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_emission.py"
SOURCE_VALUES_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_values.py"
TAU_MAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_tau_map.py"
EMITTER_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_emitter.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_emitter.json"


def main() -> int:
    subprocess.run([sys.executable, str(MAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_EMISSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_VALUES_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(TAU_MAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EMITTER_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_family_excitation_diagonal_gap_shift_emitter":
        print("wrong quark diagonal gap-shift emitter artifact id", file=sys.stderr)
        return 1
    if payload.get("tau_u_log_per_side") is not None or payload.get("tau_d_log_per_side") is not None:
        print("predictive tau slots should remain unset until emitted from OPH inputs", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "beta_u_diag_B_source_and_beta_d_diag_B_source":
        print("diagonal gap-shift emitter should point to the beta-pair amplitudes as the next predictive object", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
