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
PUBLIC_SOURCE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_d12_public_source_payload.py"
SOURCE_READBACK_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_readback.py"
SOURCE_EMISSION_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_emission.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_common_gap_shift_source_values.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_values.json"


def main() -> int:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MAP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PUBLIC_SOURCE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_EMISSION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_family_excitation_diagonal_common_gap_shift_source_values":
        print("wrong quark diagonal source-values artifact id", file=sys.stderr)
        return 1
    if payload.get("beta_u_diag_B_source") is None or payload.get("beta_d_diag_B_source") is None:
        print("quark diagonal source values should be populated on the selected public class", file=sys.stderr)
        return 1
    if payload.get("source_emission_artifact") != "oph_family_excitation_diagonal_common_gap_shift_source_emission":
        print("quark diagonal source-values artifact should consume the source-emission layer", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "source_derived_c_d_over_c_u_or_t1_value":
        print("quark diagonal source-values artifact should retain the source-derived scalar gap", file=sys.stderr)
        return 1
    if payload.get("source_readback_artifact") != "oph_family_excitation_diagonal_common_gap_shift_source_readback":
        print("quark diagonal source-values artifact should reference the source-readback layer", file=sys.stderr)
        return 1
    if payload.get("proof_status") != "comparison_only_values_from_unpromotable_source_projection":
        print("target-attached values must remain comparison-only", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
