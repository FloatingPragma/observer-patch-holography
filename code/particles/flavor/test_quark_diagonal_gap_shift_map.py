#!/usr/bin/env python3
"""Smoke-test the next quark family beyond the current closed surface."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
AUDIT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exactness_audit.py"
MAP_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_diagonal_gap_shift_map.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_map.json"


def main() -> int:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(AUDIT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MAP_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_family_excitation_diagonal_gap_shift_map":
        print("wrong diagonal gap-shift artifact id", file=sys.stderr)
        return 1
    if payload.get("surface_exhausted") is not True:
        print("diagonal gap-shift map should only appear after the current surface is exhausted", file=sys.stderr)
        return 1
    if payload.get("tau_u_log_per_side") is not None or payload.get("tau_d_log_per_side") is not None:
        print("predictive tau slots should remain unset until emitted from OPH inputs", file=sys.stderr)
        return 1
    if payload.get("B_ord") != [-1.0, 0.0, 1.0]:
        print("diagonal gap-shift map should expose the B_ord basis", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
