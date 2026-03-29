#!/usr/bin/env python3
"""Validate the scalarized quadratic-even D12 quark transport shell."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_quadratic_even_transport_scalar.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_quadratic_even_transport_scalar.json"


def main() -> int:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_quark_quadratic_even_transport_scalar":
        print("unexpected artifact id", file=sys.stderr)
        return 1
    if payload.get("next_single_residual_object") != "eta_Q_centered_value_law":
        print("quadratic shell should reduce to eta_Q_centered_value_law", file=sys.stderr)
        return 1
    if payload.get("quadratic_even_log_formula_direct") != "(eta_Q_centered / 6) * (1, -2, 1)":
        print("direct quadratic-even log formula mismatch", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
