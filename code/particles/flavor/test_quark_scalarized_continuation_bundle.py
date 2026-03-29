#!/usr/bin/env python3
"""Validate the scalarized D12 quark continuation bundle."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
QUADRATIC_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_quadratic_even_transport_scalar.py"
PHYSICAL_SCRIPT = ROOT / "particles" / "flavor" / "derive_generation_bundle_same_label_physical_invariant_bundle.py"
BUNDLE_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_scalarized_continuation_bundle.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_scalarized_continuation_bundle.json"


def main() -> int:
    subprocess.run([sys.executable, str(QUADRATIC_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PHYSICAL_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(BUNDLE_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_quark_scalarized_continuation_bundle":
        print("unexpected scalarized continuation bundle artifact", file=sys.stderr)
        return 1
    if payload.get("proof_status") != "scalarized_bundle_not_OPH_derived":
        print("scalarized continuation bundle should stay diagnostic", file=sys.stderr)
        return 1
    residuals = payload.get("honest_remaining_value_laws") or []
    if "Delta_ud_overlap_value_law" not in residuals or "same_label_left_transport_physical_invariant_value_laws" not in residuals:
        print("scalarized continuation bundle should expose the remaining D12 value laws", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
