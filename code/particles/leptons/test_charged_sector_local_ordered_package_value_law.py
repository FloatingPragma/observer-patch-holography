#!/usr/bin/env python3
"""Smoke-test the charged ordered-package value-law artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_source_emission.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_value_law.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_value_law.json"


def main() -> int:
    subprocess.run([sys.executable, str(SOURCE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_charged_sector_local_ordered_package_value_law":
        print("wrong charged ordered-package value-law artifact id", file=sys.stderr)
        return 1
    if payload.get("current_package_linear_subray_only") is not True:
        print("current charged ordered package should still sit on the linear subray", file=sys.stderr)
        return 1
    if payload.get("collapse_proven") is not True or payload.get("carrier_centered_rank") != 1:
        print("charged ordered-package value law should certify the current collapse to a rank-one centered package", file=sys.stderr)
        return 1
    if payload.get("predictive_value_law_closed") is not False:
        print("charged ordered-package value law should remain open", file=sys.stderr)
        return 1
    if payload.get("midpoint_defect_emitter_closed") is not True or payload.get("delta_midpoint_zero_on_current_support") is not True:
        print("charged ordered-package value law should certify that the same-support midpoint defect already closes to zero", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "oph_charged_sector_local_current_support_obstruction_certificate":
        print("charged ordered-package value law should now point to the current-support obstruction certificate", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
