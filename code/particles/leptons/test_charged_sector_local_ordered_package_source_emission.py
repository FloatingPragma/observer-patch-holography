#!/usr/bin/env python3
"""Smoke-test the charged source ordered-package artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_ordered_package_source_emission.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_source_emission.json"


def main() -> int:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_charged_sector_local_ordered_package_source_emission":
        print("wrong charged source ordered-package artifact id", file=sys.stderr)
        return 1
    if payload.get("proof_status") != "current_family_source_emission_closed":
        print("charged source ordered-package artifact should be closed on the current family", file=sys.stderr)
        return 1
    ordered = payload.get("source_side_ordered_package_log_per_side_emitted")
    if ordered != sorted(ordered):
        print("charged source ordered-package artifact should expose sorted package values", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
