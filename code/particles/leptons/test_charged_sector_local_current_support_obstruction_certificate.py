#!/usr/bin/env python3
"""Smoke-test the charged current-support obstruction certificate artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_current_support_obstruction_certificate.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_current_support_obstruction_certificate.json"


def main() -> int:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_charged_sector_local_current_support_obstruction_certificate":
        print("wrong charged obstruction certificate artifact id", file=sys.stderr)
        return 1
    if payload.get("same_support_exhausted") is not True:
        print("charged obstruction certificate should certify same-support exhaustion", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "oph_charged_sector_local_minimal_source_support_extension_emitter":
        print("charged obstruction certificate should advance to the minimal support extension emitter", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
