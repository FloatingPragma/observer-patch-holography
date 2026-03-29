#!/usr/bin/env python3
"""Fail if the neutrino lane silently promotes the real seed as a resolved phase output."""

from __future__ import annotations

import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
LIFT = ROOT / "particles" / "runs" / "neutrino" / "majorana_holonomy_lift.json"
MAJORANA = ROOT / "particles" / "runs" / "neutrino" / "forward_majorana_matrix.json"


def main() -> int:
    lift = json.loads(LIFT.read_text(encoding="utf-8"))
    majorana = json.loads(MAJORANA.read_text(encoding="utf-8"))
    omega = float(lift["cycle_constraint"]["omega_012"])
    if abs(omega) > 1.0e-15 and majorana.get("phase_mode") == "real_seed":
        if majorana.get("certification_status") not in {"real_seed_phase_unresolved", "real_seed_surrogate"}:
            print("real-seed mode was promoted despite a nonzero Majorana phase obstruction", file=sys.stderr)
            return 1
    print("no zero-phase shortcut in neutrino lane")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
