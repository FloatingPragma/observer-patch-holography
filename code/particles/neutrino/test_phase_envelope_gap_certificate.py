#!/usr/bin/env python3
"""Check that the residual phase envelope gates ordering promotion."""

from __future__ import annotations

import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
ENVELOPE = ROOT / "particles" / "runs" / "neutrino" / "majorana_phase_envelope.json"
SPLITTINGS = ROOT / "particles" / "runs" / "neutrino" / "forward_splittings.json"


def main() -> int:
    envelope = json.loads(ENVELOPE.read_text(encoding="utf-8"))
    splittings = json.loads(SPLITTINGS.read_text(encoding="utf-8"))
    certificate = envelope.get("gap_vs_radius_certificate") or {}
    if not certificate:
        print("missing gap_vs_radius_certificate", file=sys.stderr)
        return 1
    if str(splittings.get("ordering_theorem_status", "")).startswith("selector_"):
        if not splittings.get("ordering_phase_certified"):
            print("selector-certified ordering is missing the certified label", file=sys.stderr)
            return 1
        print("phase envelope gate bypassed legitimately by selector certification")
        return 0
    if envelope.get("ordering_phase_stable"):
        if not splittings.get("ordering_phase_certified"):
            print("ordering should be certified when the envelope says it is phase-stable", file=sys.stderr)
            return 1
    else:
        if splittings.get("ordering_phase_certified") is not None:
            print("ordering was promoted without a phase-stability certificate", file=sys.stderr)
            return 1
    if splittings.get("phase_certificate_source") != str(ENVELOPE):
        print("splittings are not pointing at the envelope artifact as the phase certificate source", file=sys.stderr)
        return 1
    print("phase envelope correctly gates ordering promotion")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
