#!/usr/bin/env python3
"""Fail if isotropic selector data are overpromoted beyond the equal-split theorem."""

from __future__ import annotations

import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
LIFT = ROOT / "particles" / "runs" / "neutrino" / "majorana_holonomy_lift.json"
PULLBACK = ROOT / "particles" / "runs" / "neutrino" / "majorana_phase_pullback_metric.json"


def main() -> int:
    lift = json.loads(LIFT.read_text(encoding="utf-8"))
    pullback = json.loads(PULLBACK.read_text(encoding="utf-8")) if PULLBACK.exists() else {}
    isotropic = bool((lift.get("edge_weight_isotropy_certificate") or {}).get("closed"))
    if not isotropic:
        print("selector isotropy gate skipped on nonisotropic data")
        return 0
    if lift.get("selector_equivalence_class") != "principal_equal_split":
        print("isotropic selector branch lost its principal_equal_split classification", file=sys.stderr)
        return 1
    if lift.get("selector_closure_reason") != "s3_fixed_point":
        print("isotropic selector branch lost its s3_fixed_point reason", file=sys.stderr)
        return 1
    if not pullback or not pullback.get("phase_action_closed", False):
        if lift.get("canonical_selector_status") != "closed_equal_split":
            print("isotropic selector branch overpromoted beyond closed_equal_split without a phase-action theorem", file=sys.stderr)
            return 1
        if lift.get("selector_law_status") != "candidate_only":
            print("selector law was promoted without a closed phase-action theorem", file=sys.stderr)
            return 1
    print("isotropic selector reason gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
