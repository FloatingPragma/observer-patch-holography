#!/usr/bin/env python3
"""Ensure selector-point certification and selector-law certification stay distinct."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_MAJORANA = ROOT / "particles" / "runs" / "neutrino" / "forward_majorana_matrix.json"
DEFAULT_SPLITTINGS = ROOT / "particles" / "runs" / "neutrino" / "forward_splittings.json"
DEFAULT_PULLBACK = ROOT / "particles" / "runs" / "neutrino" / "majorana_phase_pullback_metric.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate selector-point versus selector-law status split.")
    parser.add_argument("--majorana", default=str(DEFAULT_MAJORANA))
    parser.add_argument("--splittings", default=str(DEFAULT_SPLITTINGS))
    parser.add_argument("--pullback", default=str(DEFAULT_PULLBACK))
    args = parser.parse_args()

    majorana = json.loads(pathlib.Path(args.majorana).read_text(encoding="utf-8"))
    splittings = json.loads(pathlib.Path(args.splittings).read_text(encoding="utf-8"))
    pullback = json.loads(pathlib.Path(args.pullback).read_text(encoding="utf-8"))

    if bool(majorana.get("selector_law_certified", False)) != bool(splittings.get("selector_law_certified", False)):
        print("majorana/splittings selector_law_certified mismatch", file=sys.stderr)
        return 1
    if bool(majorana.get("selector_point_certified", False)) != bool(splittings.get("selector_point_certified", False)):
        print("majorana/splittings selector_point_certified mismatch", file=sys.stderr)
        return 1
    if bool(majorana.get("selector_law_certified", False)) and not bool(pullback.get("phase_action_closed", False)):
        print("selector_law_certified claimed without closed pullback action", file=sys.stderr)
        return 1

    print("selector point vs law status split passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
