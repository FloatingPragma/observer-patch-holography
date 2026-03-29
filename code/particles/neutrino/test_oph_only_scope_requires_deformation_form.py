#!/usr/bin/env python3
"""Block OPH-only selector-law promotion without a closed deformation form."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_PULLBACK = ROOT / "particles" / "runs" / "neutrino" / "majorana_phase_pullback_metric.json"
DEFAULT_DEFORMATION = ROOT / "particles" / "runs" / "neutrino" / "majorana_deformation_bilinear_form.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate OPH-only ambient-metric promotion rules.")
    parser.add_argument("--pullback", default=str(DEFAULT_PULLBACK), help="Input pullback-metric artifact.")
    parser.add_argument("--deformation", default=str(DEFAULT_DEFORMATION), help="Input deformation-form artifact.")
    args = parser.parse_args()

    pullback = json.loads(pathlib.Path(args.pullback).read_text(encoding="utf-8"))
    deformation = json.loads(pathlib.Path(args.deformation).read_text(encoding="utf-8")) if pathlib.Path(args.deformation).exists() else {}

    law_scope = str(pullback.get("law_closure_scope", ""))
    if law_scope == "oph_only" and str(deformation.get("oph_origin_status", "")) != "closed":
        print("oph_only selector-law closure claimed without a closed deformation bilinear form", file=sys.stderr)
        return 1

    print("oph-only ambient-metric guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
