#!/usr/bin/env python3
"""Validate basic functoriality and residual-burden rules for sector pushforward."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the sector response object.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input sector-response JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if payload.get("sector_pushforward_kind") != "family_observable_to_sector_response":
        print("unexpected sector_pushforward_kind", file=sys.stderr)
        return 1
    response_object = dict(payload.get("sector_response_object", {}))
    failures: list[str] = []

    for sector, item in response_object.items():
        functoriality = dict(item.get("functoriality_certificate", {}))
        if not bool(functoriality.get("conjugation_class_stable", False)):
            failures.append(f"{sector}: conjugation_class_stable is false")

        if sector != "nu" and not bool(item.get("orientation_preserved", False)):
            failures.append(f"{sector}: orientation should be preserved for Dirac sectors")
        if sector == "nu" and item.get("symmetry_type") != "majorana":
            failures.append("nu: expected majorana symmetry type")
        if "A" in item:
            failures.append(f"{sector}: free entrywise amplitude matrix leaked into sector response")
        if sector in {"u", "d", "e"}:
            scalarization = dict(item.get("charged_dirac_scalarization_certificate", {}))
            if not bool(scalarization):
                failures.append(f"{sector}: missing charged_dirac_scalarization_certificate")
            elif not bool(scalarization.get("conjugation_invariant", False)):
                failures.append(f"{sector}: charged_dirac_scalarization_certificate lost conjugation_invariant")

    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1

    print("sector pushforward functoriality checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
