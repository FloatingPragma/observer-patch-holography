#!/usr/bin/env python3
"""Validate residual normalization classes for sector response objects."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate residual normalization classes.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input sector-response JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    response_object = dict(payload.get("sector_response_object", {}))
    failures: list[str] = []

    for sector, item in response_object.items():
        norm_class = item.get("normalization_class")
        norms = dict(item.get("residual_norms", {}))
        if sector == "nu":
            if norm_class != "symmetric_diagonal":
                failures.append("nu: expected symmetric_diagonal normalization class")
            if "diag" not in norms:
                failures.append("nu: missing diagonal residual norms")
        elif sector == "e":
            if norm_class != "scalar_only":
                failures.append("e: expected scalar_only normalization class")
            if "mu" not in norms:
                failures.append("e: missing scalar residual norm")
        else:
            if norm_class != "left_right_diagonal":
                failures.append(f"{sector}: expected left_right_diagonal normalization class")
            if "left" not in norms or "right" not in norms:
                failures.append(f"{sector}: missing left/right residual norms")
        if sector in {"u", "d", "e"}:
            scalarization = dict(item.get("charged_dirac_scalarization_certificate", {}))
            if scalarization.get("functional_kind") != "charged_dirac_scalarization_candidate":
                failures.append(f"{sector}: missing charged-dirac scalarization candidate surface")

    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1

    print("sector residual factorization checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
