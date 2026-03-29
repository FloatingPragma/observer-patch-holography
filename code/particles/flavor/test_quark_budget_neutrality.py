#!/usr/bin/env python3
"""Ensure the odd quark split stays budget-neutral inside the shared charged branch."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_descent.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate quark budget neutrality.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input quark-sector-descent artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    certificate = dict(payload.get("budget_neutrality_certificate", {}))
    if certificate.get("status") != "closed":
        print("quark budget neutrality is not certified", file=sys.stderr)
        return 1
    if abs(float(certificate.get("u_plus_d_odd_part_sum", 1.0))) > 1.0e-12:
        print("quark odd split is not budget neutral", file=sys.stderr)
        return 1
    print("quark budget-neutrality guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
