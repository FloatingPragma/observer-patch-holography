#!/usr/bin/env python3
"""Require a noncentrality witness for quark CKM/CP readiness claims."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_descent.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the quark noncentrality witness.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input quark-sector-descent artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    witness = dict(payload.get("noncentrality_witness", {}))
    if witness.get("status") != "closed" or float(witness.get("fro_norm", 0.0)) <= 1.0e-18:
        print("missing noncentrality witness for the quark odd split", file=sys.stderr)
        return 1
    print("quark noncentrality witness guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
