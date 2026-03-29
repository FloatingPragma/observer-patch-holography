#!/usr/bin/env python3
"""Ensure the quark sector-descent artifact really separates u and d."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_descent.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate quark sector distinctness.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input quark-sector-descent artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    witness = dict(payload.get("sector_distinctness_witness", {}))
    if witness.get("status") != "closed" or float(witness.get("value", 0.0)) <= 1.0e-12:
        print("quark sector descent still leaves u and d cloned", file=sys.stderr)
        return 1
    print("quark sector nonclone guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
