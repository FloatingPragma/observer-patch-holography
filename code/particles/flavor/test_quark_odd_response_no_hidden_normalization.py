#!/usr/bin/env python3
"""Reject hidden normalization when the odd quark response law claims coefficient-free closure."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "quark_odd_response_law.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate hidden-normalization gating for the quark odd response law.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    coefficient_free = bool(payload.get("coefficient_free", False))
    lift_constants = payload.get("lift_constants")
    if coefficient_free and lift_constants not in (None, {}):
        print("coefficient-free odd response law still carries lift constants", file=sys.stderr)
        return 1
    print("quark odd-response hidden-normalization guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
