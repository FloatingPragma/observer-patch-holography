#!/usr/bin/env python3
"""Require the overlap-edge line-lift artifact to stay explicit about its open status."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_line_lift.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the overlap-edge line-lift boundary.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if payload.get("proof_status") == "closed" and payload.get("upstream_missing_object") == "missing_overlap_edge_line_lift_theorem":
        print("line-lift artifact claims closure while still naming the line-lift theorem as missing", file=sys.stderr)
        return 1
    if not bool(payload.get("raw_entry_readback_forbidden_as_closed_origin", False)):
        print("line-lift artifact does not explicitly forbid raw-entry promotion", file=sys.stderr)
        return 1
    print("overlap-edge line-lift boundary guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
