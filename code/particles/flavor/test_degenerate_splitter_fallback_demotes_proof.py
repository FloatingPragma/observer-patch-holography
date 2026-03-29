#!/usr/bin/env python3
"""Require degenerate odd-splitter fallback to demote quark descent rather than silently promote it."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_descent.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate degenerate odd-splitter fallback handling.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if bool(payload.get("degenerate_placeholder_fallback_used", False)) and str(payload.get("quark_descent_proof_status", "")) != "open":
        print("degenerate splitter fallback did not demote quark descent to open", file=sys.stderr)
        return 1
    print("degenerate splitter fallback guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
