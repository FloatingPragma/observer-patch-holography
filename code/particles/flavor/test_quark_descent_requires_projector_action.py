#!/usr/bin/env python3
"""Require projector-resolved quark action before any certified quark promotion."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate projector-resolved quark descent.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input forward Yukawa artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if bool(payload.get("forward_certified", False)) and not bool(payload.get("uses_full_projector_algebra", False)):
        print("forward-certified quark artifact lacks projector-resolved descent", file=sys.stderr)
        return 1
    print("quark projector-action guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
