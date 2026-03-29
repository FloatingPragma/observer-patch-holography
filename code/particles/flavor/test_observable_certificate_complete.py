#!/usr/bin/env python3
"""Ensure the flavor observable exports the full theorem-scope certificate surface."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "flavor_observable_artifact.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate flavor observable certificate completeness.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input flavor observable artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    certificate = dict(payload.get("persistent_projector_certificate", {}))
    required = ("conjugacy_defect_sup", "theorem_gap_gamma", "defect_gap_ratio")
    missing = [key for key in required if certificate.get(key) is None]
    if missing:
        print(f"flavor observable certificate is missing: {', '.join(missing)}", file=sys.stderr)
        return 1
    print("observable certificate completeness guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
