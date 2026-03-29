#!/usr/bin/env python3
"""Fail if the flavor edge statistics remain clamp-saturated placeholders."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate non-placeholder edge statistics.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input overlap-edge cocycle artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    certificate = dict(payload.get("non_floor_amplitude_certificate", {}))
    if certificate.get("status") != "closed":
        print("overlap-edge cocycle still uses floor-saturated amplitudes", file=sys.stderr)
        return 1
    if bool(certificate.get("all_equal_off_diagonal", True)):
        print("overlap-edge cocycle still has all-equal off-diagonal amplitudes", file=sys.stderr)
        return 1
    print("edge-statistics non-placeholder guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
