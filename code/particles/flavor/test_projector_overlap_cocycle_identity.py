#!/usr/bin/env python3
"""Validate the common-refinement projector-overlap cocycle identity surface."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate projector-overlap cocycle identity.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if not bool(payload.get("cocycle_identity_closed", False)):
        print("projector-overlap cocycle identity is not closed on the current family", file=sys.stderr)
        return 1
    if not bool(payload.get("refinement_functoriality_closed", False)):
        print("projector-overlap cocycle identity is missing refinement functoriality closure", file=sys.stderr)
        return 1
    print("projector-overlap cocycle identity guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
