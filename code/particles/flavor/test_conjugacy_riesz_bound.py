#!/usr/bin/env python3
"""Guard the explicit transported defect / spectral gap ratio for the flavor kernel."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the conjugacy-Riesz defect/gap bound.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input overlap-edge cocycle artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    gap = float(payload.get("theorem_gap_gamma", 0.0))
    ratio = payload.get("defect_gap_ratio")
    if gap <= 0.0:
        print("missing positive theorem gap gamma", file=sys.stderr)
        return 1
    if ratio is None:
        print("missing defect_gap_ratio in overlap-edge transport cocycle", file=sys.stderr)
        return 1
    if not bool(payload.get("riesz_bound_passes", False)):
        print("conjugacy-Riesz bound does not pass", file=sys.stderr)
        return 1
    print("conjugacy-Riesz bound guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
