#!/usr/bin/env python3
"""Validate the common-refinement projector-overlap certificate."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_line_lift.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate common-refinement projector-overlap invariance.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    cert = dict(payload.get("common_refinement_overlap_certificate", {}))
    if not bool(cert.get("common_refinement_invariance_closed_on_current_family", False)):
        print("common-refinement overlap invariance is not closed on the current family", file=sys.stderr)
        return 1
    if not bool(cert.get("all_edge_pairs_nondegenerate", False)):
        print("common-refinement overlap certificate has degenerate edge pairs", file=sys.stderr)
        return 1
    print("common-refinement projector-overlap invariance guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

