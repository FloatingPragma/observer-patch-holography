#!/usr/bin/env python3
"""Validate the explicit projective eigenline transport candidate surface."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_line_lift.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate projective eigenline transport candidate fields.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if payload.get("functor_kind") != "projective_polar_riesz_common_refinement":
        print("line-lift artifact is missing the explicit projective transport candidate kind", file=sys.stderr)
        return 1
    if payload.get("transport_group") != "objectwise_u1":
        print("line-lift artifact is missing the objectwise U(1) quotient structure", file=sys.stderr)
        return 1
    if payload.get("line_lift_is_readout_of") != "oph_intrinsic_generation_bundle_branch_generator":
        print("line-lift artifact is not explicitly marked as a readout of the intrinsic generation-branch generator", file=sys.stderr)
        return 1
    overlaps = list(payload.get("same_label_overlap_by_label_and_refinement_pair", []))
    if len(overlaps) != 3:
        print("line-lift artifact does not expose the three same-label diagonal overlaps", file=sys.stderr)
        return 1
    transports = list(payload.get("transport_partial_isometry_by_label_and_refinement_pair", []))
    if len(transports) != 3:
        print("line-lift artifact does not expose the three same-label transport maps", file=sys.stderr)
        return 1
    print("projective eigenline transport candidate guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
