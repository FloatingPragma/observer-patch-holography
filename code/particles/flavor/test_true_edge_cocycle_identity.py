#!/usr/bin/env python3
"""Guard against silently promoting the readback cocycle into a true closed 1-cocycle."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate edge cocycle identity closure flags.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    closed = bool(payload.get("cocycle_identity_closed", False))
    refinement_closed = bool(payload.get("refinement_functoriality_closed", False))
    gauge_class = str(payload.get("vertex_rephasing_gauge_class", ""))
    descendant = payload.get("descendant_transport_operator_by_refinement")

    if closed and not refinement_closed:
        print("cocycle identity marked closed without refinement functoriality closure", file=sys.stderr)
        return 1
    if closed and not gauge_class:
        print("cocycle identity marked closed without a vertex rephasing gauge class", file=sys.stderr)
        return 1
    if descendant is None:
        print("missing descendant transport operator map by refinement", file=sys.stderr)
        return 1
    if not closed and str(payload.get("cocycle_origin_status", "")) == "closed":
        print("readback cocycle claims closed origin while identity remains open", file=sys.stderr)
        return 1
    print("true edge cocycle identity guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
