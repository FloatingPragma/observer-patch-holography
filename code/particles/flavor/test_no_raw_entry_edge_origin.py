#!/usr/bin/env python3
"""Ensure the flavor cocycle no longer presents itself as raw matrix-entry readback."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate non-readback cocycle origin.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if payload.get("cocycle_origin_status") != "induced_from_projective_eigenline_transport":
        print("edge cocycle origin still looks like a primitive readback object rather than an induced edge object", file=sys.stderr)
        return 1
    if payload.get("upstream_missing_object") != "oph_generation_bundle_branch_generator_splitting":
        print("edge cocycle has not inherited the reduced generation-branch splitting blocker", file=sys.stderr)
        return 1
    print("no raw-entry edge-origin guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
