#!/usr/bin/env python3
"""Ensure flavor cycle holonomy is sourced from the overlap-edge cocycle."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_COCYCLE = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
DEFAULT_OBSERVABLE = ROOT / "particles" / "runs" / "flavor" / "flavor_observable_artifact.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate cycle holonomy provenance.")
    parser.add_argument("--cocycle", default=str(DEFAULT_COCYCLE), help="Input cocycle artifact.")
    parser.add_argument("--observable", default=str(DEFAULT_OBSERVABLE), help="Input flavor observable artifact.")
    args = parser.parse_args()

    cocycle = json.loads(pathlib.Path(args.cocycle).read_text(encoding="utf-8"))
    observable = json.loads(pathlib.Path(args.observable).read_text(encoding="utf-8"))
    derived = dict(cocycle.get("derived_cycle_holonomy", {}))
    cycle = dict(observable.get("cycle_phases", {}))
    if derived != cycle:
        print("observable cycle holonomy is not reconstructed from the overlap-edge cocycle", file=sys.stderr)
        return 1
    print("cycle holonomy provenance guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
