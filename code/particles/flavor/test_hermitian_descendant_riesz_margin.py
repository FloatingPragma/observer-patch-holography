#!/usr/bin/env python3
"""Require the lifted Hermitian-descendant Riesz margin to pass on the current family."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the Hermitian-descendant Riesz margin.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    margin = dict(payload.get("hermitian_descendant_riesz_margin", {}))
    if not margin:
        print("missing Hermitian-descendant Riesz margin", file=sys.stderr)
        return 1
    if not bool(margin.get("passes", False)):
        print("Hermitian-descendant Riesz margin does not pass", file=sys.stderr)
        return 1
    if float(margin.get("hermitian_descendant_norm_direct", 1.0)) >= float(margin.get("gamma_half", 0.0)):
        print("direct Hermitian-descendant bound is not below gamma/2", file=sys.stderr)
        return 1
    print("Hermitian-descendant Riesz margin guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
