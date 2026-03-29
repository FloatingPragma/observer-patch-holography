#!/usr/bin/env python3
"""Block exact-scale promotion of the Majorana action germ without the scalar evaluator."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_action_germ.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate action-germ scale gating.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if str(payload.get("scale_status", "")) == "closed" and payload.get("upstream_missing_object") == "oph_majorana_overlap_defect_scalar_evaluator":
        print("action germ claims closed scale while the scalar evaluator is still missing", file=sys.stderr)
        return 1
    print("Majorana action-germ scale gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
