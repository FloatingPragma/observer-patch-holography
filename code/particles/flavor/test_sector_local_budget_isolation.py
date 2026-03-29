#!/usr/bin/env python3
"""Block sector-local budget claims unless an isolation witness closes."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate sector-local budget isolation claims.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input charged-budget artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    witnesses = dict(payload.get("sector_isolation_witness_by_sector", {}))
    for sector, witness in witnesses.items():
        if not bool(dict(witness).get("sector_local_closed", False)):
            continue
        if dict(witness).get("value") not in {0, 0.0}:
            print(f"sector {sector} claims sector-local closure without zero isolation witness", file=sys.stderr)
            return 1

    print("sector-local budget claims remain gated by an isolation witness")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
