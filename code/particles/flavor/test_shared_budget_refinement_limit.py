#!/usr/bin/env python3
"""Fail if shared-budget closure is claimed without a real refinement-limit certificate."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate shared-budget refinement-limit closure.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input charged-budget artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    proof_status = str(payload.get("proof_status", ""))
    certificate = dict(payload.get("refinement_limit_certificate", {}))

    if proof_status == "shared_budget_closed":
        if not bool(certificate.get("refinement_stable", False)):
            print("shared_budget_closed claimed without refinement_stable certificate", file=sys.stderr)
            return 1
        if int(certificate.get("samples", 0)) < 2:
            print("shared_budget_closed claimed without at least two refinement samples", file=sys.stderr)
            return 1

    print("shared-budget refinement-limit guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
