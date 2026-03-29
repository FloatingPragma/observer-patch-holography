#!/usr/bin/env python3
"""Ensure a shared-budget artifact does not silently localize into sector-local g_e."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_channel_norm.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate shared-budget handling in the lepton channel norm.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input lepton channel-norm artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if payload.get("proof_status") != "shared_budget_only":
        print("no shared-budget-only state present; nothing to guard")
        return 0

    if payload.get("g_e") is not None or payload.get("channel_norm_closed"):
        print("shared-budget-only artifact was silently localized into g_e", file=sys.stderr)
        return 1
    if not payload.get("shared_budget_key"):
        print("shared-budget-only artifact is missing its shared budget key", file=sys.stderr)
        return 1
    if payload.get("closure_route") not in {None, "shared_charged_budget"}:
        print("shared-budget-only artifact was mislabeled as a local closure route", file=sys.stderr)
        return 1

    print("shared-budget-only state stays explicitly non-local")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
