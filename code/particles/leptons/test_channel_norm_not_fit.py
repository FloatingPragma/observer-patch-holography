#!/usr/bin/env python3
"""Ensure unresolved channel norms do not silently become fitted absolutes."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_channel_norm.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the charged-lepton channel-norm artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input channel-norm artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    closed = bool(payload.get("channel_norm_closed", False))
    g_e = payload.get("g_e")
    proof_status = str(payload.get("proof_status", "open"))

    if not closed and g_e is not None:
        print("channel norm is open but g_e is populated", file=sys.stderr)
        return 1
    if closed and g_e is None:
        print("channel norm is marked closed but g_e is missing", file=sys.stderr)
        return 1
    if closed and proof_status not in {"sector_local_closed", "shared_budget_closed"}:
        print("channel norm is marked closed without closed proof status", file=sys.stderr)
        return 1
    if proof_status == "shared_budget_closed" and payload.get("closure_route") != "shared_charged_budget":
        print("shared-budget closure is missing its explicit closure route", file=sys.stderr)
        return 1

    print("channel-norm closure state is consistent")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
