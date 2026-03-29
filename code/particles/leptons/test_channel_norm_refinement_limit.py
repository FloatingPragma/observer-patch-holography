#!/usr/bin/env python3
"""Ensure a closed channel norm carries explicit refinement-limit evidence."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_channel_norm.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate refinement-limit evidence for g_e closure.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input lepton channel-norm artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if payload.get("proof_status") != "sector_local_closed":
        print("channel norm not closed; refinement-limit guard idle")
        return 0

    stream = list(payload.get("g_e_by_refinement", []))
    certificate = dict(payload.get("channel_norm_refinement_certificate", {}))
    if len(stream) < 2:
        print("sector-local closure is missing a refinement stream", file=sys.stderr)
        return 1
    if certificate.get("status") == "snapshot_only":
        print("sector-local closure cannot rest on snapshot-only evidence", file=sys.stderr)
        return 1

    print("closed g_e carries refinement-limit evidence")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
