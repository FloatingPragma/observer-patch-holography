#!/usr/bin/env python3
"""Export a blind charged-lepton artifact without compare surfaces."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "leptons" / "forward_charged_leptons.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "blind_forward_artifact.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Export a blind charged-lepton artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input JSON path.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    blind = {
        "artifact": "oph_charged_lepton_blind_forward",
        "generated_utc": _timestamp(),
        "labels": payload.get("labels"),
        "closure_state": payload.get("closure_state"),
        "Y_e_shape": payload.get("Y_e_shape"),
        "singular_values_shape": payload.get("singular_values_shape"),
        "Y_e": payload.get("Y_e"),
        "singular_values_abs": payload.get("singular_values_abs"),
        "U_e_left": payload.get("U_e_left"),
        "U_e_right": payload.get("U_e_right"),
        "metadata": payload.get("metadata", {}),
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(blind, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
