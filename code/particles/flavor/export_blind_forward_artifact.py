#!/usr/bin/env python3
"""Strip compare surfaces from a flavor artifact and export a blind forward JSON."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "blind_forward_artifact.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Export a blind forward flavor artifact.")
    parser.add_argument("--input", required=True, help="Input JSON path from build_forward_yukawas.py.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    blind = {
        "artifact": "oph_flavor_blind_forward",
        "generated_utc": _timestamp(),
        "source_mode": payload.get("source_mode"),
        "certification_status": payload.get("certification_status"),
        "promotion_blockers": payload.get("promotion_blockers"),
        "template_amplitude_fallback_used": payload.get("template_amplitude_fallback_used"),
        "up_down_sector_distinct": payload.get("up_down_sector_distinct"),
        "Y_u": payload.get("Y_u"),
        "Y_d": payload.get("Y_d"),
        "singular_values_u": payload.get("singular_values_u"),
        "singular_values_d": payload.get("singular_values_d"),
        "U_u_left": payload.get("U_u_left"),
        "U_d_left": payload.get("U_d_left"),
        "V_CKM": payload.get("V_CKM"),
        "jarlskog": payload.get("jarlskog"),
        "commutator_invariant": payload.get("commutator_invariant"),
        "metadata": payload.get("metadata", {}),
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(blind, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
