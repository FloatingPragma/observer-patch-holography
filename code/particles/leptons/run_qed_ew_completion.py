#!/usr/bin/env python3
"""Build the charged-lepton QED/EW completion artifact with promotion gating."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from math import sqrt


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "leptons" / "forward_charged_leptons.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "qed_ew_completion.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a declared QED/EW completion artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input charged-lepton forward artifact.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output completion JSON path.")
    parser.add_argument("--vev", type=float, default=246.22, help="Electroweak vev to use for the tree-level bridge.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    closure_state = str(payload.get("closure_state", payload.get("theorem_status", "open")))
    shape_values = [float(item) for item in payload.get("singular_values_shape", [])]
    shape_only_bridge = [args.vev / sqrt(2.0) * item for item in shape_values]

    if closure_state != "absolute_scale_closed":
        artifact = {
            "artifact": "oph_lepton_qed_ew_completion",
            "generated_utc": _timestamp(),
            "labels": payload.get("labels"),
            "source_artifact": payload.get("artifact"),
            "vev_used_gev": args.vev,
            "input_closure_state": closure_state,
            "shape_only_bridge_gev": shape_only_bridge,
            "reported_masses_gev": None,
            "promotion_blocked": True,
            "completion_kind": "blocked_pending_absolute_scale",
            "metadata": {
                "note": "QED/EW completion stays blocked until the charged-lepton lane reaches absolute-scale closure.",
            },
        }
    else:
        singular_values = [float(item) for item in payload.get("singular_values_abs", [])]
        tree_level_masses = [args.vev / sqrt(2.0) * item for item in singular_values]
        artifact = {
            "artifact": "oph_lepton_qed_ew_completion",
            "generated_utc": _timestamp(),
            "labels": payload.get("labels"),
            "source_artifact": payload.get("artifact"),
            "vev_used_gev": args.vev,
            "input_closure_state": closure_state,
            "tree_level_masses_gev": tree_level_masses,
            "reported_masses_gev": tree_level_masses,
            "promotion_blocked": False,
            "completion_kind": "declared_qed_ew_completion",
            "metadata": {
                "note": "QED/EW completion applied after charged-lepton absolute-scale closure.",
            },
        }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
