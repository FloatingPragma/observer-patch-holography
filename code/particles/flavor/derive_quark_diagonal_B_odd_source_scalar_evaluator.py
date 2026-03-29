#!/usr/bin/env python3
"""Isolate the remaining odd B-mode scalars in the quark diagonal lane.

Chain role: reduce the diagonal quark residual to the smallest source-side
scalars that still move the masses on the active branch.

Mathematics: projection onto the pure-B ordered basis and scalar readback on
that one-dimensional odd direction.

OPH-derived inputs: the ordered B-mode basis emitted by the quark source
readback artifact.

Output: the `J_B_source_u` / `J_B_source_d` evaluator shell consumed by the
quark exactness audit and completion prompts.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_READBACK = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_readback.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_B_odd_source_scalar_evaluator.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_readback: dict) -> dict:
    b_ord = [float(value) for value in source_readback["B_ord"]]
    b_norm_sq = float(source_readback["B_ord_norm_sq"])
    return {
        "artifact": "oph_quark_diagonal_B_odd_source_scalar_evaluator",
        "generated_utc": _timestamp(),
        "proof_status": "pure_B_odd_projection_formula_closed_source_values_open",
        "predictive_promotion_allowed": False,
        "source_artifact": source_readback.get("artifact"),
        "B_ord": b_ord,
        "B_ord_norm_sq": b_norm_sq,
        "J_B_functional_kind": "pure_B_odd_point_separating_projection",
        "J_B_formula": "dot(v, B_ord) / B_ord_norm_sq",
        "J_B_endpoint_formula": "(v[2] - v[0]) / 2",
        "J_B_on_B_ord": 1.0,
        "J_B_on_center_vector": 0.0,
        "J_B_on_Q_ord": 0.0,
        "J_B_source_u": None,
        "J_B_source_d": None,
        "beta_u_diag_B_source_identification": "J_B_source_u",
        "beta_d_diag_B_source_identification": "J_B_source_d",
        "beta_u_diag_B_source_formula": "J_B_source_u",
        "beta_d_diag_B_source_formula": "J_B_source_d",
        "predictive_J_B_source_law_status": "missing",
        "smallest_constructive_missing_object": "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
        "next_single_residual_object": "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
        "derived_scalar_pair_after_payload_emission": "J_B_source_u_and_J_B_source_d",
        "notes": [
            "This artifact isolates the odd projector that will read back the emitted pure-B payload pair.",
            "Once the payload pair source_readback_u_log_per_side and source_readback_d_log_per_side is emitted, J_B_source_u and J_B_source_d are fixed algebraically from its endpoints.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark diagonal B-odd source scalar evaluator.")
    parser.add_argument("--source-readback", default=str(DEFAULT_SOURCE_READBACK))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_readback = json.loads(Path(args.source_readback).read_text(encoding="utf-8"))
    artifact = build_artifact(source_readback)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
