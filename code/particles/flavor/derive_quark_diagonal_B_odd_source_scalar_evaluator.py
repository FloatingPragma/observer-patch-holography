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
    j_b_u = source_readback.get("J_B_source_u")
    j_b_d = source_readback.get("J_B_source_d")
    values_present = j_b_u is not None and j_b_d is not None
    payload_closed = (
        values_present
        and source_readback.get("payload_pair_status")
        == "emitted_selected_public_class"
    )
    comparison_only = values_present and not payload_closed
    return {
        "artifact": "oph_quark_diagonal_B_odd_source_scalar_evaluator",
        "generated_utc": _timestamp(),
        "proof_status": (
            "closed_public_selected_class_pure_B_odd_source_scalar_evaluator"
            if payload_closed
            else "comparison_only_odd_scalar_projection_value_source_open"
            if comparison_only
            else "pure_B_odd_projection_formula_closed_source_values_open"
        ),
        "predictive_promotion_allowed": False,
        "source_artifact": source_readback.get("artifact"),
        "source_readback_status": source_readback.get("proof_status"),
        "source_payload_artifact": source_readback.get("source_payload_artifact"),
        "source_payload_status": source_readback.get("source_payload_status"),
        "B_ord": b_ord,
        "B_ord_norm_sq": b_norm_sq,
        "J_B_functional_kind": "pure_B_odd_point_separating_projection",
        "J_B_formula": "dot(v, B_ord) / B_ord_norm_sq",
        "J_B_endpoint_formula": "(v[2] - v[0]) / 2",
        "J_B_on_B_ord": 1.0,
        "J_B_on_center_vector": 0.0,
        "J_B_on_Q_ord": 0.0,
        "J_B_source_u": j_b_u,
        "J_B_source_d": j_b_d,
        "beta_u_diag_B_source_identification": "J_B_source_u",
        "beta_d_diag_B_source_identification": "J_B_source_d",
        "beta_u_diag_B_source": j_b_u,
        "beta_d_diag_B_source": j_b_d,
        "beta_u_diag_B_source_formula": "J_B_source_u",
        "beta_d_diag_B_source_formula": "J_B_source_d",
        "source_readback_u_log_per_side": source_readback.get("source_readback_u_log_per_side"),
        "source_readback_d_log_per_side": source_readback.get("source_readback_d_log_per_side"),
        "predictive_J_B_source_law_status": (
            "selected_public_class_closed"
            if payload_closed
            else "comparison_only_unpromotable"
            if comparison_only
            else "missing"
        ),
        "smallest_constructive_missing_object": (
            "off_canonical_pure_B_source_payload_family"
            if payload_closed
            else "source_derived_c_d_over_c_u_or_t1_value"
            if comparison_only
            else "source_readback_u_log_per_side_and_source_readback_d_log_per_side"
        ),
        "next_single_residual_object": (
            "off_canonical_pure_B_source_payload_family"
            if payload_closed
            else "source_derived_c_d_over_c_u_or_t1_value"
            if comparison_only
            else "source_readback_u_log_per_side_and_source_readback_d_log_per_side"
        ),
        "derived_scalar_pair_after_payload_emission": "J_B_source_u_and_J_B_source_d",
        "notes": [
            "This artifact isolates the odd projector that will read back the emitted pure-B payload pair.",
            "Values inherited from an unpromotable payload remain comparison-only; a source-derived c_d/c_u or t1 value is still required.",
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
