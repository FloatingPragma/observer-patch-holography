#!/usr/bin/env python3
"""Emit the quark diagonal common gap-shift source-readback artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_LAW = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_law.json"
DEFAULT_SOURCE_PAYLOAD = ROOT / "particles" / "runs" / "flavor" / "quark_d12_public_source_payload.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_readback.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_law: dict, source_payload: dict | None = None) -> dict:
    b_ord = [float(value) for value in source_law["B_ord"]]
    source_payload = source_payload or {}
    source_u = source_payload.get("source_readback_u_log_per_side")
    source_d = source_payload.get("source_readback_d_log_per_side")
    beta_u = source_payload.get("beta_u_diag_B_source")
    beta_d = source_payload.get("beta_d_diag_B_source")
    j_b_u = source_payload.get("J_B_source_u")
    j_b_d = source_payload.get("J_B_source_d")
    values_present = (
        source_u is not None
        and source_d is not None
        and beta_u is not None
        and beta_d is not None
    )
    payload_closed = (
        values_present and source_payload.get("public_promotion_allowed") is True
    )
    comparison_only = values_present and not payload_closed
    return {
        "artifact": "oph_family_excitation_diagonal_common_gap_shift_source_readback",
        "generated_utc": _timestamp(),
        "proof_status": (
            "closed_public_selected_class_source_readback"
            if payload_closed
            else "comparison_only_readback_from_unpromotable_payload"
            if comparison_only
            else "source_readback_law_closed_waiting_pure_B_payload_pair"
        ),
        "predictive_promotion_allowed": False,
        "source_artifact": source_law.get("artifact"),
        "source_payload_artifact": source_payload.get("artifact"),
        "source_payload_status": source_payload.get("proof_status"),
        "source_payload_promotion_allowed": source_payload.get(
            "public_promotion_allowed"
        ),
        "value_classification": (
            "source_emitted" if payload_closed else "comparison_only" if comparison_only else "absent"
        ),
        "B_ord": b_ord,
        "B_ord_norm_sq": sum(value * value for value in b_ord),
        "J_B_functional_kind": "pure_B_odd_point_separating_projection",
        "J_B_formula": "dot(v, B_ord) / B_ord_norm_sq",
        "J_B_endpoint_formula": "(v[2] - v[0]) / 2",
        "J_B_on_B_ord": 1.0,
        "J_B_on_center_vector": 0.0,
        "J_B_on_Q_ord": 0.0,
        "source_payload_kind": "pure_B_odd_source_readback_pair",
        "source_payload_shape": [3],
        "source_payload_constraints": {
            "center_entry": 0.0,
            "endpoint_sum": 0.0,
            "pure_B_residual": [0.0, 0.0, 0.0],
        },
        "first_data_bearing_primitive_beneath_scalar_pair": "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
        "payload_pair_status": (
            "emitted_selected_public_class"
            if payload_closed
            else "comparison_only_unpromotable"
            if comparison_only
            else "open_waiting_value_emission"
        ),
        "J_B_source_u": j_b_u,
        "J_B_source_d": j_b_d,
        "beta_u_diag_B_source": beta_u,
        "beta_d_diag_B_source": beta_d,
        "source_readback_u_log_per_side": source_u,
        "source_readback_d_log_per_side": source_d,
        "J_B_source_u_formula": "J_B(source_readback_u_log_per_side)",
        "J_B_source_d_formula": "J_B(source_readback_d_log_per_side)",
        "beta_u_diag_B_source_formula": "(source_readback_u_log_per_side[2] - source_readback_u_log_per_side[0]) / 2",
        "beta_d_diag_B_source_formula": "(source_readback_d_log_per_side[2] - source_readback_d_log_per_side[0]) / 2",
        "source_readback_u_log_per_side_formula": "beta_u_diag_B_source * B_ord",
        "source_readback_d_log_per_side_formula": "beta_d_diag_B_source * B_ord",
        "source_readback_center_residual_u": 0.0,
        "source_readback_center_residual_d": 0.0,
        "source_readback_endpoint_residual_u": 0.0,
        "source_readback_endpoint_residual_d": 0.0,
        "B_mode_residual_u": [0.0, 0.0, 0.0],
        "B_mode_residual_d": [0.0, 0.0, 0.0],
        "pure_B_source_certificate_u": 0.0,
        "pure_B_source_certificate_d": 0.0,
        "source_readback_center_residual_u_formula": "0.0",
        "source_readback_center_residual_d_formula": "0.0",
        "source_readback_endpoint_residual_u_formula": "0.0",
        "source_readback_endpoint_residual_d_formula": "0.0",
        "B_mode_residual_u_formula": "[0.0, 0.0, 0.0]",
        "B_mode_residual_d_formula": "[0.0, 0.0, 0.0]",
        "pure_B_source_certificate_u_formula": "0.0",
        "pure_B_source_certificate_d_formula": "0.0",
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
            "The source-readback law is closed: the minimal pure-B readback is uniquely [-beta, 0, +beta] on B_ord = [-1, 0, 1].",
            "Numerical values inherited from an unpromotable payload remain comparison-only even though the readback arithmetic is exact.",
            "The remaining source gap is a source-derived c_d/c_u or t1 value; once emitted, the odd projector values follow algebraically.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark diagonal common gap-shift source-readback artifact.")
    parser.add_argument("--source-law", default=str(DEFAULT_SOURCE_LAW))
    parser.add_argument("--source-payload", default=str(DEFAULT_SOURCE_PAYLOAD))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_law = json.loads(Path(args.source_law).read_text(encoding="utf-8"))
    source_payload_path = Path(args.source_payload)
    source_payload = json.loads(source_payload_path.read_text(encoding="utf-8")) if source_payload_path.exists() else None
    artifact = build_artifact(source_law, source_payload)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
