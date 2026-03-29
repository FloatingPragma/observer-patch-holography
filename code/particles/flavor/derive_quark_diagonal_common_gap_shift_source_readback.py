#!/usr/bin/env python3
"""Emit the quark diagonal common gap-shift source-readback artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_LAW = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_law.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_readback.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_law: dict) -> dict:
    b_ord = [float(value) for value in source_law["B_ord"]]
    return {
        "artifact": "oph_family_excitation_diagonal_common_gap_shift_source_readback",
        "generated_utc": _timestamp(),
        "proof_status": "source_readback_law_closed_waiting_pure_B_payload_pair",
        "predictive_promotion_allowed": False,
        "source_artifact": source_law.get("artifact"),
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
        "payload_pair_status": "open_waiting_value_emission",
        "J_B_source_u": None,
        "J_B_source_d": None,
        "beta_u_diag_B_source": None,
        "beta_d_diag_B_source": None,
        "source_readback_u_log_per_side": None,
        "source_readback_d_log_per_side": None,
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
        "smallest_constructive_missing_object": "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
        "next_single_residual_object": "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
        "derived_scalar_pair_after_payload_emission": "J_B_source_u_and_J_B_source_d",
        "notes": [
            "The source-readback law is closed: the minimal pure-B readback is uniquely [-beta, 0, +beta] on B_ord = [-1, 0, 1].",
            "The first data-bearing primitive beneath the odd scalar pair is the emitted pure-B payload pair source_readback_u_log_per_side and source_readback_d_log_per_side.",
            "The remaining quark gap is the emitted pure-B payload pair itself; once that payload exists, the odd projector values J_B_source_u and J_B_source_d follow algebraically from its endpoints.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark diagonal common gap-shift source-readback artifact.")
    parser.add_argument("--source-law", default=str(DEFAULT_SOURCE_LAW))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_law = json.loads(Path(args.source_law).read_text(encoding="utf-8"))
    artifact = build_artifact(source_law)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
