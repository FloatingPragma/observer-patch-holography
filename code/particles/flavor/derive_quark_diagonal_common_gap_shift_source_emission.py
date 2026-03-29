#!/usr/bin/env python3
"""Emit the quark diagonal common gap-shift source-emission artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_LAW = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_law.json"
DEFAULT_SOURCE_READBACK = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_readback.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_emission.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_law: dict, source_readback: dict | None = None) -> dict:
    b_ord = [float(value) for value in source_law["B_ord"]]
    x2 = float(source_law["normalized_coordinate_x2"])
    q_ord = [
        (1.0 - x2 * x2) / 3.0,
        -(2.0 * (1.0 - x2 * x2)) / 3.0,
        (1.0 - x2 * x2) / 3.0,
    ]
    source_readback = source_readback or {}
    source_u = source_readback.get("source_readback_u_log_per_side")
    source_d = source_readback.get("source_readback_d_log_per_side")
    beta_u = None if source_u is None else sum(float(u) * b for u, b in zip(source_u, b_ord)) / sum(value * value for value in b_ord)
    beta_d = None if source_d is None else sum(float(d) * b for d, b in zip(source_d, b_ord)) / sum(value * value for value in b_ord)
    residual_u = None if source_u is None or beta_u is None else [float(u) - beta_u * b for u, b in zip(source_u, b_ord)]
    residual_d = None if source_d is None or beta_d is None else [float(d) - beta_d * b for d, b in zip(source_d, b_ord)]
    return {
        "artifact": "oph_family_excitation_diagonal_common_gap_shift_source_emission",
        "generated_utc": _timestamp(),
        "proof_status": (
            "source_emission_derived_from_source_readback"
            if source_u is not None and source_d is not None
            else "source_emission_waiting_pure_B_payload_pair"
        ),
        "predictive_promotion_allowed": False,
        "source_artifact": source_law.get("artifact"),
        "source_readback_artifact": source_readback.get("artifact"),
        "source_readback_status": source_readback.get("proof_status"),
        "B_ord": b_ord,
        "B_ord_norm_sq": sum(value * value for value in b_ord),
        "Q_ord": q_ord,
        "source_readback_u_log_per_side": source_u,
        "source_readback_d_log_per_side": source_d,
        "beta_u_diag_B_source": beta_u,
        "beta_d_diag_B_source": beta_d,
        "B_mode_residual_u": residual_u,
        "B_mode_residual_d": residual_d,
        "source_emission_status": (
            "source_emission_derived_from_source_readback"
            if source_u is not None and source_d is not None
            else "source_readback_law_closed_waiting_pure_B_payload_pair"
        ),
        "J_B_source_u": source_readback.get("J_B_source_u"),
        "J_B_source_d": source_readback.get("J_B_source_d"),
        "J_B_source_u_formula": source_readback.get("J_B_source_u_formula"),
        "J_B_source_d_formula": source_readback.get("J_B_source_d_formula"),
        "beta_u_diag_B_source_formula": "dot(source_readback_u_log_per_side, B_ord) / B_ord_norm_sq",
        "beta_d_diag_B_source_formula": "dot(source_readback_d_log_per_side, B_ord) / B_ord_norm_sq",
        "B_mode_residual_u_formula": "[0.0, 0.0, 0.0]",
        "B_mode_residual_d_formula": "[0.0, 0.0, 0.0]",
        "source_emission_closed_formula": "B_mode_residual_u == 0 and B_mode_residual_d == 0",
        "non_inverse_fit_certificate": {
            "status": "required",
            "reference_values_consumed": False,
        },
        "pdg_independence_certificate": {
            "status": "required",
            "reference_values_consumed": False,
        },
        "smallest_constructive_missing_object": source_readback.get(
            "smallest_constructive_missing_object",
            "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
        ),
        "notes": [
            "This artifact is now a derived projection layer on top of the closed pure-B source-readback law.",
            "The current mean and quadratic annihilators are blind to the B-mode amplitude, so the predictive gap sits first in the emitted pure-B payload pair beneath this projection.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark diagonal common gap-shift source-emission artifact.")
    parser.add_argument("--source-law", default=str(DEFAULT_SOURCE_LAW))
    parser.add_argument("--source-readback", default=str(DEFAULT_SOURCE_READBACK))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_law = json.loads(Path(args.source_law).read_text(encoding="utf-8"))
    source_readback_path = Path(args.source_readback)
    source_readback = json.loads(source_readback_path.read_text(encoding="utf-8")) if source_readback_path.exists() else None
    artifact = build_artifact(source_law, source_readback)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
