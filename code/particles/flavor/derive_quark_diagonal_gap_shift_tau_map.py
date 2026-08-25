#!/usr/bin/env python3
"""Export the quark diagonal gap-shift tau-map artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MAP = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_map.json"
DEFAULT_SPREAD = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DEFAULT_SCALAR_EVALUATOR = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_scalar_evaluator.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_tau_map.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    diagonal_gap_shift_map: dict,
    spread_map: dict,
    scalar_evaluator: dict | None = None,
) -> dict:
    q_ord = [float(value) for value in spread_map["quadratic_residual_basis_vector"]]
    b_ord = [float(value) for value in diagonal_gap_shift_map["B_ord"]]
    scalar_evaluator = scalar_evaluator or {}
    source_value_promotion_ready = (
        scalar_evaluator.get("source_value_promotion_ready") is True
        and scalar_evaluator.get("value_classification") == "source_emitted"
    )
    return {
        "artifact": "oph_family_excitation_diagonal_gap_shift_tau_map",
        "generated_utc": _timestamp(),
        "proof_status": (
            scalar_evaluator.get("proof_status")
            if scalar_evaluator.get("artifact") == "oph_family_excitation_diagonal_gap_shift_scalar_evaluator"
            else "missing_predictive_law"
        ),
        "parent_artifact": scalar_evaluator.get("artifact", "oph_family_excitation_diagonal_gap_shift_emitter"),
        "input_kind": "ordered_branch_generator_spectral_package_plus_closed_current_surface",
        "output_kind": "sector_diagonal_gap_shift_coefficients",
        "scalar_evaluator_artifact": scalar_evaluator.get("artifact"),
        "scalar_evaluator_status": scalar_evaluator.get("proof_status"),
        "source_value_promotion_ready": source_value_promotion_ready,
        "value_classification": scalar_evaluator.get(
            "value_classification", "absent"
        ),
        "promotion_blockers": scalar_evaluator.get("promotion_blockers", []),
        "B_ord": b_ord,
        "tau_u_log_per_side": scalar_evaluator.get("tau_u_log_per_side"),
        "tau_d_log_per_side": scalar_evaluator.get("tau_d_log_per_side"),
        "tau_u_formula": scalar_evaluator.get("tau_u_formula"),
        "tau_d_formula": scalar_evaluator.get("tau_d_formula"),
        "identifiability_certificate": {
            "sum_B_ord": sum(b_ord),
            "dot_B_ord_Q_ord": sum(b * q for b, q in zip(b_ord, q_ord)),
            "mean_blind_to_tau": True,
            "quadratic_mode_blind_to_tau": True
        },
        "predictive_promotion_allowed": False,
        "smallest_constructive_missing_object": scalar_evaluator.get(
            "smallest_constructive_missing_object",
            "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
        ),
        "notes": [
            "The diagonal residual family is fixed, and the tau-map is the first data-bearing shell consumed by the active quark builder.",
            "The source-side beta difference fixes Delta_ud_overlap; the D12 overlap map then combines it with a fixed positive sigma branch to produce tau_u and tau_d.",
            "A comparison-only beta pair or sigma branch is carried for arithmetic auditing but cannot promote the tau-map.",
            "No PDG quark values are consumed here."
        ]
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark diagonal gap-shift tau-map artifact.")
    parser.add_argument("--map", default=str(DEFAULT_MAP))
    parser.add_argument("--spread-map", default=str(DEFAULT_SPREAD))
    parser.add_argument("--scalar-evaluator", default=str(DEFAULT_SCALAR_EVALUATOR))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    diagonal_gap_shift_map = json.loads(Path(args.map).read_text(encoding="utf-8"))
    spread_map = json.loads(Path(args.spread_map).read_text(encoding="utf-8"))
    scalar_evaluator_path = Path(args.scalar_evaluator)
    scalar_evaluator = (
        json.loads(scalar_evaluator_path.read_text(encoding="utf-8"))
        if scalar_evaluator_path.exists()
        else None
    )
    artifact = build_artifact(diagonal_gap_shift_map, spread_map, scalar_evaluator)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
