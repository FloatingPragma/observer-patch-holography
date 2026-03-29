#!/usr/bin/env python3
"""Emit the quark diagonal common gap-shift source-values artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_LAW = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_law.json"
DEFAULT_SOURCE_READBACK = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_readback.json"
DEFAULT_SOURCE_EMISSION = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_emission.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_values.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_law: dict, source_readback: dict | None = None, source_emission: dict | None = None) -> dict:
    b_ord = [float(value) for value in source_law["B_ord"]]
    x2 = float(source_law["normalized_coordinate_x2"])
    q_ord = [
        (1.0 - x2 * x2) / 3.0,
        -(2.0 * (1.0 - x2 * x2)) / 3.0,
        (1.0 - x2 * x2) / 3.0,
    ]
    source_readback = source_readback or {}
    source_emission = source_emission or {}
    beta_u = source_emission.get("beta_u_diag_B_source")
    beta_d = source_emission.get("beta_d_diag_B_source")
    delta_u = None if beta_u is None else [float(beta_u) * value for value in b_ord]
    delta_d = None if beta_d is None else [float(beta_d) * value for value in b_ord]
    mean_u = None if delta_u is None else sum(delta_u)
    mean_d = None if delta_d is None else sum(delta_d)
    quad_u = None if delta_u is None else sum(u * q for u, q in zip(delta_u, q_ord))
    quad_d = None if delta_d is None else sum(d * q for d, q in zip(delta_d, q_ord))
    return {
        "artifact": "oph_family_excitation_diagonal_common_gap_shift_source_values",
        "generated_utc": _timestamp(),
        "proof_status": (
            "source_values_derived_from_source_emission"
            if beta_u is not None and beta_d is not None
            else "source_values_shell_waiting_pure_B_payload_pair"
        ),
        "predictive_promotion_allowed": False,
        "source_artifact": source_emission.get("artifact", source_law.get("artifact")),
        "source_law_artifact": source_law.get("artifact"),
        "source_readback_artifact": source_readback.get("artifact"),
        "source_readback_status": source_readback.get("proof_status"),
        "source_emission_artifact": source_emission.get("artifact"),
        "source_emission_status": source_emission.get("proof_status"),
        "B_ord": b_ord,
        "Q_ord": q_ord,
        "beta_u_diag_B_source": beta_u,
        "beta_d_diag_B_source": beta_d,
        "delta_E_u_diag_log_per_side_source": delta_u,
        "delta_E_d_diag_log_per_side_source": delta_d,
        "mean_annihilator_u": mean_u,
        "quadratic_annihilator_u": quad_u,
        "mean_annihilator_d": mean_d,
        "quadratic_annihilator_d": quad_d,
        "pure_B_mode_certificate_u": (
            abs(mean_u) < 1.0e-12 and abs(quad_u) < 1.0e-12
            if mean_u is not None and quad_u is not None
            else None
        ),
        "pure_B_mode_certificate_d": (
            abs(mean_d) < 1.0e-12 and abs(quad_d) < 1.0e-12
            if mean_d is not None and quad_d is not None
            else None
        ),
        "delta_E_u_diag_log_per_side_source_formula": "beta_u_diag_B_source * B_ord",
        "delta_E_d_diag_log_per_side_source_formula": "beta_d_diag_B_source * B_ord",
        "mean_annihilator_u_formula": "0.0",
        "quadratic_annihilator_u_formula": "0.0",
        "mean_annihilator_d_formula": "0.0",
        "quadratic_annihilator_d_formula": "0.0",
        "pure_B_mode_certificate_u_formula": "true",
        "pure_B_mode_certificate_d_formula": "true",
        "smallest_constructive_missing_object": source_readback.get(
            "smallest_constructive_missing_object",
            "source_readback_u_log_per_side_and_source_readback_d_log_per_side",
        ),
        "notes": [
            "This artifact derives the B-mode source values from the closed pure-B source-readback law when the beta-pair exists.",
            "The predictive gap is the emitted pure-B payload pair beneath this derived shell; the odd projector outputs are algebraic consequences of that payload.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark diagonal common gap-shift source-values artifact.")
    parser.add_argument("--source-law", default=str(DEFAULT_SOURCE_LAW))
    parser.add_argument("--source-readback", default=str(DEFAULT_SOURCE_READBACK))
    parser.add_argument("--source-emission", default=str(DEFAULT_SOURCE_EMISSION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_law = json.loads(Path(args.source_law).read_text(encoding="utf-8"))
    source_readback_path = Path(args.source_readback)
    source_readback = json.loads(source_readback_path.read_text(encoding="utf-8")) if source_readback_path.exists() else None
    source_emission_path = Path(args.source_emission)
    source_emission = json.loads(source_emission_path.read_text(encoding="utf-8")) if source_emission_path.exists() else None
    artifact = build_artifact(source_law, source_readback, source_emission)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
