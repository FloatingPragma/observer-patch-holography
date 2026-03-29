#!/usr/bin/env python3
"""Export the quark diagonal gap-shift emitter artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MAP = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_map.json"
DEFAULT_SPREAD = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DEFAULT_SCALAR_EVALUATOR = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_scalar_evaluator.json"
DEFAULT_TAU_MAP = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_tau_map.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_emitter.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    diagonal_gap_shift_map: dict,
    spread_map: dict,
    scalar_evaluator: dict | None = None,
    tau_map: dict | None = None,
) -> dict:
    scalar_evaluator = scalar_evaluator or {}
    tau_map = tau_map or {}
    q_ord = [float(value) for value in spread_map["quadratic_residual_basis_vector"]]
    b_ord = [float(value) for value in diagonal_gap_shift_map["B_ord"]]
    return {
        "artifact": "oph_family_excitation_diagonal_gap_shift_emitter",
        "generated_utc": _timestamp(),
        "proof_status": "emitter_shell_waiting_scalar_evaluator",
        "source_artifact": diagonal_gap_shift_map.get("artifact"),
        "scalar_evaluator_artifact": scalar_evaluator.get("artifact"),
        "scalar_evaluator_status": scalar_evaluator.get("proof_status"),
        "tau_map_artifact": tau_map.get("artifact"),
        "tau_map_status": tau_map.get("proof_status"),
        "tau_emitter_status": "waiting_scalar_evaluator",
        "B_ord": b_ord,
        "orthogonality_certificate": {
            "sum_B_ord": sum(b_ord),
            "dot_B_ord_Q_ord": sum(b * q for b, q in zip(b_ord, q_ord)),
        },
        "tau_u_log_per_side": tau_map.get("tau_u_log_per_side"),
        "tau_d_log_per_side": tau_map.get("tau_d_log_per_side"),
        "tau_u_formula": tau_map.get("tau_u_formula"),
        "tau_d_formula": tau_map.get("tau_d_formula"),
        "E_u_log_effective_formula": "base_E_u_log + tau_u_log_per_side * B_ord",
        "E_d_log_effective_formula": "base_E_d_log + tau_d_log_per_side * B_ord",
        "gamma21_u_effective_formula": "gamma21_u_base + tau_u_log_per_side",
        "gamma32_u_effective_formula": "gamma32_u_base + tau_u_log_per_side",
        "gamma21_d_effective_formula": "gamma21_d_base + tau_d_log_per_side",
        "gamma32_d_effective_formula": "gamma32_d_base + tau_d_log_per_side",
        "sigma_u_effective_total_log_per_side_formula": "sigma_u_total_log_per_side + 2 * tau_u_log_per_side",
        "sigma_d_effective_total_log_per_side_formula": "sigma_d_total_log_per_side + 2 * tau_d_log_per_side",
        "promotion_rule": "predictive_promotion_allowed iff tau_map_status == 'closed'",
        "smallest_constructive_missing_object": scalar_evaluator.get(
            "smallest_constructive_missing_object",
            "beta_u_diag_B_source_and_beta_d_diag_B_source",
        ),
        "notes": [
            "The diagonal gap-shift map is the correct next-family normal form after the current surface is exhausted.",
            "The remaining predictive object is the source-side amplitude pair beta_u and beta_d beneath the closed pure-B readback law, not a broader quark family.",
            "Audit-only best-fit tau values remain quarantined and are not promoted here.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark diagonal gap-shift emitter artifact.")
    parser.add_argument("--map", default=str(DEFAULT_MAP))
    parser.add_argument("--spread-map", default=str(DEFAULT_SPREAD))
    parser.add_argument("--scalar-evaluator", default=str(DEFAULT_SCALAR_EVALUATOR))
    parser.add_argument("--tau-map", default=str(DEFAULT_TAU_MAP))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    diagonal_gap_shift_map = json.loads(Path(args.map).read_text(encoding="utf-8"))
    spread_map = json.loads(Path(args.spread_map).read_text(encoding="utf-8"))
    scalar_evaluator_path = Path(args.scalar_evaluator)
    scalar_evaluator = json.loads(scalar_evaluator_path.read_text(encoding="utf-8")) if scalar_evaluator_path.exists() else None
    tau_map_path = Path(args.tau_map)
    tau_map = json.loads(tau_map_path.read_text(encoding="utf-8")) if tau_map_path.exists() else None
    artifact = build_artifact(diagonal_gap_shift_map, spread_map, scalar_evaluator, tau_map)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
