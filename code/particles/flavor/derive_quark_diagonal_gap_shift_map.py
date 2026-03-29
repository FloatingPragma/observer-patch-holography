#!/usr/bin/env python3
"""Export the next quark family beyond the exhausted current quadratic surface."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SPREAD_MAP = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DEFAULT_AUDIT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exactness_audit.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_map.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(spread_map: dict, audit: dict) -> dict:
    gap_u = dict(spread_map["gap_pair_u"])
    gap_d = dict(spread_map["gap_pair_d"])
    return {
        "artifact": "oph_family_excitation_diagonal_gap_shift_map",
        "generated_utc": _timestamp(),
        "proof_status": "next_family_declared_after_current_surface_obstruction",
        "predictive_promotion_allowed": False,
        "surface_origin_artifact": spread_map.get("artifact"),
        "surface_obstruction_artifact": audit.get("artifact"),
        "surface_exhausted": True,
        "B_ord": [-1.0, 0.0, 1.0],
        "tau_u_log_per_side": None,
        "tau_d_log_per_side": None,
        "gamma21_u_base": float(gap_u["gamma21_log_per_side"]),
        "gamma32_u_base": float(gap_u["gamma32_log_per_side"]),
        "gamma21_d_base": float(gap_d["gamma21_log_per_side"]),
        "gamma32_d_base": float(gap_d["gamma32_log_per_side"]),
        "gamma21_u_lifted_formula": "gamma21_u_base + tau_u_log_per_side",
        "gamma32_u_lifted_formula": "gamma32_u_base + tau_u_log_per_side",
        "gamma21_d_lifted_formula": "gamma21_d_base + tau_d_log_per_side",
        "gamma32_d_lifted_formula": "gamma32_d_base + tau_d_log_per_side",
        "sigma_u_effective_total_log_per_side_formula": "sigma_u_total_log_per_side + 2 * tau_u_log_per_side",
        "sigma_d_effective_total_log_per_side_formula": "sigma_d_total_log_per_side + 2 * tau_d_log_per_side",
        "E_u_log_lifted_formula": "E_u_log + tau_u_log_per_side * B_ord",
        "E_d_log_lifted_formula": "E_d_log + tau_d_log_per_side * B_ord",
        "base_surface": {
            "sigma_u_total_log_per_side": float(spread_map["sigma_u_total_log_per_side"]),
            "sigma_d_total_log_per_side": float(spread_map["sigma_d_total_log_per_side"]),
            "E_u_log": [float(value) for value in spread_map["E_u_log"]],
            "E_d_log": [float(value) for value in spread_map["E_d_log"]],
        },
        "audit_only_best_fit": {
            "tau_u_best_fit": float(audit["diagonal_gap_shift_audit"]["tau_u_best_fit"]),
            "tau_d_best_fit": float(audit["diagonal_gap_shift_audit"]["tau_d_best_fit"]),
        },
        "smallest_constructive_missing_object": "oph_family_excitation_diagonal_gap_shift_emitter",
        "notes": [
            "This is the smallest next family after the current closed quadratic surface is exhausted.",
            "It does not reopen the old four-gap branch; it adds one diagonal gap-shift slot per sector.",
            "Best-fit tau values from the audit remain diagnostic-only and must not be promoted into the predictive builder.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the diagonal gap-shift quark family artifact.")
    parser.add_argument("--spread-map", default=str(DEFAULT_SPREAD_MAP))
    parser.add_argument("--audit", default=str(DEFAULT_AUDIT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    spread_map = json.loads(Path(args.spread_map).read_text(encoding="utf-8"))
    audit = json.loads(Path(args.audit).read_text(encoding="utf-8"))
    artifact = build_artifact(spread_map, audit)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
