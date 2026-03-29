#!/usr/bin/env python3
"""Emit the quark diagonal common gap-shift source-law artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MAP = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_map.json"
DEFAULT_SPREAD = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_common_gap_shift_source_law.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(diagonal_gap_shift_map: dict, spread_map: dict) -> dict:
    b_ord = [float(value) for value in diagonal_gap_shift_map["B_ord"]]
    x2 = float(spread_map["normalized_coordinate_x2"])
    return {
        "artifact": "oph_family_excitation_diagonal_common_gap_shift_source_law",
        "generated_utc": _timestamp(),
        "proof_status": "source_law_closed_waiting_J_B_source_pair",
        "predictive_promotion_allowed": False,
        "source_artifact": diagonal_gap_shift_map.get("artifact"),
        "input_kind": "closed_current_surface_plus_diagonal_B_mode",
        "output_kind": "sector_diagonal_common_gap_shift_source_readback",
        "B_ord": b_ord,
        "B_ord_norm_sq": sum(value * value for value in b_ord),
        "normalized_coordinate_x2": x2,
        "delta_E_u_diag_log_per_side_source": None,
        "delta_E_d_diag_log_per_side_source": None,
        "tau_u_log_per_side": None,
        "tau_d_log_per_side": None,
        "equal_gap_residual_u": None,
        "equal_gap_residual_d": None,
        "delta_gamma21_u_log_per_side": None,
        "delta_gamma32_u_log_per_side": None,
        "delta_gamma21_d_log_per_side": None,
        "delta_gamma32_d_log_per_side": None,
        "delta_sigma_u_diag_total_log_per_side": None,
        "delta_sigma_d_diag_total_log_per_side": None,
        "tau_u_formula": "dot(delta_E_u_diag_log_per_side_source, B_ord) / B_ord_norm_sq",
        "tau_d_formula": "dot(delta_E_d_diag_log_per_side_source, B_ord) / B_ord_norm_sq",
        "equal_gap_residual_u_formula": "delta_E_u_diag_log_per_side_source - tau_u_log_per_side * B_ord",
        "equal_gap_residual_d_formula": "delta_E_d_diag_log_per_side_source - tau_d_log_per_side * B_ord",
        "delta_gamma21_u_log_per_side_formula": "tau_u_log_per_side",
        "delta_gamma32_u_log_per_side_formula": "tau_u_log_per_side",
        "delta_gamma21_d_log_per_side_formula": "tau_d_log_per_side",
        "delta_gamma32_d_log_per_side_formula": "tau_d_log_per_side",
        "delta_sigma_u_diag_total_log_per_side_formula": "2 * tau_u_log_per_side",
        "delta_sigma_d_diag_total_log_per_side_formula": "2 * tau_d_log_per_side",
        "smallest_constructive_missing_object": "J_B_source_u_and_J_B_source_d",
        "notes": [
            "This artifact fixes the source-side law beneath the quark diagonal gap-shift shells.",
            "The present invariants already fix the pure-B source-readback law; they do not yet emit the odd projector values that carry that law.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark diagonal common gap-shift source-law artifact.")
    parser.add_argument("--map", default=str(DEFAULT_MAP))
    parser.add_argument("--spread-map", default=str(DEFAULT_SPREAD))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    diagonal_gap_shift_map = json.loads(Path(args.map).read_text(encoding="utf-8"))
    spread_map = json.loads(Path(args.spread_map).read_text(encoding="utf-8"))
    artifact = build_artifact(diagonal_gap_shift_map, spread_map)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
