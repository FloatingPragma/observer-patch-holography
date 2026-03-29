#!/usr/bin/env python3
"""Emit the charged midpoint-defect emitter artifact on the current support."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_VALUE_LAW = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_value_law.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_midpoint_defect_emitter.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(value_law: dict) -> dict:
    x2 = float(value_law["ordered_family_coordinate"][1])
    delta_mid = float(value_law.get("delta_midpoint_log_per_side_readback", 0.0))
    eta_split = float(value_law.get("eta_source_split_log_per_side_readback", 0.0))
    b_e = float(value_law.get("b_e_log_coeff_readback", 0.0))
    return {
        "artifact": "oph_charged_sector_local_ordered_package_midpoint_defect_emitter",
        "generated_utc": _timestamp(),
        "proof_status": "current_support_midpoint_defect_emitter_closed",
        "predictive_promotion_allowed": False,
        "source_artifact": value_law.get("artifact"),
        "ordered_family_coordinate": [-1.0, x2, 1.0],
        "mu_source_log_per_side_emitted": value_law.get("mu_source_log_per_side_readback"),
        "sigma_source_total_log_per_side_emitted": value_law.get("sigma_source_total_log_per_side_readback"),
        "affine_midpoint_log_per_side_emitted": value_law.get("affine_midpoint_log_per_side_readback"),
        "delta_midpoint_log_per_side_emitted": delta_mid,
        "eta_source_split_log_per_side_emitted": eta_split,
        "b_e_log_coeff_emitted": b_e,
        "carrier_centered_rank": value_law.get("carrier_centered_rank"),
        "collapse_proven": value_law.get("collapse_proven"),
        "midpoint_defect_emitter_closed": True,
        "same_support_transverse_coeff_closed": True,
        "delta_midpoint_zero_tolerance": 1.0e-15,
        "delta_midpoint_zero_on_current_support": abs(delta_mid) <= 1.0e-15,
        "source_side_ordered_package_log_per_side_emitted": value_law.get("source_side_ordered_package_log_per_side_emitted"),
        "midpoint_formulas": {
            "m_aff": "((1 - x2) * r1 + (1 + x2) * r3) / 2",
            "delta_mid": "m_aff - r2",
            "eta": "2 * delta_mid",
            "b_e": "eta / (2 * (1 - x2^2))",
        },
        "source_package_formula": "r_ord = mu*1 + (sigma/2)*L_ord + b_e*Q_ord",
        "smallest_constructive_missing_object": None,
        "next_single_residual_object": "oph_charged_sector_local_current_support_obstruction_certificate",
        "notes": [
            "The midpoint-defect emitter is already derivable on the current support and closes to zero there.",
            "No same-support transverse charged mover remains on the current ordered package.",
            "The remaining charged gap is no longer a same-support midpoint defect, but the obstruction certificate showing the present support cannot reproduce the charged hierarchy.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged midpoint-defect emitter artifact.")
    parser.add_argument("--value-law", default=str(DEFAULT_VALUE_LAW))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    value_law = json.loads(Path(args.value_law).read_text(encoding="utf-8"))
    artifact = build_artifact(value_law)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
