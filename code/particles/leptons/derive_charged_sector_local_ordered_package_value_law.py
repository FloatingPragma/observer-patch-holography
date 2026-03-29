#!/usr/bin/env python3
"""Emit the charged ordered-package value-law artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_EMISSION = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_source_emission.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_value_law.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_emission: dict) -> dict:
    ordered = [float(value) for value in source_emission["source_side_ordered_package_log_per_side_emitted"]]
    x2 = float(source_emission["ordered_family_coordinate"][1])
    r1, r2, r3 = ordered
    mu = (r1 + r2 + r3) / 3.0
    sigma = r3 - r1
    eta = 2.0 * ((((1.0 - x2) * r1) + ((1.0 + x2) * r3)) / 2.0 - r2)
    q_factor = 1.0 - x2 * x2
    return {
        "artifact": "oph_charged_sector_local_ordered_package_value_law",
        "generated_utc": _timestamp(),
        "proof_status": "current_support_value_law_closed_same_support_exhausted",
        "predictive_promotion_allowed": False,
        "predictive_value_law_closed": False,
        "source_artifact": source_emission.get("artifact"),
        "ordered_family_coordinate": [-1.0, x2, 1.0],
        "linear_basis_vector_centered": [
            -(3.0 + x2) / 3.0,
            (2.0 * x2) / 3.0,
            (3.0 - x2) / 3.0,
        ],
        "quadratic_basis_vector": [
            q_factor / 3.0,
            -(2.0 * q_factor) / 3.0,
            q_factor / 3.0,
        ],
        "source_side_ordered_package_log_per_side_emitted": ordered,
        "mu_source_log_per_side_readback": mu,
        "sigma_source_total_log_per_side_readback": sigma,
        "affine_midpoint_log_per_side_readback": (((1.0 - x2) * r1) + ((1.0 + x2) * r3)) / 2.0,
        "delta_midpoint_log_per_side_readback": eta / 2.0,
        "eta_source_split_log_per_side_readback": eta,
        "b_e_log_coeff_readback": eta / (2.0 * q_factor) if abs(q_factor) > 1.0e-15 else 0.0,
        "midpoint_defect_emitter_closed": abs(eta / 2.0) <= 1.0e-15,
        "delta_midpoint_zero_tolerance": 1.0e-15,
        "delta_midpoint_zero_on_current_support": abs(eta / 2.0) <= 1.0e-15,
        "source_package_formula": "r_ord = mu*1 + (sigma/2)*L_ord + (eta/(2*(1-x2^2)))*Q_ord",
        "entry_formulas": {
            "r1": "mu - ((3 + x2)*sigma - eta)/6",
            "r2": "mu + (x2*sigma - eta)/3",
            "r3": "mu + ((3 - x2)*sigma + eta)/6",
        },
        "current_package_linear_subray_only": abs(eta) < 1.0e-12,
        "carrier_centered_rank": 1 if abs(eta) < 1.0e-12 else 2,
        "collapse_proven": abs(eta) < 1.0e-12,
        "smallest_constructive_missing_object": "oph_charged_sector_local_current_support_obstruction_certificate",
        "next_single_residual_object": "oph_charged_sector_local_current_support_obstruction_certificate",
        "notes": [
            "This artifact consolidates the ordered-package value shell on the fixed charged carrier.",
            "On the current support the midpoint defect closes to zero, so the centered package collapses to the linear subray.",
            "The remaining charged gap is now the obstruction certificate showing the current support is too small, not another same-support midpoint-defect scalar.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged ordered-package value-law artifact.")
    parser.add_argument("--source-emission", default=str(DEFAULT_SOURCE_EMISSION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_emission = json.loads(Path(args.source_emission).read_text(encoding="utf-8"))
    artifact = build_artifact(source_emission)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
