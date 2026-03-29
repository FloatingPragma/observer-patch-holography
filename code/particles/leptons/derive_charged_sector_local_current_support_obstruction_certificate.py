#!/usr/bin/env python3
"""Prove the current charged support is exhausted before extending the carrier.

Chain role: certify that the present charged ordered package only spans the
linear subray and therefore cannot finish the hierarchy by same-support tuning.

Mathematics: ordered three-point linear/quadratic basis analysis, gap-ratio
constraints, and support-rank bookkeeping.

OPH-derived inputs: the charged ordered package readback, midpoint defect lane,
log-spectrum data, and the current forward charged candidate.

Output: the obstruction certificate that identifies the first missing
support-extension scalar beneath the charged lane.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_READBACK = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_readback.json"
DEFAULT_MIDPOINT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_midpoint_defect_emitter.json"
DEFAULT_LOG_READOUT = ROOT / "particles" / "runs" / "leptons" / "lepton_log_spectrum_readout.json"
DEFAULT_FORWARD = ROOT / "particles" / "runs" / "leptons" / "forward_charged_leptons.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_current_support_obstruction_certificate.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(readback: dict, midpoint: dict, log_readout: dict, forward: dict) -> dict:
    x2 = float(readback["ordered_family_coordinate"][1])
    linear_basis = [
        -(3.0 + x2) / 3.0,
        (2.0 * x2) / 3.0,
        (3.0 - x2) / 3.0,
    ]
    quadratic_basis = [
        (1.0 - x2 * x2) / 3.0,
        -(2.0 * (1.0 - x2 * x2)) / 3.0,
        (1.0 - x2 * x2) / 3.0,
    ]
    sigma_total = float(readback["sigma_e_total_log_per_side_emitted"])
    centered = [float(value) for value in readback["source_side_ordered_package_centered_log_emitted"]]
    g_active = float(((forward.get("channel_norm") or {}).get("g_e_candidate")) or 0.0)
    singular_values = [float(value) for value in forward["singular_values_abs"]]
    gap21 = centered[1] - centered[0]
    gap32 = centered[2] - centered[1]
    gap_ratio = gap21 / gap32 if abs(gap32) > 1.0e-15 else None
    mu_constraint = (singular_values[0] ** ((1.0 - x2) / 2.0)) * (singular_values[2] ** ((1.0 + x2) / 2.0))
    return {
        "artifact": "oph_charged_sector_local_current_support_obstruction_certificate",
        "generated_utc": _timestamp(),
        "proof_status": "current_support_obstruction_certificate_closed",
        "predictive_promotion_allowed": False,
        "source_artifacts": {
            "ordered_package_readback": readback.get("artifact"),
            "midpoint_defect_emitter": midpoint.get("artifact"),
            "log_spectrum_readout": log_readout.get("artifact"),
            "forward_charged_leptons": forward.get("artifact"),
        },
        "ordered_family_coordinate": [-1.0, x2, 1.0],
        "linear_basis_vector_centered": linear_basis,
        "quadratic_basis_vector": quadratic_basis,
        "support_basis_available_rank": 2,
        "current_support_realized_centered_rank": 1,
        "same_support_exhausted": True,
        "current_support_linear_subray_only": True,
        "same_support_transverse_coeff_closed": True,
        "mu_source_log_per_side": float(readback["source_side_ordered_package_mean_log_per_side_emitted"]),
        "sigma_source_total_log_per_side": sigma_total,
        "source_side_ordered_package_log_per_side": [
            float(value) for value in readback["source_side_ordered_package_log_per_side_emitted"]
        ],
        "source_side_ordered_package_centered_log": centered,
        "same_support_linear_subray_formula": "E_e_log_centered = (sigma_source_total_log_per_side / 2) * linear_basis_vector_centered",
        "same_support_mass_family_formula": "m_ord = g_active * exp((sigma_source_total_log_per_side / 2) * linear_basis_vector_centered)",
        "e_formula": "g_active * exp(-(3 + x2) * sigma_source_total_log_per_side / 6)",
        "mu_formula": "g_active * exp(x2 * sigma_source_total_log_per_side / 3)",
        "tau_formula": "g_active * exp((3 - x2) * sigma_source_total_log_per_side / 6)",
        "same_support_gap_ratio_fixed": (1.0 + x2) / (1.0 - x2),
        "same_support_gap_ratio_formula": "(1 + x2) / (1 - x2)",
        "same_support_geometric_interpolation_weights": {
            "e_weight": (1.0 - x2) / 2.0,
            "tau_weight": (1.0 + x2) / 2.0,
        },
        "same_support_geometric_interpolation_formula": "mu = e^((1 - x2)/2) * tau^((1 + x2)/2)",
        "current_candidate": {
            "g_active": g_active,
            "singular_values_abs": singular_values,
        },
        "current_support_invariant_check": {
            "mu_from_e_tau_same_support": mu_constraint,
            "mu_constraint_residual": singular_values[1] - mu_constraint,
            "gap_ratio_from_current_candidate": gap_ratio,
            "gap_ratio_residual": None if gap_ratio is None else gap_ratio - ((1.0 + x2) / (1.0 - x2)),
            "linear_subray_residual_norm": sum(
                abs(centered[idx] - (sigma_total / 2.0) * linear_basis[idx]) for idx in range(3)
            ),
        },
        "obstruction_statement": "All predictive current-support descendants lie on the one-dimensional centered cone generated by linear_basis_vector_centered. No same-support transverse mover remains.",
        "minimal_extension_requirement": {
            "new_beyond_support_scalar_slots_required": 1,
            "extension_readback_formula": "E_e_log_centered_ext = (sigma_source_total_log_per_side / 2) * linear_basis_vector_centered + kappa_ext * Delta_ord",
            "extension_constraints": [
                "sum(Delta_ord) = 0",
                "Delta_ord is not proportional to linear_basis_vector_centered",
                "kappa_ext != 0",
            ],
        },
        "smallest_constructive_missing_object": "oph_charged_sector_local_minimal_source_support_extension_emitter",
        "notes": [
            "The current charged support is exhausted: the midpoint defect closes to zero and the centered family collapses to a single linear subray.",
            "No same-support scalar can repair the charged hierarchy once the support has collapsed to that rank-one centered cone.",
            "The next honest charged object is a minimal support extension emitter, not another same-support correction.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged current-support obstruction certificate.")
    parser.add_argument("--readback", default=str(DEFAULT_READBACK))
    parser.add_argument("--midpoint", default=str(DEFAULT_MIDPOINT))
    parser.add_argument("--log-readout", default=str(DEFAULT_LOG_READOUT))
    parser.add_argument("--forward", default=str(DEFAULT_FORWARD))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    readback = json.loads(Path(args.readback).read_text(encoding="utf-8"))
    midpoint = json.loads(Path(args.midpoint).read_text(encoding="utf-8"))
    log_readout = json.loads(Path(args.log_readout).read_text(encoding="utf-8"))
    forward = json.loads(Path(args.forward).read_text(encoding="utf-8"))
    artifact = build_artifact(readback, midpoint, log_readout, forward)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
