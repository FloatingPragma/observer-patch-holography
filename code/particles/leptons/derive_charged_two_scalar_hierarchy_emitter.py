#!/usr/bin/env python3
"""Export the next charged-lepton hierarchy emitter artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_READOUT = ROOT / "particles" / "runs" / "leptons" / "lepton_log_spectrum_readout.json"
DEFAULT_AUDIT = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exactness_audit.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_two_scalar_hierarchy_emitter.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(readout: dict, audit: dict) -> dict:
    x2 = float(readout["ordered_family_coordinate"][1])
    return {
        "artifact": "oph_charged_current_family_two_scalar_hierarchy_emitter",
        "generated_utc": _timestamp(),
        "theorem_scope": "current_family_only",
        "hierarchy_emitter_status": "missing_joint_emission",
        "ordered_family_coordinate": list(readout["ordered_family_coordinate"]),
        "rho_ord": float(readout["rho_ord"]),
        "sigma_e_total_log_per_side_emitted": None,
        "eta_e_split_log_per_side_emitted": None,
        "gamma21_log_per_side_emitted_formula": "((1 + x2) * sigma_e_total_log_per_side_emitted - eta_e_split_log_per_side_emitted) / 2",
        "gamma32_log_per_side_emitted_formula": "((1 - x2) * sigma_e_total_log_per_side_emitted + eta_e_split_log_per_side_emitted) / 2",
        "E_e_log_centered_emitted_formula": [
            "-(2*gamma21_log_per_side_emitted + gamma32_log_per_side_emitted)/3",
            "(gamma21_log_per_side_emitted - gamma32_log_per_side_emitted)/3",
            "(gamma21_log_per_side_emitted + 2*gamma32_log_per_side_emitted)/3",
        ],
        "diagonal_emitter_rule": "Y_e_shape = diag(exp(E_e_log_centered_emitted))",
        "frozen_sigma_branch_impossible": True,
        "upstream_predictive_object": "oph_charged_lepton_excitation_gap_map",
        "impossibility_certificate": {
            "sigma_e_total_log_per_side_current": float(audit["centered_hierarchy_audit"]["current_sigma_total_log_per_side"]),
            "sigma_e_total_log_per_side_target": float(audit["centered_hierarchy_audit"]["target_sigma_total_log_per_side"]),
            "centered_residual_norm_current": float(audit["centered_hierarchy_audit"]["residual_norm"]),
            "centered_residual_norm_after_best_eta_on_current_sigma": float(
                audit["two_scalar_on_current_sigma_audit"]["residual_norm_after_best_eta_on_current_sigma"]
            ),
        },
        "smallest_constructive_missing_object": "oph_charged_lepton_excitation_gap_map",
        "notes": [
            "The charged lane now needs the upstream excitation-gap map that emits the two sector-even coefficients together.",
            "The current sigma-only branch is explicitly marked impossible as an exact charged-hierarchy branch.",
            "Any diagnostic eta values remain quarantined in the audit and are not promoted here.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged two-scalar hierarchy-emitter artifact.")
    parser.add_argument("--readout", default=str(DEFAULT_READOUT))
    parser.add_argument("--audit", default=str(DEFAULT_AUDIT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    readout = json.loads(Path(args.readout).read_text(encoding="utf-8"))
    audit = json.loads(Path(args.audit).read_text(encoding="utf-8"))
    artifact = build_artifact(readout, audit)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
