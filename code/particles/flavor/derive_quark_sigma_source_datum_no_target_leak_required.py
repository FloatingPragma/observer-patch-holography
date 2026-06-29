#!/usr/bin/env python3
"""Record the missing source-only quark sigma selector theorem.

This artifact is intentionally not a prediction artifact. It pins down the
strict promotion gate for the selected-class quark lane: the downstream affine
and Yukawa algebra is closed once a physical sigma datum is supplied, but the
current exact sigma datum is inherited from the current-family target surface.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
STRENGTHENED_SIGMA_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"
)
EDGE_CANDIDATE_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_edge_statistics_spread_candidate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_sigma_source_datum_no_target_leak_required.json"

MISSING_FOR_PROMOTION = [
    "QUARK_SIGMA_SOURCE_QUOTIENT",
    "QUARK_SIGMA_SOURCE_SELECTOR",
    "QUARK_EDGE_STATISTICS_CORRECTION_THEOREM",
    "QUARK_SIGMA_REFINEMENT_COMPATIBILITY",
    "NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE",
]

FORBIDDEN_ANCESTORS = [
    "particle_reference_values",
    "quark_current_family_exact_readout",
    "quark_current_family_exact_sigma_target",
    "quark_current_family_exact_pdg_theorem",
    "quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem",
    "target_centered_log_u",
    "target_centered_log_d",
    "reference_targets_u",
    "reference_targets_d",
    "exact_fit_residuals",
    "PDG",
    "CODATA",
    "compare_only",
]

CANONICAL_EDGE_TARGET_RESIDUALS = {
    "required_R_u": -0.004490377677282886,
    "required_R_d": -0.1247947151634663,
    "required_R_seed": -0.06464254642037481,
    "required_R_eta": 0.060152168743091705,
}


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(strengthened_sigma: dict[str, Any], edge_candidate: dict[str, Any]) -> dict[str, Any]:
    target = dict(strengthened_sigma["theorem_grade_physical_sigma_datum"])
    candidate_sigmas = dict(edge_candidate["candidate_sigmas"])
    ordered_inputs = dict(edge_candidate["ordered_family_inputs"])
    edge_inputs = dict(edge_candidate["edge_statistics_inputs"]["shared_charged_suppression_seed"])

    sigma_u_edge = float(candidate_sigmas["sigma_u_total_log_per_side"])
    sigma_d_edge = float(candidate_sigmas["sigma_d_total_log_per_side"])
    sigma_u = float(target["sigma_u"])
    sigma_d = float(target["sigma_d"])
    sigma_seed = float(target["sigma_seed_ud"])
    eta_ud = float(target["eta_ud"])
    required_r_u = CANONICAL_EDGE_TARGET_RESIDUALS["required_R_u"]
    required_r_d = CANONICAL_EDGE_TARGET_RESIDUALS["required_R_d"]
    required_r_seed = CANONICAL_EDGE_TARGET_RESIDUALS["required_R_seed"]
    required_r_eta = CANONICAL_EDGE_TARGET_RESIDUALS["required_R_eta"]

    return {
        "artifact": "oph_quark_sigma_source_datum_no_target_leak_required",
        "generated_utc": _timestamp(),
        "status": "missing_theorem",
        "proof_status": "missing_theorem",
        "claim_tier": "selected_class_conditional_on_source_sigma",
        "promotion_allowed": False,
        "public_promotion_allowed": False,
        "source_only_sigma_emitted": False,
        "downstream_algebra_closed": True,
        "selected_class": "f_P",
        "required_identity": "P -> f_P -> Sigma_ud_phys -> (sigma_u,sigma_d,sigma_seed_ud,eta_ud)",
        "current_exact_path_is_target_derived": [
            "quark_current_family_exact_readout",
            "quark_current_family_exact_sigma_target",
            "quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem",
            "quark_public_physical_sigma_datum_descent",
            "quark_public_exact_yukawa_end_to_end_theorem",
        ],
        "target_values_for_future_source_theorem": {
            "sigma_u": sigma_u,
            "sigma_d": sigma_d,
            "sigma_seed_ud": sigma_seed,
            "eta_ud": eta_ud,
        },
        "strongest_current_source_candidate": {
            "formula_sigma_u_edge": "S_13 + rho_ord*delta21/(1+rho_ord)",
            "formula_sigma_d_edge": "S_23 + delta21/(2*(1+rho_ord-x2^2))",
            "sigma_u_edge": sigma_u_edge,
            "sigma_d_edge": sigma_d_edge,
            "required_R_u": required_r_u,
            "required_R_d": required_r_d,
            "required_R_seed": required_r_seed,
            "required_R_eta": required_r_eta,
            "S13": float(edge_inputs["S_13"]),
            "S23": float(edge_inputs["S_23"]),
            "rho_ord": float(ordered_inputs["rho_ord"]),
            "delta21": float(ordered_inputs["delta21"]),
            "x2": float(ordered_inputs["x2"]),
            "mean_denominator": float(ordered_inputs["mean_denominator"]),
        },
        "conditional_downstream_if_closed": {
            "sigma_seed_ud": "(sigma_u + sigma_d)/2",
            "eta_ud": "(sigma_u - sigma_d)/2",
            "g_u": "g_ch*exp(-(A_ud*sigma_seed_ud - B_ud*eta_ud))",
            "g_d": "g_ch*exp(-(A_ud*sigma_seed_ud + B_ud*eta_ud))",
            "then": "ordered_three_point_readout_and_exact_forward_Yu_Yd",
        },
        "missing_for_promotion": MISSING_FOR_PROMOTION,
        "forbidden_ancestors": FORBIDDEN_ANCESTORS,
        "required_theorem_package": {
            "QΣ-0": "sigma descent is representative independence, not source selection",
            "QΣ-1": "selected-frame quark masses remain underdetermined without a source sigma selector",
            "QΣ-2": "conditional quark masses follow once a no-target source sigma datum is supplied",
            "QΣ-A": "QUARK_SIGMA_SOURCE_QUOTIENT",
            "QΣ-B": "QUARK_SIGMA_SOURCE_SELECTOR",
            "QΣ-C": "QUARK_EDGE_STATISTICS_CORRECTION_THEOREM",
            "QΣ-D": "QUARK_SIGMA_REFINEMENT_COMPATIBILITY",
            "QΣ-E": "NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE",
        },
        "notes": [
            "The edge-statistics candidate is retained as the strongest current source-side theorem target.",
            "The required corrections may not be defined as residuals against running-quark target values.",
            "This artifact blocks promotion; it does not emit a replacement sigma datum.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the missing quark sigma source theorem gate artifact.")
    parser.add_argument("--strengthened-sigma", default=str(STRENGTHENED_SIGMA_JSON))
    parser.add_argument("--edge-candidate", default=str(EDGE_CANDIDATE_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(_load_json(Path(args.strengthened_sigma)), _load_json(Path(args.edge_candidate)))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
