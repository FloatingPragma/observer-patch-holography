#!/usr/bin/env python3
"""Emit the conditional D11 Higgs/top split on its declared surface.

Chain role: replace the old target-anchored split-pair assembler with one
forward algebraic map on the unpromoted D10 repair candidate and the declared
D11 core/Jacobian surface.

Mathematics: the reference central pair sits 0.7 sigma off the one-scalar
fixed ray on the declared linear surface, and the two-coordinate split of the
D11 lane is a chosen extension. The conditional map evaluates the pair
`(pi_y, pi_lambda)` from the D10 candidate tuple
`(eta_source, beta_EW, lambda_EW, tau2_tree_exact, delta_n_tree_exact)` via
one integrated shared scalar `rho_HT = log(1 + tau2_tree_exact)` and two
declared residual selectors.  Their source uniqueness and deformation
rigidity are not proved here.

Inputs: the D10 source transport pair, the target-free-input D10 repair
candidate, the declared D11 calibration surface, and the fixed-ray point
statement.  The D10 selector, D11 core/Jacobian provenance, absolute scale, and
physical-pole attachment remain upstream gates.

Output: a machine-readable conditional split-pair artifact that is exact only
as an implication on the declared numerical surface.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import canonicalize_artifact_paths
D10_SOURCE_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
D10_REPAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json"
D11_SURFACE_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_declared_calibration_surface.json"
D11_NO_GO_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_fixed_ray_no_go_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_split_pair_theorem.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(d10_source: dict, d10_repair: dict, d11_surface: dict, no_go: dict) -> dict:
    eta_source = float(d10_source["eta_source"])
    beta_ew = float(d10_source["population_basis"]["beta_EW"])
    lambda_ew = float(d10_repair["basis"]["lambda_EW"])
    tau2_tree_exact = float(d10_repair["repair_chart"]["tau2_tree_exact"])
    delta_n_tree_exact = float(d10_repair["repair_chart"]["delta_n_tree_exact"])

    core = dict(d11_surface["core"])
    jacobian = dict(d11_surface["jacobian"])

    a_t = 1.5 + beta_ew / 4.0
    b_h = (4.0 / 3.0) - beta_ew / 54.0
    rho_ht = math.log1p(tau2_tree_exact)

    # These are the current declared residual-selector formulas on the live
    # D10 chart.  They reproduce the stored Higgs/top codomain on the current
    # float surface, but this module does not prove their source uniqueness,
    # prospective selection, or rigidity under admissible deformations.
    top_residual = (
        -tau2_tree_exact * (eta_source**2)
        + (1.0 + beta_ew / 28.0) * (eta_source**6)
        + (eta_source**8) / 14.0
        + (eta_source**9) / 27.0
    )
    higgs_residual = (
        (eta_source**5)
        - (3.0 / 25.0) * (eta_source**6)
        + lambda_ew * (eta_source**6) / 18.0
        + (eta_source**8) / (2.0 * beta_ew)
    )

    pi_y = (
        eta_source
        + a_t * rho_ht
        + top_residual
    ) / math.sqrt(math.pi)
    pi_lambda = (
        eta_source
        - b_h * rho_ht
        + higgs_residual
    ) / math.sqrt(math.pi)

    delta_y = float(core["y_t_core_mt"]) * pi_y
    delta_lambda = -(16.0 / 9.0) * float(core["lambda_core_mt"]) * pi_lambda
    mt = float(core["mt_pole_core_gev"]) + float(jacobian["d_mt_pole_d_y_t"]) * delta_y
    mh = float(core["mH_core_gev"]) + float(jacobian["d_mH_d_lambda"]) * delta_lambda

    sigma_exact = 0.5 * (pi_y + pi_lambda)
    eta_exact = 0.5 * (pi_y - pi_lambda)
    c_t_live = (a_t * (rho_ht - tau2_tree_exact) + top_residual) / delta_n_tree_exact
    c_h_live = ((-(4.0 / 3.0) * tau2_tree_exact) + b_h * rho_ht - higgs_residual) / delta_n_tree_exact
    d10_repair_gate_closed = (
        d10_repair.get("status") == "closed"
        and d10_repair.get("promotion_allowed") is True
        and str(d10_repair.get("proof_status", "")).startswith("closed")
    )
    source_surface_promotable = False
    proof_status = (
        "closed_conditional_split_implication__source_surface_open"
        if d10_repair_gate_closed
        else "conditional_on_unpromoted_d10_repair_candidate"
    )
    status = "conditional_theorem_only" if d10_repair_gate_closed else "candidate_only"

    return {
        "artifact": "oph_d11_live_exact_split_pair_theorem",
        "generated_utc": _timestamp(),
        "theorem_id": "D11DeclaredSurfaceSplitImplication",
        "legacy_theorem_id": "D11SourceSplitForwardExactness",
        "proof_status": proof_status,
        "status": status,
        "theorem_scope": "declared_d10_d11_running_matching_threshold_surface_only",
        "source_surface_promotable": source_surface_promotable,
        "public_surface_candidate_allowed": False,
        "prediction_promotion_allowed": False,
        "display_allowed_as_conditional": True,
        "upstream_promotion_gate": {
            "required_artifact": "oph_d10_ew_target_free_repair_value_law",
            "required_status": "closed",
            "required_promotion_allowed": True,
            "actual_artifact": d10_repair.get("artifact"),
            "actual_status": d10_repair.get("status"),
            "actual_proof_status": d10_repair.get("proof_status"),
            "actual_promotion_allowed": d10_repair.get("promotion_allowed"),
            "passed": d10_repair_gate_closed,
        },
        "non_circularity_status": {
            "promotion_allowed": False,
            "target_derived_or_candidate_upstream_used": not d10_repair_gate_closed,
            "missing_source_object": (
                "source_emitted_D11_core_Jacobian_residual_selectors_scale_and_pole_packet"
                if d10_repair_gate_closed
                else "closed_promotable_EWTargetFreeRepairValueLaw_D10"
            ),
            "strict_audit_label": (
                "closed_conditional_implication_nonpromoting"
                if d10_repair_gate_closed
                else "conditional_candidate"
            ),
        },
        "source_artifacts": {
            "d10_source_pair": str(D10_SOURCE_JSON),
            "d10_target_free_repair": str(D10_REPAIR_JSON),
            "d11_declared_surface": str(D11_SURFACE_JSON),
            "fixed_ray_no_go": str(D11_NO_GO_JSON),
        },
        "source_tuple": {
            "eta_source": eta_source,
            "beta_EW": beta_ew,
            "lambda_EW": lambda_ew,
            "tau2_tree_exact": tau2_tree_exact,
            "delta_n_tree_exact": delta_n_tree_exact,
        },
        "shared_split_scalar": {
            "symbol": "rho_HT",
            "formula": "log(1 + tau2_tree_exact)",
            "value": rho_ht,
        },
        "split_selector": {
            "A_T": a_t,
            "B_H": b_h,
            "top_residual_formula": "-tau2_tree_exact * eta_source^2 + (1 + beta_EW/28) * eta_source^6 + eta_source^8/14 + eta_source^9/27",
            "top_residual_value": top_residual,
            "higgs_residual_formula": "eta_source^5 - (3/25) * eta_source^6 + lambda_EW * eta_source^6 / 18 + eta_source^8 / (2 * beta_EW)",
            "higgs_residual_value": higgs_residual,
        },
        "declared_surface_exactifier_functions": {
            "c_T_live_formula": "((3/2 + beta_EW/4) * (log(1 + tau2_tree_exact) - tau2_tree_exact) + top_residual) / delta_n_tree_exact",
            "c_T_live_value": c_t_live,
            "c_H_live_formula": "(-(4/3) * tau2_tree_exact + (4/3 - beta_EW/54) * log(1 + tau2_tree_exact) - higgs_residual) / delta_n_tree_exact",
            "c_H_live_value": c_h_live,
        },
        "exact_split_pair": {
            "mH_gev": mh,
            "mt_pole_gev": mt,
            "delta_lambda_mt": delta_lambda,
            "delta_y_t_mt": delta_y,
            "pi_lambda": pi_lambda,
            "pi_y": pi_y,
            "Sigma_HT_exact": sigma_exact,
            "eta_HT_exact": eta_exact,
            "w_HT_exact": pi_y - pi_lambda,
        },
        "readout_formulas": {
            "pi_y": "(eta_source + (3/2 + beta_EW/4) * rho_HT + top_residual) / sqrt(pi)",
            "pi_lambda": "(eta_source - (4/3 - beta_EW/54) * rho_HT + higgs_residual) / sqrt(pi)",
            "rho_HT": "log(1 + tau2_tree_exact)",
            "delta_y_t_mt": "pi_y * y_t_core_mt",
            "delta_lambda_mt": "-(16/9) * pi_lambda * lambda_core_mt",
            "mt_pole_gev": "mt_pole_core_gev + d_mt_pole_d_y_t * delta_y_t_mt",
            "mH_gev": "mH_core_gev + d_mH_d_lambda * delta_lambda_mt",
        },
        "closure_logic": {
            "fixed_ray_central_pair_off_ray": no_go["fixed_ray_point_test"]["central_pair_on_fixed_ray"] is False,
            "fixed_ray_data_compatible_within_one_sigma": bool(
                no_go["fixed_ray_point_test"]["central_pair_within_one_sigma_of_fixed_ray"]
            ),
            "fixed_ray_pull_sigma": float(no_go["uncertainty_pull"]["abs_pull_sigma"]),
            "split_role": "chosen_extension",
            "fixed_ray_no_go_theorem_id": no_go["theorem_id"],
            "smallest_exact_object_above_fixed_ray": "Theta_D11_HT_source(mu_t) = (pi_y, pi_lambda)",
            "equivalent_coordinates": "(Sigma_HT_exact, eta_HT_exact)",
            "fixed_ray_value_formula": "pi_y = pi_lambda = sigma_D11_HT",
            "exact_split_value": pi_y - pi_lambda,
        },
        "proof": [
            (
                "The fixed-ray point statement records that the reference central Higgs/top pair lies off the old one-scalar branch "
                f"(w_HT = {float(no_go['exact_compare_witness']['w_HT_exact']):.5f}, "
                f"{float(no_go['uncertainty_pull']['abs_pull_sigma']):.2f} sigma on the declared linear surface); "
                "the one-scalar ray stays data-compatible and the two-coordinate split is a chosen extension."
            ),
            "The unpromoted D10 repair candidate supplies the tuple (eta_source, beta_EW, lambda_EW, tau2_tree_exact, delta_n_tree_exact) used by this local implication.",
            "The declared D11 map evaluates rho_HT = log(1 + tau2_tree_exact) and the stored top_residual and higgs_residual formulas from that tuple.",
            "Those formulas determine pi_y and pi_lambda without an inverse adapter or direct reference-mass argument at runtime; this does not establish that their historical or physical selector is source-derived.",
            "The declared D11 core and Jacobian then read out delta_y_t_mt, delta_lambda_mt, mt_pole_gev, and mH_gev by direct algebra on the stored surface.",
            "Therefore the artifact closes the algebraic implication on the declared surface only. Full source prediction remains blocked by the D10 selector, D11 provenance/rigidity, scale, RG, and physical-pole gates.",
        ],
        "notes": [
            "This artifact records a conditional local D11 Higgs/top split that is algebraically exact on the declared D10/D11 surface; it is not a full source-only mass prediction.",
            "The D10 quotient-path theorem is presently conditional on QT1-QT5, and the D11 core/Jacobian and residual-selector rigidity are not source-emitted by this artifact.",
            "The stored Higgs and top coordinates are declared-surface readouts, not certified complex propagator poles with uncertainty enclosures.",
            "It does not relabel the old one-scalar fixed ray as exact. The fixed ray remains a lower-rank companion branch beneath this split theorem.",
            "The older target-anchored Higgs-only and top-side exactifier artifacts remain on disk as supporting witness surfaces rather than as the defining live pair theorem.",
            "The quark lane carries only a separate target-audit top coordinate; no repo-wide exact public top row is claimed.",
        ],
        "strictly_not_claimed": [
            "promotion_of_the_old_fixed_ray_as_exact_pair",
            "recovered_core_upgrade_of_the_d11_lane",
            "global_uniqueness_of_the_residual_selector_beyond_the_current_emitted_surface",
            "full_source_only_W_Z_H_complex_pole_prediction",
            "source_derivation_of_the_declared_D11_core_and_Jacobian",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the conditional declared-surface D11 split artifact.")
    parser.add_argument("--d10-source", default=str(D10_SOURCE_JSON))
    parser.add_argument("--d10-repair", default=str(D10_REPAIR_JSON))
    parser.add_argument("--d11-surface", default=str(D11_SURFACE_JSON))
    parser.add_argument("--no-go", default=str(D11_NO_GO_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    d10_source = json.loads(Path(args.d10_source).read_text(encoding="utf-8"))
    d10_repair = json.loads(Path(args.d10_repair).read_text(encoding="utf-8"))
    d11_surface = json.loads(Path(args.d11_surface).read_text(encoding="utf-8"))
    no_go = json.loads(Path(args.no_go).read_text(encoding="utf-8"))
    artifact = canonicalize_artifact_paths(
        build_artifact(d10_source, d10_repair, d11_surface, no_go)
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
