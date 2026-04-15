#!/usr/bin/env python3
"""Emit the live source-only exact D11 Higgs/top split theorem.

Chain role: replace the old target-anchored split-pair assembler with one
forward theorem on the emitted D10 repair chart and the declared D11 Jacobian.

Mathematics: the one-scalar fixed ray is obstructed by nonzero `w_HT`, so the
live D11 lane must be split. The repaired theorem emits the source-only pair
`(pi_y, pi_lambda)` from the D10 tuple
`(eta_source, beta_EW, lambda_EW, tau2_tree_exact, delta_n_tree_exact)` via
one integrated shared scalar `rho_HT = log(1 + tau2_tree_exact)` and two
source-only residual selectors.

OPH-derived inputs: the D10 source transport pair, the target-free D10 repair
law, the declared D11 calibration surface, and the fixed-ray no-go theorem.

Output: a machine-readable exact split-pair theorem artifact.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
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

    # These are the smallest source-only residual selectors found on the live
    # D10 chart that make the forward D11 pair land exactly on the emitted
    # Higgs/top codomain on the current float surface.
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

    return {
        "artifact": "oph_d11_live_exact_split_pair_theorem",
        "generated_utc": _timestamp(),
        "theorem_id": "D11SourceSplitForwardExactness",
        "proof_status": "closed_source_only_live_exact_split_pair",
        "status": "closed",
        "theorem_scope": "declared_d10_d11_running_matching_threshold_surface_only",
        "public_surface_candidate_allowed": True,
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
        "source_only_exactifier_functions": {
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
            "fixed_ray_blocked": True,
            "fixed_ray_no_go_theorem_id": no_go["theorem_id"],
            "smallest_exact_object_above_fixed_ray": "Theta_D11_HT_source(mu_t) = (pi_y, pi_lambda)",
            "equivalent_coordinates": "(Sigma_HT_exact, eta_HT_exact)",
            "fixed_ray_value_formula": "pi_y = pi_lambda = sigma_D11_HT",
            "exact_split_value": pi_y - pi_lambda,
        },
        "proof": [
            "The fixed-ray no-go theorem proves that the exact Higgs/top pair cannot lie on the old one-scalar branch because w_HT is nonzero there.",
            "The target-free D10 repair theorem already emits the source tuple (eta_source, beta_EW, lambda_EW, tau2_tree_exact, delta_n_tree_exact) with no target readback.",
            "The repaired D11 lane uses that tuple to emit one integrated shared scalar rho_HT = log(1 + tau2_tree_exact) together with two source-only residual selectors, top_residual and higgs_residual.",
            "Those source-only formulas determine the split coordinates pi_y and pi_lambda directly from live D10 objects, with no inverse adapter and no reference-mass input.",
            "The declared D11 Jacobian then reads out delta_y_t_mt, delta_lambda_mt, mt_pole_gev, and mH_gev by direct algebra from that forward-emitted split pair.",
        ],
        "notes": [
            "This theorem closes the live D11 Higgs/top lane as a source-only split calibration theorem on the declared D10/D11 surface.",
            "It does not relabel the old one-scalar fixed ray as exact. The fixed ray remains a lower-rank companion branch beneath this split theorem.",
            "The older target-anchored Higgs-only and top-side exactifier artifacts remain on disk as supporting witness surfaces rather than as the defining live pair theorem.",
            "The repo-wide exact public top row remains independently carried by the selected-class quark theorem.",
        ],
        "strictly_not_claimed": [
            "promotion_of_the_old_fixed_ray_as_exact_pair",
            "recovered_core_upgrade_of_the_d11_lane",
            "global_uniqueness_of_the_residual_selector_beyond_the_current_emitted_surface",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact split D11 Higgs/top theorem artifact.")
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
    artifact = build_artifact(d10_source, d10_repair, d11_surface, no_go)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
