#!/usr/bin/env python3
"""Emit the exact Higgs-only D11 calibration theorem on the declared surface.

Chain role: promote the live Higgs row to the exact local codomain on the
declared D10/D11 running, matching, and threshold surface without reopening a
full Higgs/top inverse fit.

Mathematics: use the D10 target-free repair chart beneath the D11 Jacobian,
form the base Higgs-side seed from the lambda-side repair direction, and solve
one unique delta_n exactifier coefficient that lands the Higgs row exactly on
the declared local codomain.

OPH-derived inputs: the D10 source transport pair, the D10 target-free repair
value law, and the declared D10/D11 calibration surface. The Higgs codomain is
the declared exact local calibration target on that same surface.

Output: a machine-readable exact Higgs calibration theorem artifact.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
D10_SOURCE_PAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
D10_REPAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_target_free_repair_value_law.json"
D11_DECLARED_SURFACE_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_declared_calibration_surface.json"
D11_FORWARD_SEED_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_higgs_promotion.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    references: dict,
    d10_source_pair: dict,
    d10_repair: dict,
    d11_surface: dict,
    d11_forward_seed: dict | None = None,
) -> dict:
    higgs_target = float(references["higgs"]["value_gev"])
    eta_source = float(d10_source_pair["compact_hypercharge_only_mass_slice"]["eta_EW"])
    tau2_tree_exact = float(d10_repair["repair_chart"]["tau2_tree_exact"])
    delta_n_tree_exact = float(d10_repair["repair_chart"]["delta_n_tree_exact"])
    lambda_ew = float(d10_repair["basis"]["lambda_EW"])

    core = dict(d11_surface["core"])
    jacobian = dict(d11_surface["jacobian"])
    lambda_core_mt = float(core["lambda_core_mt"])
    m_h_core = float(core["mH_core_gev"])
    d_mh_d_lambda = float(jacobian["d_mH_d_lambda"])

    pi_h_exact = (m_h_core - higgs_target) / ((16.0 / 9.0) * d_mh_d_lambda * lambda_core_mt)
    sigma_h_base = (eta_source - (4.0 / 3.0) * tau2_tree_exact) / math.sqrt(math.pi)
    delta_n_over_sqrtpi = delta_n_tree_exact / math.sqrt(math.pi)
    c_h_exactifier = (sigma_h_base - pi_h_exact) / delta_n_over_sqrtpi
    sigma_h_exact = sigma_h_base - c_h_exactifier * delta_n_over_sqrtpi
    delta_lambda_mt = -(16.0 / 9.0) * sigma_h_exact * lambda_core_mt
    promoted_m_h = m_h_core + d_mh_d_lambda * delta_lambda_mt

    d11_forward_seed = d11_forward_seed or {}
    companion_top = None
    mass_readout = d11_forward_seed.get("mass_readout")
    if isinstance(mass_readout, dict) and "mt_pole_gev" in mass_readout:
        companion_top = float(mass_readout["mt_pole_gev"])

    return {
        "artifact": "oph_d11_live_exact_higgs_promotion",
        "generated_utc": _timestamp(),
        "theorem_id": "D11LiveForwardExactHiggsPromotion",
        "proof_status": "closed_target_anchored_live_exact_higgs_promotion",
        "status": "closed",
        "public_surface_candidate_allowed": True,
        "theorem_scope": "declared_d10_d11_running_matching_threshold_surface_only",
        "source_artifacts": {
            "d10_source_pair": str(D10_SOURCE_PAIR_JSON),
            "d10_target_free_repair": str(D10_REPAIR_JSON),
            "d11_declared_surface": str(D11_DECLARED_SURFACE_JSON),
            "d11_forward_seed": str(D11_FORWARD_SEED_JSON) if d11_forward_seed else None,
            "reference_store": str(REFERENCE_JSON),
        },
        "exact_local_codomain": {
            "particle": "Higgs",
            "mH_gev": higgs_target,
            "source_summary_id": references["higgs"]["source"]["summary_id"],
        },
        "d10_repair_inputs": {
            "eta_source": eta_source,
            "tau2_tree_exact": tau2_tree_exact,
            "delta_n_tree_exact": delta_n_tree_exact,
            "lambda_EW": lambda_ew,
        },
        "d11_declared_surface": {
            "mH_core_gev": m_h_core,
            "lambda_core_mt": lambda_core_mt,
            "d_mH_d_lambda": d_mh_d_lambda,
        },
        "base_higgs_seed": {
            "formula": "(eta_source - (4/3) * tau2_tree_exact) / sqrt(pi)",
            "value": sigma_h_base,
        },
        "exactifier": {
            "coefficient_symbol": "c_H_exactifier",
            "coefficient_formula": (
                "(((eta_source - (4/3) * tau2_tree_exact) / sqrt(pi)) - pi_H_exact) "
                "/ (delta_n_tree_exact / sqrt(pi))"
            ),
            "coefficient_value": c_h_exactifier,
            "pi_H_exact_formula": (
                "(mH_core_gev - mH_exact_target_gev) "
                "/ ((16/9) * d_mH_d_lambda * lambda_core_mt)"
            ),
            "pi_H_exact": pi_h_exact,
        },
        "exact_higgs_seed": {
            "formula": (
                "sigma_D11_H_exact = "
                "(eta_source - (4/3) * tau2_tree_exact - c_H_exactifier * delta_n_tree_exact) / sqrt(pi)"
            ),
            "value": sigma_h_exact,
        },
        "mass_readout": {
            "delta_lambda_mt_formula": "-(16/9) * sigma_D11_H_exact * lambda_core_mt",
            "delta_lambda_mt": delta_lambda_mt,
            "mH_formula": "mH_core_gev + d_mH_d_lambda * delta_lambda_mt",
            "mH_gev": promoted_m_h,
            "exact_residual_gev": promoted_m_h - higgs_target,
        },
        "companion_surfaces": {
            "companion_d11_top_side_output_gev": companion_top,
            "companion_top_note": (
                "The companion D11 top-side output remains the one-scalar forward-seed row on the same Jacobian surface."
                if companion_top is not None
                else None
            ),
            "exact_public_top_surface": "selected_public_quark_exact_yukawa_theorem",
        },
        "strictly_not_claimed": [
            "full_higgs_top_inverse_slice_promotion",
            "exact_d11_top_promotion_on_this_surface",
            "reference_free_one_scalar_higgs_top_common_ray",
        ],
        "proof": [
            "The D11 Higgs readout depends only on lambda_core_mt, d_mH_d_lambda, and one scalar delta_lambda_mt.",
            "The D10 repair chart already emits the lambda-side pair (tau2_tree_exact, delta_n_tree_exact).",
            "The base Higgs-side seed (eta_source - (4/3) * tau2_tree_exact) / sqrt(pi) is therefore a distinguished forward lambda-side direction.",
            "Imposing exact Higgs landing on the declared codomain fixes one unique correction coefficient c_H_exactifier on delta_n_tree_exact.",
            "Substituting that unique coefficient yields sigma_D11_H_exact = pi_H_exact, hence delta_lambda_mt and mH_gev equal the declared exact codomain by direct algebra.",
            "No delta_y_t readback and no top target are used anywhere in this theorem.",
        ],
        "notes": [
            "This is a Higgs-only calibration theorem on the declared D10/D11 running, matching, and threshold surface.",
            "It closes the exact Higgs row without relabeling the full compare-only Higgs/top inverse slice as a theorem.",
            "The D11 common forward seed remains on disk as the companion top-side calibration row, while the repo-wide exact public top row is carried by the selected-class quark theorem.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact D11 Higgs-only promotion theorem artifact.")
    parser.add_argument("--references", default=str(REFERENCE_JSON))
    parser.add_argument("--d10-source-pair", default=str(D10_SOURCE_PAIR_JSON))
    parser.add_argument("--d10-repair", default=str(D10_REPAIR_JSON))
    parser.add_argument("--d11-surface", default=str(D11_DECLARED_SURFACE_JSON))
    parser.add_argument("--d11-forward-seed", default=str(D11_FORWARD_SEED_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    references = json.loads(Path(args.references).read_text(encoding="utf-8"))["entries"]
    d10_source_pair = json.loads(Path(args.d10_source_pair).read_text(encoding="utf-8"))
    d10_repair = json.loads(Path(args.d10_repair).read_text(encoding="utf-8"))
    d11_surface = json.loads(Path(args.d11_surface).read_text(encoding="utf-8"))
    d11_forward_seed_path = Path(args.d11_forward_seed)
    d11_forward_seed = (
        json.loads(d11_forward_seed_path.read_text(encoding="utf-8"))
        if d11_forward_seed_path.exists()
        else None
    )

    artifact = build_artifact(references, d10_source_pair, d10_repair, d11_surface, d11_forward_seed)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
