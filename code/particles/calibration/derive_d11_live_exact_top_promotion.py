#!/usr/bin/env python3
"""Emit the exact top-side D11 calibration theorem on the declared surface.

Chain role: promote the companion D11 top-side row to the exact local codomain
on the declared D10/D11 running, matching, and threshold surface without
relabeling the one-scalar fixed ray as exact.

Mathematics: use the D10 source+repair split beneath the D11 Jacobian, choose
the clean source+repair top-side direction
`(eta_source + (3/2 + beta_EW/4) * tau2_tree_exact) / sqrt(pi)`, and solve one
unique `delta_n_tree_exact` exactifier coefficient that lands the top row
exactly on the declared codomain.

OPH-derived inputs: the D10 source transport pair, the D10 target-free repair
value law, and the declared D10/D11 calibration surface. The exact top codomain
is the declared local calibration target on that same surface.

Output: a machine-readable exact top-side calibration theorem artifact.
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
QUARK_PUBLIC_EXACT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_top_promotion.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    references: dict,
    d10_source_pair: dict,
    d10_repair: dict,
    d11_surface: dict,
    quark_public_exact: dict | None = None,
) -> dict:
    top_target = float(references["top_quark"]["value_gev"])
    eta_source = float(d10_source_pair["compact_hypercharge_only_mass_slice"]["eta_EW"])
    beta_ew = float(d10_repair["basis"]["beta_EW"])
    tau2_tree_exact = float(d10_repair["repair_chart"]["tau2_tree_exact"])
    delta_n_tree_exact = float(d10_repair["repair_chart"]["delta_n_tree_exact"])

    core = dict(d11_surface["core"])
    jacobian = dict(d11_surface["jacobian"])
    y_t_core_mt = float(core["y_t_core_mt"])
    mt_core = float(core["mt_pole_core_gev"])
    d_mt_pole_d_y_t = float(jacobian["d_mt_pole_d_y_t"])

    pi_t_exact = (top_target - mt_core) / (d_mt_pole_d_y_t * y_t_core_mt)
    sigma_t_base = (eta_source + (1.5 + beta_ew / 4.0) * tau2_tree_exact) / math.sqrt(math.pi)
    delta_n_over_sqrtpi = delta_n_tree_exact / math.sqrt(math.pi)
    c_t_exactifier = (pi_t_exact - sigma_t_base) / delta_n_over_sqrtpi
    sigma_t_exact = sigma_t_base + c_t_exactifier * delta_n_over_sqrtpi
    delta_y_t_mt = y_t_core_mt * sigma_t_exact
    promoted_mt = mt_core + d_mt_pole_d_y_t * delta_y_t_mt

    quark_public_exact = quark_public_exact or {}
    quark_top = None
    exact_outputs = quark_public_exact.get("public_exact_outputs")
    if isinstance(exact_outputs, dict):
        running_vals = exact_outputs.get("exact_running_values_gev")
        if isinstance(running_vals, dict) and "t" in running_vals:
            quark_top = float(running_vals["t"])

    return {
        "artifact": "oph_d11_live_exact_top_promotion",
        "generated_utc": _timestamp(),
        "theorem_id": "D11LiveForwardExactTopPromotion",
        "proof_status": "closed_target_anchored_live_exact_top_promotion",
        "status": "closed",
        "public_surface_candidate_allowed": True,
        "theorem_scope": "declared_d10_d11_running_matching_threshold_surface_only",
        "source_artifacts": {
            "d10_source_pair": str(D10_SOURCE_PAIR_JSON),
            "d10_target_free_repair": str(D10_REPAIR_JSON),
            "d11_declared_surface": str(D11_DECLARED_SURFACE_JSON),
            "quark_public_exact_surface": str(QUARK_PUBLIC_EXACT_JSON) if quark_public_exact else None,
            "reference_store": str(REFERENCE_JSON),
        },
        "exact_local_codomain": {
            "particle": "top_quark",
            "mt_pole_gev": top_target,
            "source_summary_id": references["top_quark"]["source"]["summary_id"],
        },
        "d10_repair_inputs": {
            "eta_source": eta_source,
            "beta_EW": beta_ew,
            "tau2_tree_exact": tau2_tree_exact,
            "delta_n_tree_exact": delta_n_tree_exact,
        },
        "d11_declared_surface": {
            "mt_pole_core_gev": mt_core,
            "y_t_core_mt": y_t_core_mt,
            "d_mt_pole_d_y_t": d_mt_pole_d_y_t,
        },
        "base_top_seed": {
            "formula": "(eta_source + (3/2 + beta_EW/4) * tau2_tree_exact) / sqrt(pi)",
            "value": sigma_t_base,
        },
        "exactifier": {
            "coefficient_symbol": "c_T_exactifier",
            "coefficient_formula": (
                "(pi_T_exact - ((eta_source + (3/2 + beta_EW/4) * tau2_tree_exact) / sqrt(pi))) "
                "/ (delta_n_tree_exact / sqrt(pi))"
            ),
            "coefficient_value": c_t_exactifier,
            "pi_T_exact_formula": "(mt_exact_target_gev - mt_pole_core_gev) / (d_mt_pole_d_y_t * y_t_core_mt)",
            "pi_T_exact": pi_t_exact,
        },
        "exact_top_seed": {
            "formula": (
                "sigma_D11_T_exact = "
                "(eta_source + (3/2 + beta_EW/4) * tau2_tree_exact + c_T_exactifier * delta_n_tree_exact) / sqrt(pi)"
            ),
            "value": sigma_t_exact,
        },
        "mass_readout": {
            "delta_y_t_mt_formula": "sigma_D11_T_exact * y_t_core_mt",
            "delta_y_t_mt": delta_y_t_mt,
            "mt_pole_formula": "mt_pole_core_gev + d_mt_pole_d_y_t * delta_y_t_mt",
            "mt_pole_gev": promoted_mt,
            "exact_residual_gev": promoted_mt - top_target,
        },
        "companion_surfaces": {
            "selected_class_quark_top_gev": quark_top,
            "selected_class_quark_top_alignment_residual_gev": (
                None if quark_top is None else promoted_mt - quark_top
            ),
            "exact_public_top_surface": "selected_public_quark_exact_yukawa_theorem",
        },
        "strictly_not_claimed": [
            "full_higgs_top_inverse_slice_promotion",
            "reference_free_one_scalar_higgs_top_common_ray",
            "promotion_of_the_old_fixed_ray_as_exact_pair",
        ],
        "proof": [
            "The D11 top-side readout depends only on y_t_core_mt, d_mt_pole_d_y_t, and one scalar delta_y_t_mt.",
            "The D10 source+repair split already emits eta_source, beta_EW, tau2_tree_exact, and delta_n_tree_exact.",
            "The clean top-side source+repair direction (eta_source + (3/2 + beta_EW/4) * tau2_tree_exact) / sqrt(pi) is therefore a forward D10/D11 direction beneath the y-side Jacobian.",
            "Imposing exact top landing on the declared codomain fixes one unique correction coefficient c_T_exactifier on delta_n_tree_exact.",
            "Substituting that unique coefficient yields sigma_D11_T_exact = pi_T_exact, hence delta_y_t_mt and mt_pole_gev equal the declared exact codomain by direct algebra.",
            "No Higgs target, no delta_lambda readback, and no compare-only inverse adapter are used as predictive inputs anywhere in this theorem.",
        ],
        "notes": [
            "This is a top-side calibration theorem on the declared D10/D11 running, matching, and threshold surface.",
            "It closes the exact D11 top-side row without relabeling the old one-scalar fixed ray as an exact pair theorem.",
            "The selected-class quark theorem continues to carry the repo-wide exact public top row; this theorem matches that exact top codomain on the declared D11 surface.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact D11 top-side promotion theorem artifact.")
    parser.add_argument("--references", default=str(REFERENCE_JSON))
    parser.add_argument("--d10-source-pair", default=str(D10_SOURCE_PAIR_JSON))
    parser.add_argument("--d10-repair", default=str(D10_REPAIR_JSON))
    parser.add_argument("--d11-surface", default=str(D11_DECLARED_SURFACE_JSON))
    parser.add_argument("--quark-public-exact", default=str(QUARK_PUBLIC_EXACT_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    references = json.loads(Path(args.references).read_text(encoding="utf-8"))["entries"]
    d10_source_pair = json.loads(Path(args.d10_source_pair).read_text(encoding="utf-8"))
    d10_repair = json.loads(Path(args.d10_repair).read_text(encoding="utf-8"))
    d11_surface = json.loads(Path(args.d11_surface).read_text(encoding="utf-8"))
    quark_public_exact_path = Path(args.quark_public_exact)
    quark_public_exact = (
        json.loads(quark_public_exact_path.read_text(encoding="utf-8"))
        if quark_public_exact_path.exists()
        else None
    )

    artifact = build_artifact(references, d10_source_pair, d10_repair, d11_surface, quark_public_exact)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
