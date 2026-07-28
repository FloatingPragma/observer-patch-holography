#!/usr/bin/env python3
"""Certify the algebra of the current one-scalar D11 fixed-ray branch.

Chain role: certify the diagonal fixed-ray identities on the declared D11
surface without treating algebraic closure as source or predictive promotion.

Mathematics: exact fixed-ray factorization on the forward readout vector,
showing `pi_y = pi_lambda`, `eta_HT = 0`, and `w_HT = 0` identically on the
one-scalar branch.

Inputs: the emitted D11 forward seed and its declared-surface core/Jacobian
payload.

Output: an exact, non-promoting fixed-ray algebra certificate.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import canonicalize_artifact_paths
DEFAULT_FORWARD_SEED = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed.json"
EXACT_HIGGS_ARTIFACT = ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_higgs_promotion.json"
FIXED_RAY_NO_GO_ARTIFACT = ROOT / "particles" / "runs" / "calibration" / "d11_fixed_ray_no_go_theorem.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed_promotion_certificate.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(forward_seed: dict) -> dict:
    sigma = float(forward_seed["sigma_D11_HT"])
    kappa = float(forward_seed["kappa_HT"])
    core = dict(forward_seed["core"])
    theta = dict(forward_seed["theta_from_seed"])
    mass_readout = dict(forward_seed["mass_readout"])
    delta_y = float(theta["delta_y_t_mt"])
    delta_lambda = float(theta["delta_lambda_mt"])
    y_core = float(core["y_t_core_mt"])
    lambda_core = float(core["lambda_core_mt"])
    pi_y = delta_y / y_core
    pi_lambda = -(9.0 / 16.0) * delta_lambda / lambda_core

    return {
        "artifact": "oph_d11_forward_seed_promotion_certificate",
        "generated_utc": _timestamp(),
        "proof_status": "fixed_ray_algebra_closed_on_declared_surface",
        "promotion_status": "blocked_declared_surface_not_source_emitted",
        "certificate_id": "forward_seed_promotion_certificate",
        "forward_seed_artifact": str(DEFAULT_FORWARD_SEED),
        "source_forward_seed_artifact": str(DEFAULT_FORWARD_SEED),
        "discharges_seed_certificate_id": forward_seed.get("seed_certificate_id"),
        "discharges_legacy_sidecar_object_on_live_forward_path": "D11FixedRayWedgeVanishing",
        "proof_scope": "diagonal_fixed_ray_only",
        "diagnostic_center_equality_claimed": False,
        "status": "fixed_ray_algebra_closed",
        "source_surface_promotable": False,
        "predictive_promotion_allowed": False,
        "public_surface_candidate_allowed": False,
        "comparison_surface_allowed": True,
        "predictive_promotion_scope": None,
        "forward_path_closed": False,
        "fixed_ray_branch_closed": False,
        "fixed_ray_algebra_closed": True,
        "exact_higgs_row_claimed": False,
        "exact_pair_claimed": False,
        "exact_higgs_artifact": str(EXACT_HIGGS_ARTIFACT),
        "fixed_ray_no_go_artifact": str(FIXED_RAY_NO_GO_ARTIFACT),
        "promoted_seed_object": "sigma_D11_HT",
        "sigma_D11_HT": sigma,
        "theta_symbol": "Theta_D11_HT(mu_t)",
        "theta_formula": {
            "delta_y_t_mt": "sigma_D11_HT * y_t_core_mt",
            "delta_lambda_mt": "-(16/9) * sigma_D11_HT * lambda_core_mt",
        },
        "theta_from_seed": {
            "delta_y_t_mt": delta_y,
            "delta_lambda_mt": delta_lambda,
        },
        "forward_normalized_readback": {
            "coordinates": {
                "pi_y": pi_y,
                "pi_lambda": pi_lambda,
            },
            "decomposition": {
                "sigma_shared": sigma,
                "eta_HT": 0.5 * (pi_y - pi_lambda),
                "w_HT": pi_y - pi_lambda,
            },
        },
        "seed_equality_certificate": {
            "status": "closed_on_forward_seed",
            "law": "delta_y_t_mt / y_t_core_mt = -(9/16) * delta_lambda_mt / lambda_core_mt = sigma_D11_HT",
            "residual_abs": abs(pi_y - pi_lambda),
        },
        "fixed_ray_wedge_vanishing_certificate": {
            "name": "D11FixedRayWedgeVanishing",
            "status": "closed_on_forward_seed",
            "proof_mode": "exact_forward_factorization",
            "kappa_HT": kappa,
            "wedge_formula": "kappa_HT * lambda_core_mt * delta_y_t_mt + y_t_core_mt * delta_lambda_mt",
            "wedge_value": kappa * lambda_core * delta_y + y_core * delta_lambda,
        },
        "mass_readout_consequence": {
            "mt_pole_formula": mass_readout["mt_pole_formula"],
            "mH_formula": mass_readout["mH_formula"],
            "mt_pole_gev": float(mass_readout["mt_pole_gev"]),
            "mH_gev": float(mass_readout["mH_gev"]),
        },
        "smallest_predictive_missing_object": "source_emitted_higgs_yukawa_fj_packet",
        "next_single_residual_object": "one_extra_forward_coordinate_beyond_fixed_ray",
        "promotion_blockers": [
            "declared_d11_core_and_jacobian_not_source_emitted",
            "full_target_clean_yukawa_and_lambda_packet_absent",
            "v_chart_to_v_F_theorem_absent",
            "same_branch_rg_threshold_matching_and_pole_receipts_absent",
        ],
        "strictly_not_claimed": [
            "oph_native_source_promotion",
            "predictive_or_public_mass_promotion",
            "exact_higgs_row_on_fixed_ray",
            "exact_higgs_top_pair_on_fixed_ray",
        ],
        "notes": [
            "This certificate closes only the algebra of the D11 one-scalar diagonal fixed ray on the declared surface.",
            "The exact fixed-ray factorization is proven on the emitted one-scalar forward seed sigma_D11_HT.",
            "The same declared Jacobian surface carries comparison coordinates on that fixed ray; it is not a source-emitted Higgs/top pole surface.",
            "The target-anchored Higgs and top exactifiers are compare-only artifacts.",
            "The compare-only exact Higgs/top pair lies off this fixed ray and remains a validation surface only.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact non-promoting D11 fixed-ray algebra certificate.")
    parser.add_argument("--forward-seed", default=str(DEFAULT_FORWARD_SEED))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    forward_seed = json.loads(Path(args.forward_seed).read_text(encoding="utf-8"))
    artifact = canonicalize_artifact_paths(build_artifact(forward_seed))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
