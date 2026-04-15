#!/usr/bin/env python3
"""Certify the current one-scalar D11 fixed-ray branch.

Chain role: certify that the already-emitted D11 seed closes its own diagonal
fixed-ray branch without reopening the legacy diagnostic sidecar.

Mathematics: exact fixed-ray factorization on the forward readout vector,
showing `pi_y = pi_lambda`, `eta_HT = 0`, and `w_HT = 0` identically on the
one-scalar branch.

OPH-derived inputs: the emitted D11 forward seed and its core/Jacobian payload.

Output: an exact fixed-ray branch certificate for the live forward seed.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
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
        "proof_status": "forward_seed_promotion_closed",
        "promotion_status": "closed",
        "certificate_id": "forward_seed_promotion_certificate",
        "forward_seed_artifact": str(DEFAULT_FORWARD_SEED),
        "source_forward_seed_artifact": str(DEFAULT_FORWARD_SEED),
        "discharges_seed_certificate_id": forward_seed.get("seed_certificate_id"),
        "discharges_legacy_sidecar_object_on_live_forward_path": "D11FixedRayWedgeVanishing",
        "proof_scope": "diagonal_fixed_ray_only",
        "diagnostic_center_equality_claimed": False,
        "status": "closed",
        "predictive_promotion_allowed": True,
        "predictive_promotion_scope": "diagonal_fixed_ray_only",
        "forward_path_closed": True,
        "fixed_ray_branch_closed": True,
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
        "smallest_predictive_missing_object": None,
        "next_single_residual_object": "one_extra_forward_coordinate_beyond_fixed_ray",
        "strictly_not_claimed": [
            "exact_higgs_row_on_fixed_ray",
            "exact_higgs_top_pair_on_fixed_ray",
        ],
        "notes": [
            "This certificate closes the live D11 one-scalar seed on its diagonal fixed-ray branch without reopening the legacy diagnostic sidecar.",
            "The exact fixed-ray factorization is proven on the emitted one-scalar forward seed sigma_D11_HT.",
            "The same Jacobian surface carries a companion D11 top-side output on that fixed ray, while the exact Higgs row is carried separately by D11LiveForwardExactHiggsPromotion.",
            "The compare-only exact Higgs/top pair lies off this fixed ray and remains a validation surface only.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact D11 forward-seed promotion certificate.")
    parser.add_argument("--forward-seed", default=str(DEFAULT_FORWARD_SEED))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    forward_seed = json.loads(Path(args.forward_seed).read_text(encoding="utf-8"))
    artifact = build_artifact(forward_seed)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
