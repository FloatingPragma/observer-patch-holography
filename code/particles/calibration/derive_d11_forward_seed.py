#!/usr/bin/env python3
"""Propagate the D10 gauge core into the compact D11 Higgs/top seed.

Chain role: reuse the calibrated D10 gauge carrier to produce the current
forward seed for the Higgs and top branch.

Mathematics: a one-scalar perturbation of the D11 core followed by linear
Jacobian readout for `m_t` and `m_H`.

OPH-derived inputs: `alpha_u`, `alphaY_mz`, and `alpha2_mz` from the D10 core,
plus the D11 critical-surface core and Jacobian payload.

Output: the public `/particles` Higgs/top candidate and the remaining D11
promotion certificate gap.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "calibration" / "d11_critical_surface_readout.json"
DEFAULT_D10_SOURCE = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"
DEFAULT_PROMOTION_CERTIFICATE = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed_promotion_certificate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(payload: dict, d10_source: dict, promotion_certificate: dict | None = None) -> dict:
    core = dict(payload["core"])
    jacobian = dict(payload["jacobian"])
    d10_core = dict(d10_source["core_source"])
    alpha_y = 3.0 * float(d10_core["alpha1_mz"]) / 5.0
    alpha2 = float(d10_core["alpha2_mz"])
    alpha_u = float(d10_core["alpha_u"])
    cos2_theta_w0 = (alpha2 - alpha_y) / (alpha2 + alpha_y)
    sigma = alpha_u * cos2_theta_w0 / math.sqrt(math.pi)
    kappa = 16.0 / 9.0
    delta_y = sigma * float(core["y_t_core_mt"])
    delta_lambda = -kappa * sigma * float(core["lambda_core_mt"])
    mt_pole = float(core["mt_pole_core_gev"]) + float(jacobian["d_mt_pole_d_y_t"]) * delta_y
    m_h = float(core["mH_core_gev"]) + float(jacobian["d_mH_d_lambda"]) * delta_lambda

    promotion_certificate = promotion_certificate or {}
    promotion_closed = promotion_certificate.get("status") == "closed"

    return {
        "artifact": "oph_d11_forward_seed",
        "generated_utc": _timestamp(),
        "forward_seed_object": "sigma_D11_HT",
        "seed_status": "promoted_exact_forward_seed" if promotion_closed else "predictive_forward_law_candidate",
        "predictive_promotion_allowed": promotion_closed,
        "public_surface_candidate_allowed": True,
        "public_surface_candidate_scope": ["mH_gev", "mt_pole_gev"],
        "public_surface_policy": "best_available_reference_free_forward_seed_candidate",
        "smallest_predictive_missing_object": None if promotion_closed else "forward_seed_promotion_certificate",
        "closure_state": "forward_closed_exact_on_live_seed" if promotion_closed else "forward_seed_open_waiting_certificate",
        "forward_seed_promotion_certificate": str(DEFAULT_PROMOTION_CERTIFICATE) if promotion_closed else None,
        "promotion_certificate_artifact": str(DEFAULT_PROMOTION_CERTIFICATE) if promotion_closed else None,
        "source_artifact": payload.get("artifact"),
        "source_predictive_status": "reference_free_compact_seed_law",
        "kappa_HT": kappa,
        "sigma_D11_HT": sigma,
        "source_seed_law": {
            "sigma_D11_HT_formula": "alpha_u * cos(2*theta_W0) / sqrt(pi)",
            "cos2_theta_W0_formula": "(alpha2_mz - alphaY_mz) / (alpha2_mz + alphaY_mz)",
            "alphaY_mz": alpha_y,
            "alpha2_mz": alpha2,
            "alpha_u": alpha_u,
            "cos2_theta_W0": cos2_theta_w0,
        },
        "theta_from_seed": {
            "delta_y_t_mt": delta_y,
            "delta_lambda_mt": delta_lambda,
        },
        "seed_equality_law": "delta_y_t_mt / y_t_core_mt = -(9/16) * delta_lambda_mt / lambda_core_mt = sigma_D11_HT",
        "core": {
            "y_t_core_mt": float(core["y_t_core_mt"]),
            "lambda_core_mt": float(core["lambda_core_mt"]),
            "mt_pole_core_gev": float(core["mt_pole_core_gev"]),
            "mH_core_gev": float(core["mH_core_gev"]),
        },
        "mass_readout": {
            "mt_pole_gev": mt_pole,
            "mH_gev": m_h,
            "mt_pole_formula": "mt_pole_core_gev + d_mt_pole_d_y_t * delta_y_t_mt",
            "mH_formula": "mH_core_gev + d_mH_d_lambda * delta_lambda_mt",
            "jacobian": {
                "d_mt_pole_d_y_t": float(jacobian["d_mt_pole_d_y_t"]),
                "d_mH_d_lambda": float(jacobian["d_mH_d_lambda"]),
            },
        },
        "seed_certificate_id": "D11SeedCommonProvenanceCertificate",
        "diagnostic_literals_forbidden": True,
        "compare_only_blocks": [
            "diagnostic_required_to_current_refs",
            "legacy_diagnostic_readback",
            "center_witness_forward_consequence",
        ],
        "notes": [
            "Compact D11 forward seed extracted from the current critical-surface artifact.",
            "This artifact now emits one reference-free compact seed law from the D10 gauge core: sigma_D11_HT = alpha_u * cos(2*theta_W0) / sqrt(pi).",
            "The resulting one-scalar chain sigma_D11_HT -> Theta_D11_HT -> (m_t, m_H) is a forward candidate and no longer depends on the diagnostic delta_y_t / delta_lambda readback values.",
            (
                "The live forward D11 path is closed by the forward-seed promotion certificate on the emitted one-scalar seed."
                if promotion_closed
                else "No larger Higgs/top family is missing here; the remaining work is promotion/certification of this one compact seed, not a wider Higgs/top residual program."
            ),
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the compact D11 forward-seed artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--d10-source", default=str(DEFAULT_D10_SOURCE))
    parser.add_argument("--promotion-certificate", default=str(DEFAULT_PROMOTION_CERTIFICATE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    d10_source = json.loads(Path(args.d10_source).read_text(encoding="utf-8"))
    promotion_certificate_path = Path(args.promotion_certificate)
    promotion_certificate = (
        json.loads(promotion_certificate_path.read_text(encoding="utf-8"))
        if promotion_certificate_path.exists()
        else None
    )
    artifact = build_artifact(payload, d10_source, promotion_certificate)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
