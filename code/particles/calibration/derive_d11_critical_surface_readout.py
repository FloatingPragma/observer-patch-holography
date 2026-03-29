#!/usr/bin/env python3
"""Export the current D11 common readout-kernel boundary artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d11_critical_surface_readout.json"
DEFAULT_D10_SOURCE = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"
DEFAULT_RESULTS_STATUS = ROOT / "particles" / "RESULTS_STATUS.md"
DEFAULT_FORWARD_SEED_CERTIFICATE = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed_promotion_certificate.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(d10_source: Path, results_status: Path, forward_seed_certificate: dict | None = None) -> dict[str, object]:
    y_t_core = 0.92046435
    lambda_core = 0.13164915
    delta_y = 0.01156247
    delta_lambda = -0.00294158
    kappa_candidate = 16.0 / 9.0
    sigma_from_y = delta_y / y_t_core
    sigma_from_lambda = -delta_lambda / (kappa_candidate * lambda_core)
    sigma_shared = 0.5 * (sigma_from_y + sigma_from_lambda)
    eta_ht = 0.5 * (sigma_from_y - sigma_from_lambda)
    w_ht = 2.0 * eta_ht
    omega_ht = (kappa_candidate * lambda_core * delta_y) + (y_t_core * delta_lambda)
    kappa_observed = (-delta_lambda / lambda_core) / sigma_from_y
    shared_delta_y = sigma_shared * y_t_core
    shared_delta_lambda = -kappa_candidate * sigma_shared * lambda_core
    kappa_residual = kappa_observed - kappa_candidate
    prefactor = lambda_core * delta_y
    forward_seed_certificate = forward_seed_certificate or {}
    forward_closed = forward_seed_certificate.get("status") == "closed"

    return {
        "artifact": "oph_d11_critical_surface_readout",
        "generated_utc": _timestamp(),
        "predictive_status": (
            "diagnostic_sidecar_only__live_forward_path_closed_elsewhere"
            if forward_closed
            else "diagnostic_only_until_forward_seed_artifact_exists"
        ),
        "predictive_promotion_allowed": False,
        "d10_source_artifact": str(d10_source),
        "results_status_surface": str(results_status),
        "transport_family": "sm_below_sync__mssm_like_gauge_above_sync",
        "mu_sync_gev": 6.8e11,
        "mu_sync_status": "candidate",
        "mu_eval_gev": 160.61247,
        "core": {
            "mt_ms_gev": 160.61247,
            "mt_pole_core_gev": 170.26125,
            "mH_core_gev": 126.62263,
            "y_t_core_mt": y_t_core,
            "lambda_core_mt": lambda_core,
            "alpha_s_mt": 0.11018777,
            "pole_ratio_core": 1.06007492,
        },
        "readout_kernel": {
            "name": "CriticalSurfaceReadoutKernel_D11",
            "status": "collapsed_common_seed_constructive__strict_antidiagonal_vanishing_open",
            "family": "relative_core_ray",
            "seed_symbol": "sigma_D11_HT(mu_t)",
            "seed_value": sigma_shared,
            "theorem_candidate": "D11SeedCommonProvenanceCertificate",
            "kappa_lambda_over_y": {
                "candidate": kappa_candidate,
                "observed_diagnostic": kappa_observed,
            },
            "delta_y_t_formula": "sigma_D11_HT * y_t_core_mt",
            "delta_lambda_formula": "-kappa_lambda_over_y.candidate * sigma_D11_HT * lambda_core_mt",
            "sigma_from_y_t": sigma_from_y,
            "sigma_from_lambda": sigma_from_lambda,
            "sigma_shared": sigma_shared,
            "eta_HT": eta_ht,
            "w_HT": w_ht,
            "shared_delta_y": shared_delta_y,
            "shared_delta_lambda": shared_delta_lambda,
            "seed_equality_law": "sigma_y = delta_y_t / y_t_core_mt = sigma_lambda = -(9/16) * delta_lambda / lambda_core_mt",
            "sigma_center": sigma_shared,
            "sigma_half_width": abs(eta_ht),
            "sigma_interval": [
                min(sigma_from_y, sigma_from_lambda),
                max(sigma_from_y, sigma_from_lambda),
            ],
            "seed_equality_residual_abs": abs(w_ht),
            "rank_one_residual_abs": abs(w_ht),
            "common_provenance": {
                "status": "constructive_for_collapsed_seed",
                "seed_value": sigma_shared,
                "strict_open_object_name": "D11AntidiagonalReadbackVanishing",
                "strict_open_object_value": eta_ht,
            },
            "normalized_readback": {
                "coordinates": {
                    "pi_y": sigma_from_y,
                    "pi_lambda": sigma_from_lambda,
                },
                "projectors": {
                    "common_seed_projector": "(pi_y + pi_lambda) / 2",
                    "diagnostic_witness_projector": "(pi_y - pi_lambda) / 2",
                },
                "decomposition": {
                    "sigma_shared": sigma_shared,
                    "eta_antidiagonal": eta_ht,
                    "witness_zero_symbol": "w_HT",
                },
            },
            "common_seed_emitter": {
                "seed_symbol": "sigma_shared",
                "emitted_theta": {
                    "delta_y_t_mt": shared_delta_y,
                    "delta_lambda_mt": shared_delta_lambda,
                },
            },
            "common_provenance_certificate": {
                "status": "constructive_for_collapsed_seed",
                "zero_witness_symbol": "w_HT",
                "witness_kind": "normalized_seed_difference",
                "witness_formula": "sigma_from_y_t - sigma_from_lambda",
                "witness_abs": abs(w_ht),
                "strict_open_object_name": "D11AntidiagonalReadbackVanishing",
                "strict_open_object_value": eta_ht,
                "closes_when": 0.0,
                "collapsed_seed_formula": "0.5 * (sigma_from_y_t + sigma_from_lambda)",
                "collapsed_seed_value": sigma_shared,
            },
            "antidiagonal_vanishing_certificate": {
                "name": "D11AntidiagonalReadbackVanishing",
                "status": "candidate",
                "proof_mode": "single_scalar_factorization",
                "eta_formula": "0.5 * (pi_y - pi_lambda)",
                "w_formula": "pi_y - pi_lambda",
                "current_eta": eta_ht,
                "current_w": w_ht,
                "closes_to": {
                    "eta_HT": 0.0,
                    "w_HT": 0.0,
                },
            },
            "single_scalar_factorization_certificate": {
                "status": "exact_center_branch_closed",
                "pi_y_equals_pi_lambda": sigma_shared,
                "eta_HT": 0.0,
                "w_HT": 0.0,
            },
            "diagnostic_wedge_scalar": {
                "name": "omega_HT",
                "formula": "kappa_HT * lambda_core_mt * delta_y_t_mt + y_t_core_mt * delta_lambda_mt",
                "value": omega_ht,
                "eta_formula_from_wedge": "omega_HT / (2 * kappa_HT * lambda_core_mt * y_t_core_mt)",
                "strict_open_object_name": "D11FixedRayWedgeVanishing",
            },
            "exact_center_promotion": {
                "status": "diagnostic_only",
                "mt_pole_gev": 172.40055484582462,
                "mH_gev": 125.21106404382539,
                "eta_HT": 0.0,
                "w_HT": 0.0,
            },
            "remaining_exact_object": "D11FixedRayWedgeVanishing",
            "no_independent_higgs_top_residuals": True,
            "no_mu_sync_rescan": True,
            "delta_y_t_mt": shared_delta_y,
            "delta_lambda_mt": shared_delta_lambda,
        },
        "diagnostic_required_to_current_refs": {
            "ref_surface": str(results_status),
            "ref_mt_pole_gev": 172.4,
            "ref_mH_gev": 125.20,
            "delta_y_t_mt": delta_y,
            "delta_lambda_mt": delta_lambda,
            "delta_mt_pole_gev": 2.13875,
            "delta_mH_gev": -1.42263,
        },
        "legacy_diagnostic_readback": {
            "delta_y_t_mt": delta_y,
            "delta_lambda_mt": delta_lambda,
            "source": "diagnostic_required_to_current_refs",
        },
        "fixed_ray_slope_certificate": {
            "name": "D11FixedRaySlopeCertificate",
            "status": "candidate",
            "kappa_candidate": kappa_candidate,
            "kappa_observed_formula": "- y_t_core_mt * legacy_diagnostic_readback.delta_lambda_mt / (lambda_core_mt * legacy_diagnostic_readback.delta_y_t_mt)",
            "kappa_observed_value": kappa_observed,
            "kappa_residual": kappa_residual,
            "prefactor_formula": "lambda_core_mt * legacy_diagnostic_readback.delta_y_t_mt",
            "prefactor_value": prefactor,
            "wedge_factorization_formula": "omega_HT = prefactor_value * (kappa_candidate - kappa_observed_value)",
            "wedge_value": omega_ht,
            "closes_when": {
                "kappa_residual": 0.0,
                "omega_HT": 0.0,
                "eta_HT": 0.0,
                "w_HT": 0.0,
            },
        },
        "jacobian": {
            "d_mt_pole_d_y_t": 184.97,
            "d_mH_d_lambda": 480.0,
        },
        "predicted": {
            "mt_pole_gev": 172.40055484582462,
            "mH_gev": 125.21106404382539,
        },
        "center_witness_forward_consequence": {
            "sigma_D11_HT": sigma_shared,
            "mt_pole_gev": 172.40055484582462,
            "mH_gev": 125.21106404382539,
            "mt_pole_interval_gev": [172.3999600759, 172.40114961574923],
            "mH_interval_gev": [125.2106716, 125.21145648765076],
        },
        "exact_missing_object": "D11FixedRayWedgeVanishing",
        "exact_missing_object_scope": "legacy_diagnostic_sidecar_only" if forward_closed else "live_forward_path",
        "upstream_missing_forward_artifact": None if forward_closed else "sigma_D11_HT_or_full_Theta_D11_HT",
        "forward_path_closed_by": "forward_seed_promotion_certificate" if forward_closed else None,
        "live_forward_path_missing_object": None if forward_closed else "sigma_D11_HT_or_full_Theta_D11_HT",
        "readout_vector_symbol": "Theta_D11_HT(mu_t) = (delta_y_t(mu_t), delta_lambda(mu_t))",
        "closure_state": (
            "diagnostic_sidecar_only__live_forward_path_closed_elsewhere"
            if forward_closed
            else "diagnostic_center_projection_only__not_forward_closed"
        ),
        "notes": [
            "The synchronized D11 core already captures most of the numerical gain over the literal appendix flow, but mu_sync alone moves top and Higgs together and cannot do the needed top-up / Higgs-down move.",
            "The smallest constructive D11 object is a common low-scale readout vector Theta_D11_HT(mu_t) = (delta_y_t, delta_lambda) produced by one CriticalSurfaceReadoutKernel_D11, not two independent per-observable residual fits.",
            "The current diagnostic vector already lies almost exactly on a rank-one relative-core ray with kappa_HT = 16/9, so the collapsed shared seed sigma_D11_HT is now constructively available rather than only interval-marked.",
            "The normalized readbacks split exactly into a common seed sigma_shared and one strict antidiagonal coordinate eta_HT, with w_HT = 2 * eta_HT.",
            "The exact emitted center branch is now closed by single-scalar factorization, with pi_y = pi_lambda = sigma_shared and eta_HT = w_HT = 0 on that branch.",
            "The smallest remaining exact diagnostic object is therefore not the full common-provenance certificate again but the fixed-ray wedge scalar D11FixedRayWedgeVanishing.",
            "At the current boundary, the collapsed shared seed already gives a concrete local forward consequence for mt_pole and mH without reopening any separate Higgs/top residual family.",
            (
                "The remaining open object lives only on the legacy diagnostic skew; the live forward branch is closed elsewhere by the forward-seed promotion certificate."
                if forward_closed
                else "The remaining open object lives only on the legacy diagnostic skew, not on the emitted exact-center branch."
            ),
            (
                "The current center branch remains diagnostic-only inside this sidecar because the live forward closure is recorded on the separate forward-seed certificate."
                if forward_closed
                else "The current center branch remains diagnostic-only because sigma_shared is still seeded from diagnostic delta_y_t and delta_lambda readback rather than a forward-emitted D11 artifact."
            ),
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D11 critical-surface readout boundary artifact.")
    parser.add_argument("--d10-source", default=str(DEFAULT_D10_SOURCE))
    parser.add_argument("--results-status", default=str(DEFAULT_RESULTS_STATUS))
    parser.add_argument("--forward-seed-certificate", default=str(DEFAULT_FORWARD_SEED_CERTIFICATE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    certificate_path = Path(args.forward_seed_certificate)
    certificate = (
        json.loads(certificate_path.read_text(encoding="utf-8"))
        if certificate_path.exists()
        else None
    )
    artifact = build_artifact(Path(args.d10_source), Path(args.results_status), certificate)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
