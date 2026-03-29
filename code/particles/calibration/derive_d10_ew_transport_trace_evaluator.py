#!/usr/bin/env python3
"""Export the D10 fixed-eta transport-trace evaluator artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_READOUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_readout.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_transport_trace_evaluator.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_pair: dict, readout: dict) -> dict:
    fixed_eta = dict(readout["fixed_eta_slice"])
    return {
        "artifact": "oph_d10_ew_transport_trace_evaluator",
        "generated_utc": _timestamp(),
        "object_id": "EWTransportTraceEvaluator_D10",
        "status": "open",
        "family_source_id": readout.get("family_source_id"),
        "scheme_id": readout.get("scheme_id"),
        "origin_kernel_id": readout.get("origin_kernel_id"),
        "fixed_eta_slice": {
            "alphaY_0": float(fixed_eta["alphaY_0"]),
            "alpha2_0": float(fixed_eta["alpha2_0"]),
            "v_report": float(readout["base_running_quintet"]["v_report"]),
            "sin2_thetaW0": float(fixed_eta["sin2_thetaW0"]),
            "cos2_thetaW0": float(fixed_eta["cos2_thetaW0"]),
            "eta_EW": float(fixed_eta["eta_EW"]),
        },
        "sigma_mode": {
            "symbol": "sigma_EW",
            "role": "common_trace_mode",
            "formula_from_source_pair": "0.5 * (tau_Y + tau_2)",
            "formula_from_transport_trace": "0.5 * (Pi_AA + Pi_ZZ)",
            "formula_from_W_entry": "Pi_WW - eta_EW",
            "physical_domain": "sigma_EW > eta_EW - 1",
        },
        "transport_entries_from_sigma": {
            "tau_Y": "sigma_EW - eta_EW",
            "tau_2": "sigma_EW + eta_EW",
            "Pi_AA": "sigma_EW - cos2_thetaW0 * eta_EW",
            "Pi_AZ": "2 * sin(theta_W0) * cos(theta_W0) * eta_EW",
            "Pi_ZZ": "sigma_EW + cos2_thetaW0 * eta_EW",
            "Pi_WW": "sigma_EW + eta_EW",
        },
        "shared_scalar_image_from_sigma": {
            "delta_rho": "sqrt(1 + sigma_EW + cos2_thetaW0 * eta_EW) - 1",
            "delta_rW": "sqrt(1 + sigma_EW + eta_EW) - sqrt(1 + sigma_EW + cos2_thetaW0 * eta_EW)",
            "delta_kappa": "((1 + sigma_EW - eta_EW) / (1 + sigma_EW + cos2_thetaW0 * eta_EW)) - 1",
            "delta_alpha": "((1/(alphaY_0*(1 + sigma_EW - eta_EW)) + 1/(alpha2_0*(1 + sigma_EW + eta_EW))) / (1/alphaY_0 + 1/alpha2_0)) - 1",
        },
        "quintet_from_sigma": {
            "alphaY_star": "alphaY_0 * (1 + sigma_EW - eta_EW)",
            "alpha2_star": "alpha2_0 * (1 + sigma_EW + eta_EW)",
            "mW": "v_report * sqrt(pi * alpha2_star)",
            "mZ": "v_report * sqrt(pi * (alphaY_star + alpha2_star))",
            "alpha_em_inv": "(alphaY_star + alpha2_star) / (alphaY_star * alpha2_star)",
            "sin2_thetaW": "alphaY_star / (alphaY_star + alpha2_star)",
        },
        "insufficiency_certificate": {
            "tree_identity_holds_for_all_sigma": True,
            "sigma_not_selected_by_current_family": True,
            "running_weak_angle_outside_slice_image": True,
            "running_alpha_conflicts_with_near_coherent_mass_pair": True,
            "forbidden_live_inverse_witnesses": [
                "sigma_from_mW",
                "sigma_from_mZ",
                "sigma_from_alpha_em_inv",
                "sigma_from_sin2_thetaW",
            ],
        },
        "current_family_anchor": {
            "source_pair_artifact": source_pair.get("artifact"),
            "readout_artifact": readout.get("artifact"),
            "compact_point": readout.get("current_compact_point"),
            "compact_emitted_quintet": readout.get("current_compact_emitted_quintet"),
        },
        "diagnostic_only": True,
        "smaller_live_object": "EWSinglePostTransportTreeIdentity_D10",
        "smallest_constructive_missing_object": "EWSinglePostTransportTreeIdentity_D10",
        "notes": [
            "The fixed-eta one-sigma family is coherent but underdetermined as a predictive D10 closure branch.",
            "The running weak-angle target lies outside the fixed-eta one-sigma slice image, with alpha_em^-1 also drifting away from the near-coherent W/Z sigma.",
            "This artifact is diagnostic only: the carrier-level selector and split exact neutral closure are already closed, and the live D10 gap is now the stronger unsplit post-transport tree identity.",
            "No inverse witness from W, Z, alpha_em^-1, or sin^2(theta_W) is permitted here.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 transport-trace evaluator artifact.")
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--readout", default=str(DEFAULT_READOUT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_pair = json.loads(Path(args.source_pair).read_text(encoding="utf-8"))
    readout = json.loads(Path(args.readout).read_text(encoding="utf-8"))
    artifact = build_artifact(source_pair, readout)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
