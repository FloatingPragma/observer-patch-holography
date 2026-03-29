#!/usr/bin/env python3
"""Emit the fixed-eta affine germ beneath the open D10 unsplit tree identity.

Chain role: expose the strongest exact smaller primitive beneath the open
single-tree identity on the selected D10 carrier.

Mathematics: fixed-`eta_EW` affine probe family around the selected carrier
anchor, with `tau_Y` and `sigma_EW` determined linearly by one probe scalar.

OPH-derived inputs: the selected D10 source transport pair and the closed
population selector point on that carrier.

Output: the diagnostic affine germ that reproduces the selected public `W/Z`
pair while keeping full electroweak exact closure open.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_POPULATION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_population_evaluator.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_fixed_eta_post_transport_affine_germ.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_pair: dict, population: dict) -> dict:
    source_slots = dict(source_pair["source_pair"])
    selected_point = dict(population.get("selected_population_point", {}))
    if not selected_point:
        raise ValueError("selected population point is required")

    alpha_y0 = float(source_slots["alphaY_mz"])
    alpha2_0 = float(source_slots["alpha2_mz"])
    v_value = float(source_slots["v_inherited"])
    beta_ew = (alpha2_0 - alpha_y0) / (alpha2_0 + alpha_y0)
    eta_tree_anchor = float(selected_point["eta_EW"])
    sigma_anchor = float(selected_point["sigma_EW"])
    tau_y_anchor = float(selected_point["tau_Y"])
    tau_2_anchor = float(selected_point["tau_2"])

    return {
        "artifact": "oph_d10_ew_fixed_eta_post_transport_affine_germ",
        "generated_utc": _timestamp(),
        "object_id": "EWFixedEtaPostTransportAffineGerm_D10",
        "status": "closed_smaller_primitive",
        "proof_status": "exact_fixed_eta_affine_germ_extracted_from_selected_carrier",
        "strictly_smaller_than": "EWSinglePostTransportTreeIdentity_D10",
        "diagnostic_only": True,
        "why_not_promoted": "fixed_eta_slice_remains_diagnostic_and_nonexact_for_full_electroweak_target_surface",
        "source_transport_pair_artifact": source_pair.get("artifact"),
        "population_evaluator_artifact": population.get("artifact"),
        "alphaY_mz": alpha_y0,
        "alpha2_mz": alpha2_0,
        "v_inherited": v_value,
        "beta_EW": beta_ew,
        "beta_EW_formula": "(alpha2_mz - alphaY_mz) / (alpha2_mz + alphaY_mz)",
        "eta_tree_anchor": eta_tree_anchor,
        "eta_tree_anchor_formula": "alpha_u_from_seed * beta_EW",
        "selected_anchor_point": {
            "sigma_EW": sigma_anchor,
            "eta_EW": eta_tree_anchor,
            "tau_Y": tau_y_anchor,
            "tau_2": tau_2_anchor,
        },
        "forward_probe_family": {
            "parameter_symbol": "tau2_probe",
            "eta_EW_formula": "eta_tree_anchor",
            "sigma_EW_formula": "tau2_probe - eta_tree_anchor",
            "tau_2_formula": "tau2_probe",
            "tau_Y_formula": "tau2_probe - 2 * eta_tree_anchor",
            "alphaY_star_formula": "alphaY_mz * (1 + tau2_probe - 2 * eta_tree_anchor)",
            "alpha2_star_formula": "alpha2_mz * (1 + tau2_probe)",
            "u_EW_formula": "1 + tau2_probe",
            "n_EW_formula": "1 + (alphaY_mz * (tau2_probe - 2 * eta_tree_anchor) + alpha2_mz * tau2_probe) / (alphaY_mz + alpha2_mz)",
            "MW_formula": "v_inherited * sqrt(pi * alpha2_mz * (1 + tau2_probe))",
            "MZ_formula": "v_inherited * sqrt(pi * (alphaY_mz * (1 + tau2_probe - 2 * eta_tree_anchor) + alpha2_mz * (1 + tau2_probe)))",
        },
        "next_residual_object": "EWSinglePostTransportTreeIdentity_D10",
        "notes": [
            "This affine germ is the strongest exact smaller primitive visible on the selected carrier without promoting the open unsplit tree identity.",
            "At tau2_probe = 0 it reproduces the selected public W/Z pair exactly.",
            "Full electroweak exactness remains open because the fixed-eta affine germ is only diagnostic, not the final tree identity.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 fixed-eta post-transport affine germ.")
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--population", default=str(DEFAULT_POPULATION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_pair = json.loads(Path(args.source_pair).read_text(encoding="utf-8"))
    population = json.loads(Path(args.population).read_text(encoding="utf-8"))
    artifact = build_artifact(source_pair, population)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
