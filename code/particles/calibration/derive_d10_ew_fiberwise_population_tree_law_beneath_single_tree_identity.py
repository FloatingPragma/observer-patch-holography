#!/usr/bin/env python3
"""Emit the fixed-fiber D10 tree law beneath the unsplit exact W/Z shell.

Chain role: reduce the open unsplit D10 tree identity to a one-variable law on
the existing selected carrier, so only the charged-leg scalar `tau2_tree_exact`
remains open.

Mathematics: minimize the closed D10 population functional on each fixed
`tau_2` fiber to derive the unique `tau_Y = Phi_EW_fiber(tau_2)` law.

OPH-derived inputs: the selected D10 carrier point, the closed D10 population
functional, and the live exact-W/Z coordinate shell.

Output: the smaller fiberwise tree-law primitive that removes the placeholder
`Phi_EW_tree(tau_2)` without inventing `tau2_tree_exact`.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_POPULATION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_population_evaluator.json"
DEFAULT_EXACT_WZ = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_wz_coordinate_beyond_single_tree_identity.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_pair: dict, population: dict, exact_wz_coordinate: dict | None = None) -> dict:
    source_slots = dict(source_pair["source_pair"])
    compact_mass_slice = dict(source_pair["compact_hypercharge_only_mass_slice"])
    selected_point = dict(population.get("selected_population_point", {}))
    if not selected_point:
        raise ValueError("selected population point is required")

    alpha_y0 = float(source_slots["alphaY_mz"])
    alpha2_0 = float(source_slots["alpha2_mz"])
    v_value = float(source_slots["v_inherited"])
    eta_source = float(compact_mass_slice["eta_EW"])
    beta_ew = (alpha2_0 - alpha_y0) / (alpha2_0 + alpha_y0)
    tau2_anchor = float(selected_point["tau_2"])
    tauy_anchor = -((tau2_anchor + 2.0 * eta_source) / (1.0 + 4.0 * tau2_anchor * tau2_anchor))
    sigma_anchor = ((2.0 * tau2_anchor**3) - eta_source) / (1.0 + 4.0 * tau2_anchor * tau2_anchor)
    eta_anchor = (eta_source + tau2_anchor + 2.0 * tau2_anchor**3) / (1.0 + 4.0 * tau2_anchor * tau2_anchor)

    return {
        "artifact": "oph_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity",
        "generated_utc": _timestamp(),
        "object_id": "EWFiberwisePopulationTreeLaw_D10",
        "status": "closed_smaller_primitive",
        "proof_status": "fiberwise_unique_J_pop_minimizer_on_fixed_tau2",
        "strictly_smaller_than": "EWSinglePostTransportTreeIdentity_D10",
        "derived_from_artifacts": [
            population.get("artifact"),
            None if exact_wz_coordinate is None else exact_wz_coordinate.get("artifact"),
        ],
        "coordinate_symbol": "tau2_tree_exact",
        "tree_law_symbol": "Phi_EW_fiber",
        "carrier_basis_scalar": {
            "alphaY_mz": alpha_y0,
            "alpha2_mz": alpha2_0,
            "beta_EW": beta_ew,
            "beta_EW_formula": "(alpha2_mz - alphaY_mz) / (alpha2_mz + alphaY_mz)",
            "v_inherited": v_value,
        },
        "eta_source": eta_source,
        "eta_source_formula": "alpha_u_from_seed * beta_EW",
        "fiber_population_functional_formula": "J_pop_EW(tauY,tau2) = tau2^2 + (tauY*tau2)^2 + (0.5*(tauY + tau2) + eta_source)^2",
        "fiber_stationarity_formula": "(1 + 4*tau2_tree_exact^2)*tau_Y + tau2_tree_exact + 2*eta_source = 0",
        "fiber_second_derivative_formula": "1/2 + 2*tau2_tree_exact^2",
        "tauY_formula": "-(tau2_tree_exact + 2*eta_source) / (1 + 4*tau2_tree_exact^2)",
        "sigma_formula": "(2*tau2_tree_exact^3 - eta_source) / (1 + 4*tau2_tree_exact^2)",
        "eta_formula": "(eta_source + tau2_tree_exact + 2*tau2_tree_exact^3) / (1 + 4*tau2_tree_exact^2)",
        "u_EW_formula": "1 + tau2_tree_exact",
        "n_EW_formula": "1 + (alphaY_mz * tau_Y + alpha2_mz * tau2_tree_exact) / (alphaY_mz + alpha2_mz)",
        "sigma_from_coordinate_formula": "0.5 * (tau_Y + tau2_tree_exact)",
        "eta_from_coordinate_formula": "0.5 * (tau2_tree_exact - tau_Y)",
        "MW_formula": "v_inherited * sqrt(pi * alpha2_mz * (1 + tau2_tree_exact))",
        "MZ_formula": "v_inherited * sqrt(pi * (alphaY_mz * (1 + tau_Y) + alpha2_mz * (1 + tau2_tree_exact)))",
        "anchor_point": {
            "tau2_tree_exact": tau2_anchor,
            "tau_Y": tauy_anchor,
            "sigma_EW": sigma_anchor,
            "eta_EW": eta_anchor,
        },
        "next_single_residual_object": "tau2_tree_exact",
        "notes": [
            "The fixed-tau2 fibers of J_pop_EW have a unique tauY minimizer.",
            "This removes the placeholder unsplit tree law without inventing tau2_tree_exact.",
            "Exact W and exact Z become a one-variable forward law on the existing carrier once this primitive is emitted.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 fiberwise population tree law beneath the unsplit exact-W/Z shell.")
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--population", default=str(DEFAULT_POPULATION))
    parser.add_argument("--exact-wz-coordinate", default=str(DEFAULT_EXACT_WZ))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_pair = json.loads(Path(args.source_pair).read_text(encoding="utf-8"))
    population = json.loads(Path(args.population).read_text(encoding="utf-8"))
    exact_wz_path = Path(args.exact_wz_coordinate)
    exact_wz_coordinate = json.loads(exact_wz_path.read_text(encoding="utf-8")) if exact_wz_path.exists() else None
    artifact = build_artifact(source_pair, population, exact_wz_coordinate)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
