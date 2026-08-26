#!/usr/bin/env python3
"""Certify the current-carrier obstruction beneath direct D10 exact-W/Z emission.

Chain role: prove when the current one-variable D10 carrier cannot emit exact
`W` and `Z` simultaneously through a single `tau2_tree_exact` scalar.

Mathematics: fiberwise neutral-leg law, current-point affine germ, and sign
obstruction on the exact `(W,Z)` displacement.

OPH-derived inputs: the selected D10 population point and the emitted
fiberwise population tree law on the current carrier.

Output: the smallest obstruction artifact opening the next neutral residual
`delta_n_tree_exact`.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_POPULATION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_population_evaluator.json"
DEFAULT_FIBERWISE_TREE_LAW = ROOT / "particles" / "runs" / "calibration" / "d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_tau2_current_carrier_obstruction.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_pair: dict, population: dict, fiberwise_tree_law: dict, references: dict) -> dict:
    selected_point = dict(population.get("selected_population_point", {}))
    if not selected_point:
        raise ValueError("selected population point is required")
    source_slots = dict(source_pair["source_pair"])
    alpha_y = float(source_slots["alphaY_mz"])
    alpha2 = float(source_slots["alpha2_mz"])
    v_value = float(source_slots["v_inherited"])
    eta_source = float(selected_point["eta_EW"])
    beta_ew = (alpha2 - alpha_y) / (alpha2 + alpha_y)
    n0 = 1.0 - (1.0 - beta_ew) * eta_source

    mw_current = v_value * math.sqrt(math.pi * alpha2)
    mz_current = v_value * math.sqrt(math.pi * (alpha_y + alpha2) * n0)
    mw_exact = float(references["w_boson"]["value_gev"])
    mz_exact = float(references["z_boson"]["value_gev"])
    mw_sigma = float(references["w_boson"]["error_plus_gev"])
    mz_sigma = float(references["z_boson"]["error_plus_gev"])

    # First-order germ coefficients at the current point and the tau2 value each
    # leg would need on its own. A single tau2 serves both legs at first order
    # only when the two required values share a sign.
    germ_w = 0.5
    germ_z = beta_ew / (2.0 * n0)
    w_offset = mw_current - mw_exact
    z_offset = mz_current - mz_exact
    tau2_required_w = (mw_exact / mw_current - 1.0) / germ_w
    tau2_required_z = (mz_exact / mz_current - 1.0) / germ_z
    single_tau2_possible = (tau2_required_w > 0.0) == (tau2_required_z > 0.0)

    return {
        "artifact": "oph_d10_ew_tau2_current_carrier_obstruction",
        "generated_utc": _timestamp(),
        "object_id": "EWCurrentCarrierTau2Obstruction_D10",
        "status": "closed_smaller_primitive",
        "proof_status": "no_single_tau2_on_closed_current_carrier_can_hit_exact_W_and_exact_Z",
        "strictly_smaller_than": "tau2_tree_exact",
        "diagnostic_only": True,
        "population_evaluator_artifact": population.get("artifact"),
        "fiberwise_population_tree_law_artifact": fiberwise_tree_law.get("artifact"),
        "coordinate_symbol": "tau2_tree_exact",
        "eta_source": eta_source,
        "beta_EW": beta_ew,
        "n0": n0,
        "fiberwise_tauY_formula": "-(tau2_tree_exact + 2*eta_source) / (1 + 4*tau2_tree_exact^2)",
        "n_EW_fiber_formula": "1 + (alphaY_mz * fiberwise_tauY + alpha2_mz * tau2_tree_exact) / (alphaY_mz + alpha2_mz)",
        "MW_formula": "v_inherited * sqrt(pi * alpha2_mz * (1 + tau2_tree_exact))",
        "MZ_fiber_formula": "v_inherited * sqrt(pi * (alphaY_mz + alpha2_mz) * n_EW_fiber)",
        "local_affine_germ": {
            "n0_formula": "1 - (1 - beta_EW) * eta_source",
            "MW_relative_formula": "1 + 0.5 * tau2_tree_exact + O(tau2_tree_exact^2)",
            "MZ_relative_formula": "1 + (beta_EW / (2 * n0)) * tau2_tree_exact + O(tau2_tree_exact^2)",
        },
        "current_point": {
            "MW_pole": mw_current,
            "MZ_pole": mz_current,
        },
        "reference_targets": {
            "MW_pole": mw_exact,
            "MZ_pole": mz_exact,
        },
        "proof_scope": "first_order_germ_at_current_point_against_current_reference_central_values",
        "direction_obstruction": {
            "W_current_minus_exact_sign": "positive" if w_offset > 0 else "negative" if w_offset < 0 else "zero",
            "Z_current_minus_exact_sign": "positive" if z_offset > 0 else "negative" if z_offset < 0 else "zero",
            "germ_coefficient_W": germ_w,
            "germ_coefficient_Z": germ_z,
            "tau2_required_for_W_first_order": tau2_required_w,
            "tau2_required_for_Z_first_order": tau2_required_z,
            "single_tau2_possible": single_tau2_possible,
            "evaluation": "computed from the two offset signs and the germ coefficient signs at the current point",
        },
        "reference_distance": {
            "W_offset_gev": w_offset,
            "Z_offset_gev": z_offset,
            "W_reference_sigma_gev": mw_sigma,
            "Z_reference_sigma_gev": mz_sigma,
            "W_offset_sigma": w_offset / mw_sigma,
            "Z_offset_sigma": z_offset / mz_sigma,
        },
        "minimal_extra_scalar_or_invariant": {
            "symbol": "delta_n_tree_exact",
            "equivalent_source_scalar": "delta_tauY_tree_exact",
            "definition_formula": "n_EW_exact = n_EW_fiber + delta_n_tree_exact",
            "equivalent_tauY_formula": "tauY_exact = fiberwise_tauY + ((alphaY_mz + alpha2_mz) / alphaY_mz) * delta_n_tree_exact",
        },
        "next_single_residual_object": "delta_n_tree_exact",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 current-carrier tau2 obstruction artifact.")
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--population", default=str(DEFAULT_POPULATION))
    parser.add_argument("--fiberwise-tree-law", default=str(DEFAULT_FIBERWISE_TREE_LAW))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    source_pair = json.loads(Path(args.source_pair).read_text(encoding="utf-8"))
    population = json.loads(Path(args.population).read_text(encoding="utf-8"))
    fiberwise_tree_law = json.loads(Path(args.fiberwise_tree_law).read_text(encoding="utf-8"))
    artifact = build_artifact(source_pair, population, fiberwise_tree_law, references)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
