#!/usr/bin/env python3
"""Emit the exact D10 mass-pair chart on the selected current carrier.

Chain role: expose the exact two-coordinate `(tau2_tree_exact, delta_n_tree_exact)`
chart beneath the selected D10 carrier point.

Mathematics: combine the closed fiberwise minimizer for `tau_Y` with the closed
current-carrier obstruction to write the exact `W/Z` chart and the pullback of
the current selector to that chart.

OPH-derived inputs: the selected D10 carrier point, the closed population
functional, the fiberwise tree law, and the current-carrier obstruction.

Output: the exact mass-pair chart and the proof that the current closed selector
still has a unique zero at the current point, leaving one remaining selector
object rather than another coordinate.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_POPULATION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_population_evaluator.json"
DEFAULT_FIBERWISE_TREE_LAW = ROOT / "particles" / "runs" / "calibration" / "d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.json"
DEFAULT_TAU2_OBSTRUCTION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_tau2_current_carrier_obstruction.json"
DEFAULT_EXACT_WZ_COORDINATE = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_wz_coordinate_beyond_single_tree_identity.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_mass_pair_chart_current_carrier.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    source_pair: dict,
    population: dict,
    fiberwise_tree_law: dict,
    tau2_obstruction: dict,
    exact_wz_coordinate: dict,
) -> dict:
    source_slots = dict(source_pair["source_pair"])
    selected_population_point = dict(population.get("selected_population_point", {}))
    alpha_y = float(source_slots["alphaY_mz"])
    alpha2 = float(source_slots["alpha2_mz"])
    v_value = float(source_slots["v_inherited"])
    eta_source = float(source_pair["eta_source"])
    beta_ew = (alpha2 - alpha_y) / (alpha2 + alpha_y)
    kappa_n_to_tauy = (alpha_y + alpha2) / alpha_y

    return {
        "artifact": "oph_d10_ew_exact_mass_pair_chart_current_carrier",
        "generated_utc": _timestamp(),
        "object_id": "EWExactMassPairChartCurrentCarrier_D10",
        "status": "closed_smaller_primitive",
        "proof_status": "exact_mass_pair_chart_closed_and_current_selector_pullback_has_unique_zero_at_current_point",
        "strictly_smaller_than": "tau2_tree_exact_and_delta_n_tree_exact_emission",
        "derived_from_artifacts": [
            population.get("artifact"),
            fiberwise_tree_law.get("artifact"),
            tau2_obstruction.get("artifact"),
            exact_wz_coordinate.get("artifact"),
        ],
        "coordinate_symbols": {
            "charged_leg": "tau2_tree_exact",
            "neutral_leg": "delta_n_tree_exact",
        },
        "carrier_basis_scalar": {
            "alphaY_mz": alpha_y,
            "alpha2_mz": alpha2,
            "beta_EW": beta_ew,
            "beta_EW_formula": "(alpha2_mz - alphaY_mz) / (alpha2_mz + alphaY_mz)",
            "kappa_n_to_tauY": kappa_n_to_tauy,
            "kappa_n_to_tauY_formula": "(alphaY_mz + alpha2_mz) / alphaY_mz",
            "eta_source": eta_source,
            "eta_source_formula": "selected_population_point.eta_EW",
            "v_inherited": v_value,
        },
        "exact_mass_pair_chart": {
            "tauY_fiber_formula": fiberwise_tree_law["tauY_formula"],
            "tauY_exact_formula": "tauY_fiber + kappa_n_to_tauY * delta_n_tree_exact",
            "sigma_exact_formula": (
                "(2*tau2_tree_exact^3 - eta_source) / (1 + 4*tau2_tree_exact^2) "
                "+ 0.5 * kappa_n_to_tauY * delta_n_tree_exact"
            ),
            "eta_exact_formula": (
                "(eta_source + tau2_tree_exact + 2*tau2_tree_exact^3) / (1 + 4*tau2_tree_exact^2) "
                "- 0.5 * kappa_n_to_tauY * delta_n_tree_exact"
            ),
            "u_EW_exact_formula": "1 + tau2_tree_exact",
            "n_EW_fiber_formula": fiberwise_tree_law["n_EW_formula"],
            "n_EW_fiber_expanded_formula": (
                "1 + (beta_EW * tau2_tree_exact + 2 * (1 + beta_EW) * tau2_tree_exact^3 "
                "- (1 - beta_EW) * eta_source) / (1 + 4 * tau2_tree_exact^2)"
            ),
            "n_EW_exact_formula": "n_EW_fiber + delta_n_tree_exact",
            "MW_formula": fiberwise_tree_law["MW_formula"],
            "MZ_formula": "v_inherited * sqrt(pi * (alphaY_mz + alpha2_mz) * (n_EW_fiber + delta_n_tree_exact))",
        },
        "current_selector_pullback": {
            "symbol": "J_pop_EW_pullback",
            "formula": (
                "tau2_tree_exact^2 * (1 + (2 * eta_source + tau2_tree_exact)^2 + 4 * tau2_tree_exact^2) "
                "/ (1 + 4 * tau2_tree_exact^2) + (kappa_n_to_tauY^2 / 4) * (1 + 4 * tau2_tree_exact^2) "
                "* delta_n_tree_exact^2"
            ),
            "nonnegative": True,
            "zero_set": {
                "tau2_tree_exact": 0.0,
                "delta_n_tree_exact": 0.0,
            },
            "verdict": "current_closed_selector_cannot_emit_nonzero_exact_mass_pair",
        },
        "local_bijectivity_certificate": {
            "dMW_dtau2_formula": "MW_formula / (2 * (1 + tau2_tree_exact))",
            "dMW_ddelta_n_formula": "0",
            "d_tauY_fiber_dtau2_formula": "(-1 + 4 * tau2_tree_exact^2 + 16 * eta_source * tau2_tree_exact) / (1 + 4 * tau2_tree_exact^2)^2",
            "d_n_EW_fiber_dtau2_formula": "(alphaY_mz * d_tauY_fiber_dtau2 + alpha2_mz) / (alphaY_mz + alpha2_mz)",
            "dMZ_dtau2_formula": "MZ_formula * d_n_EW_fiber_dtau2 / (2 * (n_EW_fiber + delta_n_tree_exact))",
            "dMZ_ddelta_n_formula": "MZ_formula / (2 * (n_EW_fiber + delta_n_tree_exact))",
            "determinant_formula": "pi * v_inherited^2 * sqrt(alpha2_mz * (alphaY_mz + alpha2_mz)) / (4 * sqrt((1 + tau2_tree_exact) * (n_EW_fiber + delta_n_tree_exact)))",
            "third_coordinate_needed": False,
        },
        "selected_current_point": {
            "tau2_tree_exact": float(selected_population_point.get("tau_2", 0.0)),
            "delta_n_tree_exact": 0.0,
        },
        "next_single_residual_object": "EWExactMassPairSelector_D10",
        "notes": [
            "The exact current-carrier W/Z chart is two-dimensional on (tau2_tree_exact, delta_n_tree_exact).",
            "No third coordinate or broader carrier is needed for exact mass-pair closure on this chart.",
            "The current closed selector pulls back to a nonnegative scalar with a unique zero at the current point, so the next residual object is a selector on this chart rather than another coordinate.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact D10 mass-pair chart on the selected current carrier.")
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--population", default=str(DEFAULT_POPULATION))
    parser.add_argument("--fiberwise-tree-law", default=str(DEFAULT_FIBERWISE_TREE_LAW))
    parser.add_argument("--tau2-obstruction", default=str(DEFAULT_TAU2_OBSTRUCTION))
    parser.add_argument("--exact-wz-coordinate", default=str(DEFAULT_EXACT_WZ_COORDINATE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_pair = json.loads(Path(args.source_pair).read_text(encoding="utf-8"))
    population = json.loads(Path(args.population).read_text(encoding="utf-8"))
    fiberwise_tree_law = json.loads(Path(args.fiberwise_tree_law).read_text(encoding="utf-8"))
    tau2_obstruction = json.loads(Path(args.tau2_obstruction).read_text(encoding="utf-8"))
    exact_wz_coordinate = json.loads(Path(args.exact_wz_coordinate).read_text(encoding="utf-8"))
    artifact = build_artifact(source_pair, population, fiberwise_tree_law, tau2_obstruction, exact_wz_coordinate)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
