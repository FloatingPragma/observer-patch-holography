#!/usr/bin/env python3
"""Expose the first mass-moving D10 coordinate beneath exact W/Z closure.

Chain role: isolate the smallest charged-leg post-transport coordinate that can
move the public `W/Z` pair once the unsplit D10 tree identity is available.

Mathematics: one-scalar charged-leg transport coordinate above the selected
carrier point, with the neutral leg determined by the still-open unsplit tree
identity.

OPH-derived inputs: the selected D10 carrier point, the reopened two-scalar
transport pair, and the carrier selector that closes the current public point.

Output: the exact-W/Z coordinate shell together with the current-carrier
obstruction that opens the next neutral residual `delta_n_tree_exact`.
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
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_wz_coordinate_beyond_single_tree_identity.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    source_pair: dict,
    population: dict,
    fiberwise_tree_law: dict | None = None,
    tau2_obstruction: dict | None = None,
) -> dict:
    source_slots = dict(source_pair["source_pair"])
    selected_point = dict(population.get("selected_population_point", {}))
    if not selected_point:
        raise ValueError("selected population point is required")

    alpha_y0 = float(source_slots["alphaY_mz"])
    alpha2_0 = float(source_slots["alpha2_mz"])
    v_value = float(source_slots["v_inherited"])
    beta_ew = (alpha2_0 - alpha_y0) / (alpha_y0 + alpha2_0)
    selected_tau_2 = float(selected_point["tau_2"])
    fiberwise_tree_law = fiberwise_tree_law or {}
    fiber_law_closed = fiberwise_tree_law.get("status") == "closed_smaller_primitive"
    if fiber_law_closed:
        depends_on_object = fiberwise_tree_law.get("object_id")
        next_residual = fiberwise_tree_law.get("next_single_residual_object", "tau2_tree_exact")
        tauy_formula = fiberwise_tree_law["tauY_formula"]
        sigma_formula = fiberwise_tree_law["sigma_formula"]
        eta_formula = fiberwise_tree_law["eta_formula"]
        proof_status = "one_scalar_beyond_fiberwise_population_tree_law_is_necessary_and_sufficient_for_exact_wz"
        status = "open_waiting_tau2_tree_exact"
    else:
        depends_on_object = "EWSinglePostTransportTreeIdentity_D10"
        next_residual = "EWSinglePostTransportTreeIdentity_D10"
        tauy_formula = "Phi_EW_tree(tau2_tree_exact)"
        sigma_formula = "0.5 * (Phi_EW_tree(tau2_tree_exact) + tau2_tree_exact)"
        eta_formula = "0.5 * (tau2_tree_exact - Phi_EW_tree(tau2_tree_exact))"
        proof_status = "one_scalar_beyond_single_tree_identity_is_necessary_and_sufficient_for_exact_wz"
        status = "open_waiting_single_tree_identity"

    if tau2_obstruction and tau2_obstruction.get("status") == "closed_smaller_primitive":
        proof_status = tau2_obstruction["proof_status"]
        status = "open_current_carrier_insufficient"
        next_residual = tau2_obstruction.get("next_single_residual_object", "delta_n_tree_exact")

    return {
        "artifact": "oph_d10_ew_exact_wz_coordinate_beyond_single_tree_identity",
        "generated_utc": _timestamp(),
        "object_id": "EWExactWZCoordinateBeyondSingleTreeIdentity_D10",
        "status": status,
        "proof_status": proof_status,
        "completion_kind": "strictly_smaller_coordinate_beneath_exact_wz_completion",
        "depends_on_object": depends_on_object,
        "fiberwise_population_tree_law_artifact": fiberwise_tree_law.get("artifact"),
        "fiberwise_population_tree_law_status": fiberwise_tree_law.get("status"),
        "current_carrier_obstruction_artifact": tau2_obstruction.get("artifact") if tau2_obstruction else None,
        "direct_tau2_emission_blocked": bool(tau2_obstruction and tau2_obstruction.get("status") == "closed_smaller_primitive"),
        "minimal_extra_scalar_or_invariant": (
            tau2_obstruction.get("minimal_extra_scalar_or_invariant", {}).get("symbol")
            if tau2_obstruction
            else None
        ),
        "population_evaluator_artifact": population.get("artifact"),
        "source_transport_pair_artifact": source_pair.get("artifact"),
        "coordinate_symbol": "tau2_tree_exact",
        "coordinate_kind": "charged_leg_post_transport_mass_coordinate",
        "tau2_tree_exact": None,
        "tauY_from_single_tree_identity_formula": tauy_formula,
        "u_EW_formula": "1 + tau2_tree_exact",
        "n_EW_formula": "1 + (alphaY_mz * tau_Y + alpha2_mz * tau2_tree_exact) / (alphaY_mz + alpha2_mz)",
        "sigma_from_coordinate_formula": sigma_formula,
        "eta_from_coordinate_formula": eta_formula,
        "eta_from_u_n_formula": "(u_EW - n_EW) / (1 - beta_EW)",
        "sigma_from_u_n_formula": "(n_EW - 1 - beta_EW * (u_EW - 1)) / (1 - beta_EW)",
        "MW_formula": "v_inherited * sqrt(pi * alpha2_mz * (1 + tau2_tree_exact))",
        "MZ_formula": "v_inherited * sqrt(pi * (alphaY_mz * (1 + tau_Y) + alpha2_mz * (1 + tau2_tree_exact)))",
        "selected_current_tau2": selected_tau_2,
        "carrier_basis_scalar": {
            "beta_EW": beta_ew,
            "beta_EW_formula": "(alpha2_mz - alphaY_mz) / (alpha2_mz + alphaY_mz)",
            "alphaY_mz": alpha_y0,
            "alpha2_mz": alpha2_0,
            "v_inherited": v_value,
        },
        "minimality_certificate": {
            "current_selected_tau2": selected_tau_2,
            "zero_extra_scalars_fail": True,
            "why_zero_extra_scalars_fail": "any downstream object preserving tau2=0 leaves MW fixed at the current public D10 value",
            "one_scalar_suffices_after_single_tree_identity": not bool(tau2_obstruction and tau2_obstruction.get("status") == "closed_smaller_primitive"),
            "why_one_scalar_suffices_after_single_tree_identity": (
                "the one-variable tree law determines tauY from tau2, so one scalar fixes both MW and MZ"
                if not (tau2_obstruction and tau2_obstruction.get("status") == "closed_smaller_primitive")
                else "the current one-variable carrier moves W and Z with the same local sign, so one extra neutral scalar is required beyond tau2_tree_exact"
            ),
        },
        "next_residual_object_if_open": next_residual,
        "notes": [
            "The selected current carrier closes the split exact law but freezes the charged-leg factor at tau2 = 0.",
            "The emitted fiberwise tree law isolates tau2_tree_exact as the charged-leg mover on the current carrier.",
            "The current-carrier obstruction shows that exact W/Z closure still needs one extra neutral scalar delta_n_tree_exact beyond direct tau2 emission.",
            "Once tau2_tree_exact and delta_n_tree_exact are both emitted, W and Z move independently on the exact mass-pair surface.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 exact-W/Z coordinate shell beyond the unsplit tree identity.")
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--population", default=str(DEFAULT_POPULATION))
    parser.add_argument("--fiberwise-tree-law", default=str(DEFAULT_FIBERWISE_TREE_LAW))
    parser.add_argument("--tau2-obstruction", default=str(DEFAULT_TAU2_OBSTRUCTION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_pair = json.loads(Path(args.source_pair).read_text(encoding="utf-8"))
    population = json.loads(Path(args.population).read_text(encoding="utf-8"))
    fiberwise_tree_law_path = Path(args.fiberwise_tree_law)
    fiberwise_tree_law = (
        json.loads(fiberwise_tree_law_path.read_text(encoding="utf-8"))
        if fiberwise_tree_law_path.exists()
        else None
    )
    tau2_obstruction_path = Path(args.tau2_obstruction)
    tau2_obstruction = (
        json.loads(tau2_obstruction_path.read_text(encoding="utf-8"))
        if tau2_obstruction_path.exists()
        else None
    )
    artifact = build_artifact(source_pair, population, fiberwise_tree_law, tau2_obstruction)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
