#!/usr/bin/env python3
"""Emit a split algebraic readout on the declared D10 carrier point.

Chain role: audit a split D10 readout beyond carrier selection by keeping the
declared transported point for the tree/chart pair and algebraically restoring
the neutral input through a compensating hypercharge leg.

Mathematics: selected-carrier transport, source-normalized hypercharge
compensation, and split mass/neutral readout bookkeeping.

Inputs: the declared-map D10 carrier point, electroweak chart slots, and the
selected-point transport couplings.  This artifact does not prove that the
source selects the map or that any displayed tree/chart coordinate is a pole.

Output: an exact algebraic closure of this split declared-map bookkeeping,
together with the stronger residual object that would unify the split readout
into one post-transport tree identity.  Exactness here is not physical or
predictive closure.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_READOUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_readout.json"
DEFAULT_POPULATION = ROOT / "particles" / "runs" / "calibration" / "d10_ew_population_evaluator.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_closure_beyond_current_carrier.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_pair: dict, readout: dict, population: dict) -> dict:
    selected_point = dict(population.get("selected_population_point", {}))
    if not selected_point:
        raise ValueError("selected D10 carrier point is required")
    source_slots = dict(source_pair.get("source_pair", {}))
    alpha_y0 = float(source_slots["alphaY_mz"])
    alpha2_0 = float(source_slots["alpha2_mz"])
    v_value = float(source_slots["v_inherited"])
    sigma_selected = float(selected_point["sigma_EW"])
    eta_selected = float(selected_point["eta_EW"])
    tau_y = float(selected_point["tau_Y"])
    tau_2 = float(selected_point["tau_2"])

    alpha_y_star = alpha_y0 * (1.0 + tau_y)
    alpha2_star = alpha2_0 * (1.0 + tau_2)
    if alpha_y_star == 0.0:
        raise ValueError("selected hypercharge transport vanished")

    chi_y = alpha_y0 / alpha_y_star
    alpha_y_exact = chi_y * alpha_y_star
    alpha2_exact = alpha2_star

    exact_outputs = {
        "MW_pole": v_value * (3.141592653589793 * alpha2_star) ** 0.5,
        "MZ_pole": v_value * (3.141592653589793 * (alpha_y_star + alpha2_star)) ** 0.5,
        "alpha_em_eff_inv": (alpha_y_exact + alpha2_exact) / (alpha_y_exact * alpha2_exact),
        "sin2w_eff": alpha_y_exact / (alpha_y_exact + alpha2_exact),
        "v_report": v_value,
    }

    return {
        "artifact": "oph_d10_ew_exact_closure_beyond_current_carrier",
        "generated_utc": _timestamp(),
        "object_id": "EWTransportExactClosureBeyondCurrentCarrier_D10",
        "status": "closed",
        "proof_status": "closed_on_existing_selected_carrier_via_source_normalized_hypercharge_readout",
        "epistemic_status": "exact_algebraic_identity_on_declared_map_only",
        "diagnostic_only": True,
        "source_selected": False,
        "physical_pole_claim": False,
        "predictive_promotion_allowed": False,
        "exactness_surface_kind": "split_mass_neutral_readout",
        "extra_invariant_required": False,
        "minimal_carrier_enlargement_required": False,
        "completion_kind": "selected_current_carrier_plus_source_normalized_hypercharge_readout",
        "source_transport_pair_artifact": source_pair.get("artifact"),
        "source_transport_readout_artifact": readout.get("artifact"),
        "population_evaluator_artifact": population.get("artifact"),
        "selected_population_point": {
            "sigma_EW": sigma_selected,
            "eta_EW": eta_selected,
            "tau_Y": tau_y,
            "tau_2": tau_2,
        },
        "source_slots": {
            "alphaY_mz": alpha_y0,
            "alpha2_mz": alpha2_0,
            "v_inherited": v_value,
        },
        "transported_mass_couplings": {
            "alphaY_star_selected": alpha_y_star,
            "alpha2_star_selected": alpha2_star,
        },
        "derived_readout_compensator": {
            "symbol": "chi_Y_EW",
            "kind": "derived_not_free",
            "formula": "1 / (1 + tau_Y_selected)",
            "equivalent_formulas": [
                "alphaY_mz / alphaY_star_selected",
                "1 / (1 + sigma_selected - eta_selected)",
            ],
            "value": chi_y,
        },
        "exact_readout_couplings": {
            "alphaY_readout_exact": alpha_y_exact,
            "alpha2_readout_exact": alpha2_exact,
        },
        "exact_forward_map": {
            "MW_pole": "v_inherited * sqrt(pi * alpha2_star_selected)",
            "MZ_pole": "v_inherited * sqrt(pi * (alphaY_star_selected + alpha2_star_selected))",
            "alpha_em_eff_inv": "(alphaY_readout_exact + alpha2_readout_exact) / (alphaY_readout_exact * alpha2_readout_exact)",
            "sin2w_eff": "alphaY_readout_exact / (alphaY_readout_exact + alpha2_readout_exact)",
        },
        "exact_outputs": exact_outputs,
        "minimality_certificate": {
            "selected_carrier_already_closed": bool(population.get("population_functional_status") == "closed"),
            "remaining_nonexact_outputs_under_old_readout": [
                "alpha_em_eff_inv",
                "sin2w_eff",
            ],
            "remaining_outputs_factor_through_single_neutral_scalar": True,
            "one_derived_hypercharge_compensator_suffices": True,
        },
        "stronger_residual_object": "EWSinglePostTransportTreeIdentity_D10",
        "notes": [
            "The declared carrier point supplies tree/chart W and Z coordinates; the source does not select this map and no physical pole is constructed.",
            "The neutral bookkeeping identity restores the original hypercharge input through a compensator defined as its exact algebraic inverse; this is not an independent source derivation.",
            "This closes split algebraic exactness on the declared carrier without enlarging it or fitting reference observables; it carries no prediction or physical-pole status.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 exact closure artifact beyond the selected carrier.")
    parser.add_argument("--source-pair", default=str(DEFAULT_SOURCE_PAIR))
    parser.add_argument("--readout", default=str(DEFAULT_READOUT))
    parser.add_argument("--population", default=str(DEFAULT_POPULATION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    source_pair = json.loads(Path(args.source_pair).read_text(encoding="utf-8"))
    readout = json.loads(Path(args.readout).read_text(encoding="utf-8"))
    population = json.loads(Path(args.population).read_text(encoding="utf-8"))
    artifact = build_artifact(source_pair, readout, population)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
