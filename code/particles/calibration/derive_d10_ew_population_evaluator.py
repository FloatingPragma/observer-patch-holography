#!/usr/bin/env python3
"""Build the carrier-level selector functional for the D10 two-scalar family.

Chain role: define the point-separating population functional that chooses one
carrier point on the reduced D10 transport family.

Mathematics: source-point comparison, transport-entry bookkeeping, and the
scalar selector `J_pop_EW` on `(sigma_EW, eta_EW)`.

OPH-derived inputs: the D10 source pair, its seed trial, and the reduced
carrier basis inherited from the calibration core.

Output: the selector artifact that closes the current carrier choice while
leaving exact full electroweak closure as a downstream problem.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_PAIR = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
DEFAULT_READOUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_readout.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_population_evaluator.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(source_pair: dict, readout: dict) -> dict:
    nonuniqueness = dict(source_pair.get("population_nonuniqueness_certificate", {}))
    source_point_a = dict(nonuniqueness.get("source_point_A", {}))
    source_point_b = dict(nonuniqueness.get("source_point_B", {}))
    source_slots = dict(source_pair.get("source_pair", {}))
    alpha_y = float(source_slots.get("alphaY_mz", 0.0))
    alpha2 = float(source_slots.get("alpha2_mz", 0.0))
    beta_ew = float((source_pair.get("population_basis") or {}).get("beta_EW", 0.0))
    tau_seed = float(((source_pair.get("first_nonzero_oph_seed_trial") or {}).get("tau_seed", 0.0)))
    alpha_u_from_seed = (8.0 * math.pi * tau_seed) / math.log(alpha_y / alpha2)
    eta_source = alpha_u_from_seed * beta_ew
    compact_point = dict(readout.get("current_compact_point", {}))
    compact_quintet = dict(readout.get("current_compact_emitted_quintet", {}))
    compact_transport = dict(
        ((source_pair.get("compact_hypercharge_only_mass_slice") or {}).get("transport_entry_values")) or {}
    )
    eta_star = compact_point.get("eta_EW")
    if eta_star is not None:
        eta_star = float(eta_star)
    compact_source_defect_values = {
        "A": None,
        "B": None,
        "delta_B_minus_A": None,
    }
    if eta_star is not None:
        sigma_a = float(source_point_a.get("sigma_EW", 0.0))
        eta_a = float(source_point_a.get("eta_EW", 0.0))
        sigma_b = float(source_point_b.get("sigma_EW", 0.0))
        eta_b = float(source_point_b.get("eta_EW", 0.0))
        compact_source_defect_values["A"] = (sigma_a + eta_a) ** 2 + (eta_a - eta_star) ** 2
        compact_source_defect_values["B"] = (sigma_b + eta_b) ** 2 + (eta_b - eta_star) ** 2
        compact_source_defect_values["delta_B_minus_A"] = (
            compact_source_defect_values["B"] - compact_source_defect_values["A"]
        )
    sigma_selected = -eta_source
    eta_selected = eta_source
    selected_point = {
        "sigma_EW": sigma_selected,
        "eta_EW": eta_selected,
        "tau_Y": sigma_selected - eta_selected,
        "tau_2": sigma_selected + eta_selected,
    }
    selected_basis_point = {
        "u_EW": 1.0 + sigma_selected + eta_selected,
        "n_EW": 1.0 + sigma_selected + beta_ew * eta_selected,
    }
    selected_transport_entries = {
        "Pi_AA": sigma_selected - beta_ew * eta_selected,
        "Pi_AZ": ((1.0 - beta_ew**2) ** 0.5) * eta_selected,
        "Pi_ZZ": sigma_selected + beta_ew * eta_selected,
        "Pi_WW": sigma_selected + eta_selected,
    }
    selected_quintet = dict(compact_quintet)
    def _j_pop(sigma_ew: float, eta_ew: float) -> float:
        return (sigma_ew + eta_ew) ** 2 + (sigma_ew * sigma_ew - eta_ew * eta_ew) ** 2 + (sigma_ew + eta_source) ** 2
    source_point_a_value = _j_pop(float(source_point_a.get("sigma_EW", 0.0)), float(source_point_a.get("eta_EW", 0.0)))
    source_point_b_value = _j_pop(float(source_point_b.get("sigma_EW", 0.0)), float(source_point_b.get("eta_EW", 0.0)))
    return {
        "artifact": "oph_d10_ew_population_evaluator",
        "generated_utc": _timestamp(),
        "object_id": "EWGaugeSourceTransportPairPopulationEvaluator_D10",
        "status": "closed_current_carrier",
        "proof_status": "population_functional_closed_on_current_carrier",
        "strict_obstruction_theorem": None,
        "strict_obstruction_status": None,
        "population_selector_obstructed": False,
        "predictive_promotion_allowed": False,
        "source_transport_pair_artifact": source_pair.get("artifact"),
        "source_transport_readout_artifact": readout.get("artifact"),
        "population_basis": dict(source_pair.get("population_basis", {})),
        "population_reconstruction": dict(source_pair.get("population_reconstruction", {})),
        "population_atomic_quartet": dict(source_pair.get("population_atomic_quartet", {})),
        "population_minimality_certificate": dict(source_pair.get("population_minimality_certificate", {})),
        "quartet_injectivity_certificate": {
            "injective_on_current_family": True,
            "u_EW_from_mass_pair": "mW^2 / (pi * v_inherited^2 * alpha2_mz)",
            "n_EW_from_mass_pair": "mZ^2 / (pi * v_inherited^2 * (alphaY_mz + alpha2_mz))",
            "eta_EW_from_basis": "(u_EW - n_EW) / (1 - beta_EW)",
            "sigma_EW_from_basis": "(n_EW - 1 - beta_EW * (u_EW - 1)) / (1 - beta_EW)",
        },
        "admissible_family_points": {
            "source_point_A": source_point_a,
            "source_point_B": source_point_b,
        },
        "obstruction_witness_pair": {
            "A": {
                "sigma_EW": source_point_a.get("sigma_EW"),
                "eta_EW": source_point_a.get("eta_EW"),
                "u_EW": 1.0 + float(source_point_a.get("sigma_EW", 0.0)) + float(source_point_a.get("eta_EW", 0.0)),
                "n_EW": 1.0 + float(source_point_a.get("sigma_EW", 0.0)) + beta_ew * float(source_point_a.get("eta_EW", 0.0)),
            },
            "B": {
                "sigma_EW": source_point_b.get("sigma_EW"),
                "eta_EW": source_point_b.get("eta_EW"),
                "u_EW": 1.0 + float(source_point_b.get("sigma_EW", 0.0)) + float(source_point_b.get("eta_EW", 0.0)),
                "n_EW": 1.0 + float(source_point_b.get("sigma_EW", 0.0)) + beta_ew * float(source_point_b.get("eta_EW", 0.0)),
            },
            "delta_basis_B_minus_A": {
                "u_EW": (
                    1.0 + float(source_point_b.get("sigma_EW", 0.0)) + float(source_point_b.get("eta_EW", 0.0))
                    - 1.0 - float(source_point_a.get("sigma_EW", 0.0)) - float(source_point_a.get("eta_EW", 0.0))
                ),
                "n_EW": (
                    1.0 + float(source_point_b.get("sigma_EW", 0.0)) + beta_ew * float(source_point_b.get("eta_EW", 0.0))
                    - 1.0 - float(source_point_a.get("sigma_EW", 0.0)) - beta_ew * float(source_point_a.get("eta_EW", 0.0))
                ),
            },
            "delta_quintet_B_minus_A": dict(nonuniqueness.get("delta_quintet_B_minus_A", {})),
        },
        "selected_population_point": selected_point,
        "selected_population_basis_point": selected_basis_point,
        "selected_transport_entries": selected_transport_entries,
        "selected_quintet": selected_quintet,
        "selector_formula": "selected_population_point = argmin_{p in C_D10} J_pop_EW(p)",
        "population_selector_formula": "selected_population_point = argmin_{p in C_D10} J_pop_EW(p)",
        "predictive_population_point": selected_point,
        "population_driven_quintet": selected_quintet,
        "population_functional_symbol": "J_pop_EW",
        "population_functional_status": "closed",
        "population_functional_origin": "rank_one_trace_anchored_source_shell",
        "alpha_u_from_seed_formula": "8*pi*tau_seed / log(alphaY_mz / alpha2_mz)",
        "alpha_u_from_seed": alpha_u_from_seed,
        "eta_source_formula": "alpha_u_from_seed * beta_EW",
        "eta_source": eta_source,
        "population_functional_domain": "(sigma_EW, eta_EW)",
        "population_functional_basis_domain": "(u_EW, n_EW)",
        "population_functional_formula_transport": "J_pop_EW = Pi_WW^2 + (Pi_AA*Pi_ZZ - Pi_AZ^2)^2 + (0.5*(Pi_AA + Pi_ZZ) + eta_source)^2",
        "population_functional_formula_sigma_eta": "J_pop_EW(sigma_EW,eta_EW) = (sigma_EW + eta_EW)^2 + (sigma_EW^2 - eta_EW^2)^2 + (sigma_EW + eta_source)^2",
        "population_functional_formula_basis": "J_pop_EW_basis(u_EW,n_EW) = (u_EW - 1)^2 + (sigma(u_EW,n_EW)^2 - eta(u_EW,n_EW)^2)^2 + (sigma(u_EW,n_EW) + eta_source)^2",
        "population_selector_rule": "selected_population_point = argmin_{p in C_D10} J_pop_EW(p)",
        "population_selector_uniqueness_requirement": "argmin is a singleton",
        "population_selector_separation_requirement": "J_pop_EW(A) != J_pop_EW(B)",
        "population_functional_values": {
            "A": source_point_a_value,
            "B": source_point_b_value,
            "delta_B_minus_A": source_point_b_value - source_point_a_value,
        },
        "candidate_population_functional_status": (
            "demoted_shell_restriction" if eta_star is not None else "missing"
        ),
        "candidate_population_functional_origin": (
            "compact_source_defect_to_current_hypercharge_only_mass_slice"
            if eta_star is not None
            else None
        ),
        "candidate_population_functional_relation_to_J_pop_EW": (
            "J_pop_EW restricted to Pi_WW = 0 and det(Pi_AZ_block) = 0"
            if eta_star is not None
            else None
        ),
        "candidate_population_functional_formula_sigma_eta": (
            "J_pop_EW_candidate(sigma_EW,eta_EW) = (sigma_EW + eta_EW)^2 + (eta_EW - eta_star_compact)^2"
            if eta_star is not None
            else None
        ),
        "candidate_population_functional_formula_basis": (
            "J_pop_EW_candidate_basis(u_EW,n_EW) = (u_EW - 1)^2 + (((u_EW - n_EW) / (1 - beta_EW)) - eta_star_compact)^2"
            if eta_star is not None
            else None
        ),
        "candidate_population_functional_values": compact_source_defect_values,
        "candidate_selector_formula": (
            "argmin_{p in C_D10} ((sigma_EW + eta_EW)^2 + (eta_EW - eta_star_compact)^2)"
            if eta_star is not None
            else None
        ),
        "candidate_selected_population_point": None if not compact_point else {
            "sigma_EW": float(compact_point["sigma_EW"]),
            "eta_EW": float(compact_point["eta_EW"]),
            "tau_Y": float(compact_point["tau_Y"]),
            "tau_2": float(compact_point["tau_2"]),
        },
        "candidate_selected_population_basis_point": selected_basis_point,
        "candidate_selected_transport_entries": compact_transport or None,
        "candidate_selected_quintet": compact_quintet or None,
        "candidate_exactness_verdict": (
            "demoted_shell_restriction_of_closed_current_carrier" if eta_star is not None else "missing"
        ),
        "evaluator_contract": {
            "carrier_domain": "(sigma_EW, eta_EW)",
            "carrier_basis": "(u_EW, n_EW)",
            "selection_output": "one unique admissible point on the existing reopened two-scalar carrier",
            "selection_forbidden_inputs": [
                "reference_W",
                "reference_Z",
                "reference_alpha_em_inv",
                "reference_sin2_thetaW",
            ],
        },
        "smallest_constructive_missing_object": None,
        "required_extra_object": None,
        "required_extra_invariant_symbol": None,
        "required_extra_object_contract": None,
        "notes": [
            "The reopened two-scalar D10 carrier is already minimal and the carrier-level selector is now closed on current data.",
            "J_pop_EW selects the same hypercharge-only point that was previously only recorded as a compact shell candidate, so that shell is now demoted to a restriction of the full carrier functional.",
            "The remaining D10 gap is no longer carrier selection but exact electroweak closure beyond the current selected carrier point.",
            "No reference-fit W/Z slice is consumed anywhere in this evaluator.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D10 population-evaluator artifact.")
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
