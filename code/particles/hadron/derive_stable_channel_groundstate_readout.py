#!/usr/bin/env python3
"""Aggregate per-ensemble stable-channel outputs into hadron mass readouts.

Chain role: gather the evaluated `pi_iso` and `N_iso` sequence families and
expose the first hadron mass/ratio readout surface once windows populate.

Mathematics: forward-window aggregation, effective-mass family bookkeeping, and
ensemble-by-ensemble ground-state/ratio extraction.

OPH-derived inputs: the seeded unquenched family together with the stable
channel population/evaluation artifacts and baryon contraction plan.

Output: the stable-channel ground-state readout that can feed public hadron rows
once convergence certificates exist.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AUDIT = ROOT / "particles" / "runs" / "hadron" / "current_hadron_lane_audit.json"
DEFAULT_CONTRACTION_PLAN = ROOT / "particles" / "runs" / "hadron" / "proton_contraction_plan.json"
DEFAULT_FULL_UNQUENCHED = ROOT / "particles" / "runs" / "hadron" / "full_unquenched_correlator.json"
DEFAULT_SEQUENCE_POPULATION = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_population.json"
DEFAULT_SEQUENCE_EVALUATION = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_evaluation.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "stable_channel_groundstate_readout.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _series(value: object) -> list[float]:
    if not isinstance(value, list):
        return []
    return [float(item) for item in value]


def _per_ensemble_channel_family(
    channel_key: str,
    ensembles: list[dict[str, object]],
) -> list[dict[str, object]]:
    family: list[dict[str, object]] = []
    for ensemble in ensembles:
        payload = ensemble.get(channel_key) or {}
        entry = {
            "ensemble_id": ensemble.get("ensemble_id"),
            "family_index": ensemble.get("family_index"),
            "T": ensemble.get("T"),
            "aLambda_msbar3": ensemble.get("aLambda_msbar3"),
            "corr_t": payload.get("corr_t", []),
            "am_eff_t": payload.get("am_eff_t", []),
            "forward_window_t": payload.get("forward_window_t", []),
            "forward_window_runs": payload.get("forward_window_runs", []),
            "selected_forward_window": payload.get("selected_forward_window", []),
            "selected_forward_window_cardinality": payload.get("selected_forward_window_cardinality", 0),
            "forward_window_limit_exists": payload.get("forward_window_limit_exists", False),
            "tail_drop_t": payload.get("tail_drop_t", []),
            "tail_drop_t_stderr": payload.get("tail_drop_t_stderr", []),
            "log_convexity_residual_t": payload.get("log_convexity_residual_t", []),
            "log_convexity_residual_t_stderr": payload.get("log_convexity_residual_t_stderr", []),
            "mirror_tail_indicator_t": payload.get("mirror_tail_indicator_t", []),
            "mirror_tail_indicator_t_stderr": payload.get("mirror_tail_indicator_t_stderr", []),
            "plateau_noise_floor_t": payload.get("plateau_noise_floor_t", []),
            "log_convexity_certified_t": payload.get("log_convexity_certified_t", []),
            "monotone_tail_t": payload.get("monotone_tail_t", []),
            "plateau_flat_t": payload.get("plateau_flat_t", []),
            "mirror_suppressed_t": payload.get("mirror_suppressed_t", []),
            "forward_certificate_t": payload.get("forward_certificate_t", []),
            "am_ground": payload.get("am_ground_candidate", payload.get("am_ground")),
            "am_ground_err": payload.get("am_ground_candidate_err"),
            "ratio_to_lambda_msbar3": payload.get("ratio_to_lambda_msbar3_candidate", payload.get("ratio_to_lambda_msbar3")),
            "mass_gev": payload.get("mass_gev_candidate", payload.get("mass_gev")),
            "convergence_status": payload.get("convergence_status", payload.get("sequence_status")),
        }
        if channel_key == "N_iso":
            entry.update(
                {
                    "corr_direct_t": payload.get("corr_direct_t", []),
                    "corr_exchange_t": payload.get("corr_exchange_t", []),
                    "corr_sign_t": payload.get("corr_sign_t", []),
                    "sign_stable_t": payload.get("sign_stable_t", []),
                    "direct_minus_exchange_residual_t": payload.get("direct_minus_exchange_residual_t", []),
                    "direct_minus_exchange_consistent_t": payload.get("direct_minus_exchange_consistent_t", []),
                }
            )
        family.append(entry)
    return family


def build_artifact(
    audit: dict[str, object],
    contraction_plan: dict[str, object] | None = None,
    full_unquenched: dict[str, object] | None = None,
    sequence_population: dict[str, object] | None = None,
    sequence_evaluation: dict[str, object] | None = None,
) -> dict[str, object]:
    contraction_plan = contraction_plan or {}
    full_unquenched = full_unquenched or {}
    sequence_population = sequence_population or {}
    sequence_evaluation = sequence_evaluation or {}
    lane = audit.get("pipeline_classification") or {}
    full_channels = full_unquenched.get("channel_payloads") or {}
    sequence_ensembles = sequence_population.get("ensemble_sequences") or []
    evaluation_ensembles = sequence_evaluation.get("ensemble_evaluations") or []
    active_sequence = sequence_ensembles[0] if sequence_ensembles else {}
    active_evaluation = evaluation_ensembles[0] if evaluation_ensembles else {}
    pi_payload = active_evaluation.get("pi_iso") or active_sequence.get("pi_iso") or full_channels.get("pi_iso") or {}
    n_payload = active_evaluation.get("N_iso") or active_sequence.get("N_iso") or full_channels.get("N_iso") or {}
    qcd_inputs = full_unquenched.get("qcd_inputs") or {}
    ensemble_family_inputs = evaluation_ensembles if evaluation_ensembles else sequence_ensembles

    pi_corr_t = _series(pi_payload.get("corr_t"))
    pi_am_eff_t = _series(pi_payload.get("am_eff_t"))
    n_corr_direct_t = _series(n_payload.get("corr_direct_t"))
    n_corr_exchange_t = _series(n_payload.get("corr_exchange_t"))
    n_corr_t = _series(n_payload.get("corr_t"))
    n_am_eff_t = _series(n_payload.get("am_eff_t"))
    contraction_status = contraction_plan.get("status")
    contraction_rule_closed = contraction_status in {"closed", "formula_closed"}

    raw_arrays_present = any(
        seq
        for seq in [pi_corr_t, n_corr_direct_t, n_corr_exchange_t, n_corr_t]
    )
    am_eff_present = any(seq for seq in [pi_am_eff_t, n_am_eff_t])
    why_not = []
    if not raw_arrays_present:
        why_not.append("The seeded ensemble family has not emitted stable-channel correlator sequences yet.")
    else:
        why_not.append("Stable-channel correlator sequences are present but not yet converged.")
    if not contraction_rule_closed:
        why_not.append("The direct-minus-exchange nucleon contraction rule is still open.")
    if not am_eff_present:
        why_not.append("Long-time effective-mass sequences are still unset.")

    return {
        "artifact": "oph_hadron_stable_channel_groundstate_readout",
        "generated_utc": _timestamp(),
        "proof_status": "candidate_only",
        "theorem_candidate": "StableChannelForwardWindowConvergence",
        "source_lane_status": lane.get("lane_status"),
        "upstream": {
            "hadron_lane_audit_artifact": str(DEFAULT_AUDIT),
            "full_unquenched_correlator_artifact": full_unquenched.get("artifact"),
            "full_unquenched_correlator_status": full_unquenched.get("status"),
            "stable_channel_sequence_population_artifact": sequence_population.get("artifact"),
            "stable_channel_sequence_population_status": sequence_population.get("status"),
            "stable_channel_sequence_evaluation_artifact": sequence_evaluation.get("artifact"),
            "stable_channel_sequence_evaluation_status": sequence_evaluation.get("status"),
            "ensemble_population_status": ((full_unquenched.get("population_contract") or {}).get("status")),
            "lambda_msbar_3_gev": qcd_inputs.get("Lambda_MSbar_3_gev"),
        },
        "ensemble_inputs": {
            "sea_flavor_target": full_unquenched.get("sea_flavor_target"),
            "light_mass_rule": full_unquenched.get("light_mass_rule"),
            "ensemble_family": full_unquenched.get("ensemble_family"),
        },
        "isospin_symmetric_inputs": {
            "qed_off": True,
            "sea_flavors": full_unquenched.get("sea_flavor_target"),
            "light_mass_rule": "m_l := (m_u + m_d)/2",
        },
        "operator_fixing_rule": "lowest_local_color_singlet_by_channel_quantum_numbers",
        "readout_fixing_rule": "full_correlator_sequence_plus_forward_window_log_ratio_limit",
        "data_availability": {
            "raw_correlator_arrays_present": raw_arrays_present,
            "effective_mass_sequences_present": am_eff_present,
            "full_baryon_contractions_present": contraction_rule_closed,
            "full_baryon_contractions_artifact": contraction_plan.get("artifact"),
            "full_baryon_contractions_status": contraction_status,
            "full_unquenched_correlator_artifact": full_unquenched.get("artifact"),
            "full_unquenched_correlator_status": full_unquenched.get("status", "missing"),
            "why_not": why_not,
        },
        "channels": {
            "pi_iso": {
                "quantum_numbers": {"I": 1, "JPC": "0-+"},
                "operator_id": pi_payload.get("operator_id", "P_local_isovector"),
                "operator_formula": pi_payload.get("operator_formula", "sum_x qbar gamma5 tau^a q"),
                "contraction_kind": "connected_isovector",
                "corr_t": pi_corr_t,
                "am_eff_t": pi_am_eff_t,
                "am_ground": pi_payload.get("am_ground"),
                "ratio_to_lambda_msbar3": pi_payload.get("ratio_to_lambda_msbar3"),
                "mass_gev": pi_payload.get("mass_gev"),
                "forward_window_t": pi_payload.get("forward_window_t", []),
                "forward_window_runs": pi_payload.get("forward_window_runs", []),
                "selected_forward_window": pi_payload.get("selected_forward_window", []),
                "selected_forward_window_cardinality": pi_payload.get("selected_forward_window_cardinality", 0),
                "forward_window_limit_exists": pi_payload.get("forward_window_limit_exists", False),
                "tail_drop_t": pi_payload.get("tail_drop_t", []),
                "tail_drop_t_stderr": pi_payload.get("tail_drop_t_stderr", []),
                "log_convexity_residual_t": pi_payload.get("log_convexity_residual_t", []),
                "log_convexity_residual_t_stderr": pi_payload.get("log_convexity_residual_t_stderr", []),
                "mirror_tail_indicator_t": pi_payload.get("mirror_tail_indicator_t", []),
                "mirror_tail_indicator_t_stderr": pi_payload.get("mirror_tail_indicator_t_stderr", []),
                "plateau_noise_floor_t": pi_payload.get("plateau_noise_floor_t", []),
                "log_convexity_certified_t": pi_payload.get("log_convexity_certified_t", []),
                "monotone_tail_t": pi_payload.get("monotone_tail_t", []),
                "plateau_flat_t": pi_payload.get("plateau_flat_t", []),
                "mirror_suppressed_t": pi_payload.get("mirror_suppressed_t", []),
                "forward_certificate_t": pi_payload.get("forward_certificate_t", []),
                "per_ensemble": _per_ensemble_channel_family("pi_iso", ensemble_family_inputs),
                "groundstate_limit_formula": "a m_pi_iso := lim_{t up to tail(W_pi)} log(C_pi(t) / C_pi(t+1))",
                "ratio_to_lambda_msbar3_formula": "R_pi_iso = (a m_pi_iso) / (a Lambda_MSbar_3)",
                "algebraic_mass_rule": "mass_gev = ratio_to_lambda_msbar3 * Lambda_MSbar_3_gev",
                "promoted_channel_fields": [
                    "corr_t",
                    "am_eff_t",
                    "am_ground",
                    "ratio_to_lambda_msbar3",
                    "mass_gev",
                ],
            },
            "N_iso": {
                "quantum_numbers": {"I": 0.5, "JP": "1/2+"},
                "operator_id": n_payload.get("operator_id", "N_local_Cgamma5"),
                "operator_formula": n_payload.get("operator_formula", "sum_x eps^{abc}(u^T C gamma5 d)u"),
                "parity_projector": "(1+gamma_0)/2",
                "contraction_kind": (
                    "direct_minus_exchange_formula_closed"
                    if contraction_rule_closed
                    else "full_baryon_contractions_required"
                ),
                "corr_direct_t": n_corr_direct_t,
                "corr_exchange_t": n_corr_exchange_t,
                "corr_t": n_corr_t,
                "am_eff_t": n_am_eff_t,
                "am_ground": n_payload.get("am_ground"),
                "ratio_to_lambda_msbar3": n_payload.get("ratio_to_lambda_msbar3"),
                "mass_gev": n_payload.get("mass_gev"),
                "corr_sign_t": n_payload.get("corr_sign_t", []),
                "sign_stable_t": n_payload.get("sign_stable_t", []),
                "forward_window_t": n_payload.get("forward_window_t", []),
                "forward_window_runs": n_payload.get("forward_window_runs", []),
                "selected_forward_window": n_payload.get("selected_forward_window", []),
                "selected_forward_window_cardinality": n_payload.get("selected_forward_window_cardinality", 0),
                "forward_window_limit_exists": n_payload.get("forward_window_limit_exists", False),
                "tail_drop_t": n_payload.get("tail_drop_t", []),
                "tail_drop_t_stderr": n_payload.get("tail_drop_t_stderr", []),
                "log_convexity_residual_t": n_payload.get("log_convexity_residual_t", []),
                "log_convexity_residual_t_stderr": n_payload.get("log_convexity_residual_t_stderr", []),
                "mirror_tail_indicator_t": n_payload.get("mirror_tail_indicator_t", []),
                "mirror_tail_indicator_t_stderr": n_payload.get("mirror_tail_indicator_t_stderr", []),
                "plateau_noise_floor_t": n_payload.get("plateau_noise_floor_t", []),
                "log_convexity_certified_t": n_payload.get("log_convexity_certified_t", []),
                "monotone_tail_t": n_payload.get("monotone_tail_t", []),
                "plateau_flat_t": n_payload.get("plateau_flat_t", []),
                "mirror_suppressed_t": n_payload.get("mirror_suppressed_t", []),
                "forward_certificate_t": n_payload.get("forward_certificate_t", []),
                "direct_minus_exchange_residual_t": n_payload.get("direct_minus_exchange_residual_t", []),
                "direct_minus_exchange_consistent_t": n_payload.get("direct_minus_exchange_consistent_t", []),
                "per_ensemble": _per_ensemble_channel_family("N_iso", ensemble_family_inputs),
                "groundstate_limit_formula": "a m_N_iso := lim_{t up to tail(W_N)} log(|C_N(t)| / |C_N(t+1)|)",
                "ratio_to_lambda_msbar3_formula": "R_N_iso = (a m_N_iso) / (a Lambda_MSbar_3)",
                "algebraic_mass_rule": "mass_gev = ratio_to_lambda_msbar3 * Lambda_MSbar_3_gev",
                "promoted_channel_fields": [
                    "corr_direct_t",
                    "corr_exchange_t",
                    "corr_t",
                    "am_eff_t",
                    "am_ground",
                    "ratio_to_lambda_msbar3",
                    "mass_gev",
                ],
                "full_baryon_contractions_artifact": contraction_plan.get("artifact"),
                "full_baryon_contractions_status": contraction_status,
            },
        },
        "rho_extension": {
            "status": "separate_scattering_problem",
            "next_artifact": "oph_hadron_rho_scattering_readout",
        },
        "promotion_gate": audit.get("promotion_gate", "StableChannelForwardWindowConvergence"),
        "promotion_gate_payload": {
            "name": "StableChannelForwardWindowConvergence",
            "scope": ["pi_iso", "N_iso"],
            "excludes": ["rho_scattering"],
            "requires": [
                "full_unquenched_correlator",
                "full_baryon_contractions",
                "stable_channel_sequence_evaluation",
                "forward_window_log_ratio_limit_exists",
            ],
        },
        "minimal_closure_frontier": audit.get("minimal_closure_frontier", []),
        "notes": [
            "This artifact tracks the stable-channel readout boundary used by the hadron pipeline.",
            "Stable-channel masses remain unset until the unquenched producer emits correlator sequences and the long-time limits converge.",
            "The nucleon branch also requires the full baryon contraction payload on the same boundary.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the stable-channel hadron ground-state readout artifact.")
    parser.add_argument("--audit", default=str(DEFAULT_AUDIT))
    parser.add_argument("--contraction-plan", default=str(DEFAULT_CONTRACTION_PLAN))
    parser.add_argument("--full-unquenched", default=str(DEFAULT_FULL_UNQUENCHED))
    parser.add_argument("--sequence-population", default=str(DEFAULT_SEQUENCE_POPULATION))
    parser.add_argument("--sequence-evaluation", default=str(DEFAULT_SEQUENCE_EVALUATION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    audit_path = Path(args.audit)
    contraction_plan_path = Path(args.contraction_plan)
    full_unquenched_path = Path(args.full_unquenched)
    sequence_population_path = Path(args.sequence_population)
    sequence_evaluation_path = Path(args.sequence_evaluation)
    out_path = Path(args.output)

    contraction_plan = _load(contraction_plan_path) if contraction_plan_path.exists() else None
    full_unquenched = _load(full_unquenched_path) if full_unquenched_path.exists() else None
    sequence_population = _load(sequence_population_path) if sequence_population_path.exists() else None
    sequence_evaluation = _load(sequence_evaluation_path) if sequence_evaluation_path.exists() else None
    artifact = build_artifact(_load(audit_path), contraction_plan, full_unquenched, sequence_population, sequence_evaluation)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
