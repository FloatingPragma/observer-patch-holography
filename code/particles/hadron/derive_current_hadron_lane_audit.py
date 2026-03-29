#!/usr/bin/env python3
"""Summarize the live hadron frontier and the next missing constructive object.

Chain role: classify the active hadron path so the project can tell whether the
remaining blocker is ensemble seeding, cfg/source realization, convergence, or
resonance extraction.

Mathematics: rule-based pipeline classification only; this file does not derive
hadron masses itself.

OPH-derived inputs: the statuses of the active `/particles` hadron artifacts.

Output: the current hadron frontier report consumed by the public surface and
completion planning.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
DEFAULT_FULL_UNQUENCHED = ROOT / "particles" / "runs" / "hadron" / "full_unquenched_correlator.json"
DEFAULT_CONTRACTION_PLAN = ROOT / "particles" / "runs" / "hadron" / "proton_contraction_plan.json"
DEFAULT_SEQUENCE_POPULATION = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_population.json"
DEFAULT_CFG_SOURCE_PAYLOAD = ROOT / "particles" / "runs" / "hadron" / "stable_channel_cfg_source_measure_payload.json"
DEFAULT_SEQUENCE_EVALUATION = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_evaluation.json"
DEFAULT_SURROGATE_BRIDGE = ROOT / "particles" / "runs" / "hadron" / "hadron_surrogate_execution_bridge_status.json"
DEFAULT_GEOMETRY_SUMMARY = ROOT / "particles" / "runs" / "hadron" / "production_geometry_summary.json"
DEFAULT_CLOSURE_REPORT = ROOT / "particles" / "runs" / "hadron" / "hadron_production_closure_validation_report.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "current_hadron_lane_audit.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def classify_pipeline(
    full_unquenched: dict[str, Any],
    contraction_plan: dict[str, Any] | None,
    sequence_population: dict[str, Any] | None,
    cfg_source_payload: dict[str, Any] | None,
    sequence_evaluation: dict[str, Any] | None,
    surrogate_bridge: dict[str, Any] | None,
    geometry_summary: dict[str, Any] | None,
    closure_report: dict[str, Any] | None,
) -> dict[str, Any]:
    blockers: list[str] = []
    stable_channels_ready = bool((closure_report or {}).get("public_unsuppression_ready"))
    next_missing_object = full_unquenched.get("next_missing_object")
    population_status = ((full_unquenched.get("population_contract") or {}).get("status"))
    sequence_status = None if sequence_population is None else sequence_population.get("status")
    payload_status = None if cfg_source_payload is None else cfg_source_payload.get("status")
    evaluation_status = None if sequence_evaluation is None else sequence_evaluation.get("status")
    contraction_status = None if contraction_plan is None else contraction_plan.get("status")
    contraction_ready = contraction_status in {"closed", "formula_closed"}

    if stable_channels_ready:
        blockers.extend(["finite_volume_resonance_and_spectrum_readout"])
        effective_next_missing_object = "rho_resonance_extraction"
    elif evaluation_status == "awaiting_measure_evaluation" or payload_status == "law_closed_waiting_measure_realization":
        blockers.extend(
            [
                "stable_channel_sequence_evaluation",
                "stable_channel_cfg_source_measure_payload",
            ]
        )
        schedule = (cfg_source_payload or {}).get("support_realization_schedule") or {}
        runtime_receipt_emitted = cfg_source_payload is not None and (
            cfg_source_payload.get("runtime_receipt_artifact") == "runtime_schedule_receipt_N_therm_and_N_sep"
            or schedule.get("runtime_receipt_artifact") == "runtime_schedule_receipt_N_therm_and_N_sep"
        )
        closure_residual = None if closure_report is None else closure_report.get("smallest_live_residual_object")
        effective_next_missing_object = (
            str(closure_residual)
            if closure_residual
            else "backend_correlator_dump.production.json from real production RHMC/HMC execution on the theorem-emitted seeded family"
            if runtime_receipt_emitted
            else cfg_source_payload.get("smallest_constructive_missing_object")
            if cfg_source_payload is not None
            else "runtime_schedule_receipt_N_therm_and_N_sep"
        )
    elif sequence_status == "law_closed_waiting_measure_evaluation":
        blockers.append("stable_channel_sequence_evaluation")
        effective_next_missing_object = "stable_channel_sequence_evaluation"
    elif next_missing_object:
        blockers.append(str(next_missing_object))
        effective_next_missing_object = str(next_missing_object)
    else:
        blockers.append("stable_channel_sequence_population")
        effective_next_missing_object = "stable_channel_sequence_population"
    if not contraction_ready:
        blockers.append("full_baryon_contractions")
    if not stable_channels_ready:
        blockers.extend(
            [
                "StableChannelForwardWindowConvergence",
                "finite_volume_resonance_and_spectrum_readout",
            ]
        )

    return {
        "lane_status": (
            "stable_channel_public_ready"
            if stable_channels_ready
            else "assembly_in_progress"
            if blockers
            else "closure_candidate"
        ),
        "blockers": blockers,
        "current_frontier": blockers[:2],
        "smallest_constructive_missing_object": effective_next_missing_object,
        "summary": {
            "full_unquenched_status": full_unquenched.get("status"),
            "population_status": population_status,
            "next_missing_object": effective_next_missing_object,
            "sequence_population_status": sequence_status,
            "cfg_source_payload_status": payload_status,
            "sequence_evaluation_status": evaluation_status,
            "lambda_msbar3_status": (
                (full_unquenched.get("predictive_input_status") or {}).get("Lambda_MSbar_3_gev")
            ),
            "contraction_plan_status": contraction_status,
            "stable_theorem_after_population": full_unquenched.get("next_theorem_after_population"),
            "surrogate_execution_bridge_status": (
                surrogate_bridge.get("status") if surrogate_bridge is not None else None
            ),
            "production_geometry_summary_artifact": (
                geometry_summary.get("artifact") if geometry_summary is not None else None
            ),
            "production_closure_grade": (
                closure_report.get("closure_grade") if closure_report is not None else None
            ),
        },
    }


def build_audit(
    full_unquenched: dict[str, Any],
    contraction_plan: dict[str, Any] | None,
    sequence_population: dict[str, Any] | None,
    cfg_source_payload: dict[str, Any] | None,
    sequence_evaluation: dict[str, Any] | None,
    surrogate_bridge: dict[str, Any] | None,
    geometry_summary: dict[str, Any] | None,
    closure_report: dict[str, Any] | None,
    references: dict[str, Any],
) -> dict[str, Any]:
    classification = classify_pipeline(
        full_unquenched,
        contraction_plan,
        sequence_population,
        cfg_source_payload,
        sequence_evaluation,
        surrogate_bridge,
        geometry_summary,
        closure_report,
    )
    seeded = (full_unquenched.get("population_contract") or {}).get("status") == "predictive_ensemble_seeded_candidate"
    stable_channels_ready = bool((closure_report or {}).get("public_unsuppression_ready"))
    sequence_status = None if sequence_population is None else sequence_population.get("status")
    payload_status = None if cfg_source_payload is None else cfg_source_payload.get("status")
    evaluation_status = None if sequence_evaluation is None else sequence_evaluation.get("status")
    measure_realization_open = (
        evaluation_status == "awaiting_measure_evaluation"
        or payload_status == "law_closed_waiting_measure_realization"
        or sequence_status == "law_closed_waiting_measure_evaluation"
    )
    schedule_emitted = bool((cfg_source_payload or {}).get("support_realization_schedule"))
    return {
        "artifact": "oph_current_hadron_lane_audit",
        "generated_utc": _timestamp(),
        "current_frontier": list(classification.get("current_frontier", [])),
        "smallest_constructive_missing_object": classification.get("smallest_constructive_missing_object"),
        "public_status_surface_policy": {
            "include_hadrons_by_default": stable_channels_ready,
            "reason": (
                "stable_channels_ready_rho_pending"
                if stable_channels_ready
                else "hadron_pipeline_not_closed"
            ),
        },
        "reference_targets_gev": {
            "proton": float(references["proton"]["value_gev"]),
            "neutron": float(references["neutron"]["value_gev"]),
            "neutral_pion": float(references["neutral_pion"]["value_gev"]),
            "rho_770_0": float(references["rho_770_0"]["value_gev"]),
        },
        "pipeline_classification": classification,
        "promotion_verdict": (
            "stable_channels_ready_rho_pending"
            if stable_channels_ready
            else "suppress_from_public_surface"
        ),
        "surrogate_execution_bridge": (
            None
            if surrogate_bridge is None
            else {
                "artifact": surrogate_bridge.get("artifact"),
                "status": surrogate_bridge.get("status"),
                "public_promotion_allowed": surrogate_bridge.get("public_promotion_allowed"),
                "surrogate_execution": surrogate_bridge.get("surrogate_execution"),
                "surrogate_finest_ensemble": surrogate_bridge.get("surrogate_finest_ensemble"),
                "surrogate_error_summary": surrogate_bridge.get("surrogate_error_summary"),
            }
        ),
        "production_geometry_summary": (
            None
            if geometry_summary is None
            else {
                "artifact": geometry_summary.get("artifact"),
                "totals": geometry_summary.get("totals"),
            }
        ),
        "strict_missing_program": [
            "nonperturbative_hadron_closure_law",
            classification.get("smallest_constructive_missing_object"),
            "full_baryon_contractions" if not seeded else "StableChannelForwardWindowConvergence",
            "rho_resonance_extraction",
        ],
        "recommended_next_predictive_artifact": {
            "name": (
                "rho_phase_shift_fit.json"
                if stable_channels_ready
                else
                "backend_correlator_dump.production.json"
                if measure_realization_open
                else ("oph_hadron_stable_channel_sequence_population" if seeded else "oph_hadron_full_unquenched_correlator")
            ),
            "scope": (
                ["rho_scattering"]
                if stable_channels_ready
                else ["pi_iso", "N_iso"] if seeded else ["pi_iso", "N_iso", "rho_scattering"]
            ),
            "why": (
                "The stable channels are promotable on the executed seeded family, so the remaining hadron work is the separate rho finite-volume scattering readout."
                if stable_channels_ready
                else
                "The sequence-emission, cfg/source jackknife law, runtime receipt, and frozen execution schema are all explicit, so the next productive move is the real backend correlator dump on the seeded 2+1 family before forward-window convergence and production systematics."
                if measure_realization_open and schedule_emitted
                else
                "The sequence-emission and cfg/source jackknife laws are fixed on the seeded ensemble family, so the next productive move is realizing the per-ensemble cfg/source arrays on that family before the forward-window convergence theorem."
                if measure_realization_open
                else
                "The ensemble family is seeded, so the next productive move is a dedicated per-ensemble sequence-population artifact "
                "for stable channels before the forward-window convergence theorem."
                if seeded
                else
                "The stable-channel and rho readout surfaces already exist. The remaining work is "
                "to populate the upstream unquenched correlator producer and its ensemble seed."
            ),
        },
        "minimal_closure_frontier": [
            "full_unquenched_correlator",
            "stable_channel_sequence_population",
            "stable_channel_sequence_evaluation",
            "stable_channel_groundstate_readout",
            "finite_volume_resonance_and_spectrum_readout",
        ],
        "smallest_missing_theorem_after_full_unquenched": "StableChannelForwardWindowConvergence",
        "latest_external_hints": [
            {
                "source": "e131_hadron_binding_emitter_push.response.md",
                "summary": (
                    "The minimal hadron closure object is a stable-channel ground-state "
                    "readout for pi_iso and N_iso built from fixed color-singlet Euclidean "
                    "two-point functions, with rho deferred to a later scattering readout."
                ),
            },
            {
                "source": "e132_hadron_resonance_unquenching_push.response.md",
                "summary": (
                    "The smallest hadron closure move is the jump to a full unquenched correlator together with "
                    "finite-volume resonance and spectrum readout."
                ),
            },
            {
                "source": "e145_hadron_exactness_bundle_push.response.md",
                "summary": (
                    "The next missing object is an explicit full unquenched correlator "
                    "feeding pi_iso, N_iso direct-minus-exchange, and rho-channel correlation matrices."
                ),
            },
            {
                "source": "e152_hadron_stable_convergence_push.response.md",
                "summary": (
                    "After the full unquenched correlator exists, the stable-channel theorem payload "
                    "is convergence of the long-time log-ratio readouts for pi_iso and N_iso, "
                    "which then promotes am_ground, ratio_to_lambda_msbar3, and mass_gev."
                ),
            },
            {
                "source": "e154_hadron_unquenched_producer_push.response.md",
                "summary": (
                    "The producer must run in the upstream direction: declare a predictive-only 2+1-flavor "
                    "QCD measure pushforward on OPH quark/QCD descendants, then let stable-channel and rho "
                    "readouts consume that correlator family."
                ),
            },
            {
                "source": "e158_hadron_unquenched_ensemble_population_push.response.md",
                "summary": (
                    "The next producer-side bridge is a Lambda-based 2+1 ensemble-population law: "
                    "consume Lambda_MSbar^(3) plus OPH quark descendants, form rho_l/rho_s, and refine "
                    "a predictive unquenched ensemble family in Lambda-units before stable-channel convergence."
                ),
            },
            {
                "source": "e170_hadron_sequence_emission_or_obstruction_push.response.md",
                "summary": (
                    "The stable-channel sequence-emission law is already fixed by the seeded 2+1 measure, the stable-channel operator rule, "
                    "and the closed nucleon direct-minus-exchange contraction; the remaining local step is evaluation of per-ensemble sequences before convergence."
                ),
            },
            {
                "source": "e174_hadron_groundstate_convergence_push.response.md",
                "summary": (
                    "The correct finite-T theorem is StableChannelForwardWindowConvergence; the real next local object is the sequence evaluator that emits corr_t, am_eff_t, and forward-window certificates."
                ),
            },
        ],
        "promotion_gate": "StableChannelForwardWindowConvergence",
        "notes": [
            "Exact quark masses are not sufficient by themselves because hadron masses depend dominantly on nonperturbative QCD binding energy.",
            "The current hadron pipeline is anchored on the unquenched correlator producer, the stable-channel readout, and the rho scattering readout.",
            (
                "The seeded stable-channel branch now has executed cfg/source arrays, forward-window certificates, and published statistical/systematic fields, so `pi_iso` and `N_iso` are promotable while `rho` remains a separate scattering task."
                if stable_channels_ready
                else
                "The ensemble seed is now fixed from Lambda_MSbar^(3) and the quark descendants; the sequence-emission and jackknife evaluation laws are fixed too, so the next local move is realizing the per-ensemble stable-channel cfg/source payload arrays on that family."
                if measure_realization_open and not schedule_emitted
                else "The ensemble seed is now fixed from Lambda_MSbar^(3) and the quark descendants; the sequence-emission and jackknife evaluation laws are fixed too, so the next local move is executing the fixed RHMC/HMC support schedule and writing the resulting per-ensemble stable-channel cfg/source payload arrays."
                if measure_realization_open
                else
                "The ensemble seed is now fixed from Lambda_MSbar^(3) and the quark descendants; the next producer-side bridge is a dedicated stable-channel sequence-population artifact on that seeded family."
                if seeded
                else "The remaining producer-side bridge is ensemble population from Lambda_MSbar^(3) and the quark descendants."
            ),
            (
                "Stable-channel public unsuppression is now driven by the closure validator rather than another theorem-only artifact."
                if stable_channels_ready
                else
                "The runtime receipt `(N_therm, N_sep)` is the emitted execution-and-systematics contract. With that contract explicit, the first honest live residual object is the production backend correlator dump on the seeded 2+1 family."
                if measure_realization_open
                else "The stable-channel cfg/source arrays are not yet the live blocker on this branch."
            ),
            (
                "A separate surrogate HMC/RHMC execution bridge is now recorded too: it closes the receipt/writeback/evaluation/convergence/systematics path on the emitted schema, but it is diagnostic only and does not replace production unquenched execution."
                if surrogate_bridge is not None
                else "No surrogate execution-bridge diagnostic is attached to this audit."
            ),
            "The runtime receipt should now be treated as an execution-and-systematics contract: once it is filled, the honest remaining work is executed unquenching, array writeback, and published continuum/volume/chiral/statistical budgets.",
            "After the evaluator is populated, the next stable-channel theorem is forward-window convergence of the ground-state extraction.",
            "The public table should stay suppressed until that pipeline emits stable-channel masses on its own path.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit the current local hadron pipeline.")
    parser.add_argument("--full-unquenched", default=str(DEFAULT_FULL_UNQUENCHED))
    parser.add_argument("--contraction-plan", default=str(DEFAULT_CONTRACTION_PLAN))
    parser.add_argument("--sequence-population", default=str(DEFAULT_SEQUENCE_POPULATION))
    parser.add_argument("--cfg-source-payload", default=str(DEFAULT_CFG_SOURCE_PAYLOAD))
    parser.add_argument("--sequence-evaluation", default=str(DEFAULT_SEQUENCE_EVALUATION))
    parser.add_argument("--surrogate-bridge", default=str(DEFAULT_SURROGATE_BRIDGE))
    parser.add_argument("--geometry-summary", default=str(DEFAULT_GEOMETRY_SUMMARY))
    parser.add_argument("--closure-report", default=str(DEFAULT_CLOSURE_REPORT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    full_unquenched = json.loads(pathlib.Path(args.full_unquenched).read_text(encoding="utf-8"))
    contraction_plan_path = pathlib.Path(args.contraction_plan)
    contraction_plan = (
        json.loads(contraction_plan_path.read_text(encoding="utf-8"))
        if contraction_plan_path.exists()
        else None
    )
    sequence_population_path = pathlib.Path(args.sequence_population)
    sequence_evaluation_path = pathlib.Path(args.sequence_evaluation)
    sequence_population = (
        json.loads(sequence_population_path.read_text(encoding="utf-8"))
        if sequence_population_path.exists()
        else None
    )
    cfg_source_payload_path = pathlib.Path(args.cfg_source_payload)
    cfg_source_payload = (
        json.loads(cfg_source_payload_path.read_text(encoding="utf-8"))
        if cfg_source_payload_path.exists()
        else None
    )
    sequence_evaluation = (
        json.loads(sequence_evaluation_path.read_text(encoding="utf-8"))
        if sequence_evaluation_path.exists()
        else None
    )
    surrogate_bridge_path = pathlib.Path(args.surrogate_bridge)
    surrogate_bridge = (
        json.loads(surrogate_bridge_path.read_text(encoding="utf-8"))
        if surrogate_bridge_path.exists()
        else None
    )
    geometry_summary_path = pathlib.Path(args.geometry_summary)
    geometry_summary = (
        json.loads(geometry_summary_path.read_text(encoding="utf-8"))
        if geometry_summary_path.exists()
        else None
    )
    closure_report_path = pathlib.Path(args.closure_report)
    closure_report = (
        json.loads(closure_report_path.read_text(encoding="utf-8"))
        if closure_report_path.exists()
        else None
    )
    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))["entries"]
    audit = build_audit(
        full_unquenched,
        contraction_plan,
        sequence_population,
        cfg_source_payload,
        sequence_evaluation,
        surrogate_bridge,
        geometry_summary,
        closure_report,
        references,
    )

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
