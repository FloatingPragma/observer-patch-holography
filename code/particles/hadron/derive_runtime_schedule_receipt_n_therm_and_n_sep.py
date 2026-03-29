#!/usr/bin/env python3
"""Emit the runtime receipt shell for hadron fixed-schedule execution.

Chain role: make the externally supplied `(N_therm, N_sep)` receipt explicit so
the hadron lane stops searching for a theorem-only selector on that surface.

Mathematics: deterministic schedule bookkeeping only; no hadron masses are
derived here.

OPH-derived inputs: the seeded cfg/source payload shell and its fixed RHMC/HMC
execution contract.

Output: the non-theorem runtime receipt artifact consumed by payload and
evaluation writeback.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.hadron.production_execution_support import fill_runtime_receipt

DEFAULT_PAYLOAD = ROOT / "particles" / "runs" / "hadron" / "stable_channel_cfg_source_measure_payload.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "runtime_schedule_receipt_N_therm_and_N_sep.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(payload: dict) -> dict:
    schedule = dict(payload.get("support_realization_schedule", {}))
    required = dict(schedule.get("required_schedule_scalars", {}))
    ensemble_schedule = []
    for item in schedule.get("ensemble_schedule", []):
        ensemble_schedule.append(
            {
                "ensemble_id": item["ensemble_id"],
                "cfg_ids": list(item.get("cfg_ids", [])),
                "trajectory_stop_by_cfg": dict(item.get("trajectory_stop_by_cfg", {})),
                "trajectory_stop_by_cfg_formula": dict(item.get("trajectory_stop_by_cfg_formula", {})),
                "n_src_per_cfg": item.get("n_src_per_cfg"),
                "t_extent": item.get("t_extent", item.get("T")),
            }
        )
    return {
        "artifact": "runtime_schedule_receipt_N_therm_and_N_sep",
        "generated_utc": _timestamp(),
        "status": "awaiting_external_runtime_inputs",
        "receipt_kind": "non_theorem_external_runtime_receipt",
        "completion_mode": "execution_and_systematics_contract",
        "execution_artifact": "oph_hadron_fixed_schedule_rhmc_hmc_execution",
        "source_artifact": payload.get("artifact"),
        "required_schedule_scalars": {
            "N_therm": required.get("N_therm"),
            "N_sep": required.get("N_sep"),
        },
        "execution_and_systematics_contract": {
            "runtime_inputs_are_external": True,
            "execution_receipt_required_before_public_hadron_rows": True,
            "published_systematics_required_after_execution": True,
            "statistical_contract": {
                "sigma_stat_formula": "JKstderr(am_ground_X)",
                "effective_cfg_count_formula": "n_eff_cfg = n_cfg / (2 * tau_int_cfg)",
            },
            "systematics_contract": {
                "sigma_sys_formula": "sqrt(delta_cont^2 + delta_vol^2 + delta_chi^2)",
                "required_components": [
                    "delta_cont",
                    "delta_vol",
                    "delta_chi",
                ],
            },
            "provenance_contract": {
                "required_execution_fields": [
                    "required_schedule_scalars.N_therm",
                    "required_schedule_scalars.N_sep",
                    "execution_contract.ensemble_schedule[*].trajectory_stop_by_cfg",
                    "execution_contract.writeback_targets",
                ],
                "required_publication_fields": [
                    "stable_channel_sequence_evaluation.ensemble_evaluations[*].pi_iso.published_systematics",
                    "stable_channel_sequence_evaluation.ensemble_evaluations[*].pi_iso.published_statistical_error",
                    "stable_channel_sequence_evaluation.ensemble_evaluations[*].N_iso.published_systematics",
                    "stable_channel_sequence_evaluation.ensemble_evaluations[*].N_iso.published_statistical_error",
                ],
            },
        },
        "execution_contract": {
            "kernel_id": schedule.get("kernel_id"),
            "initial_configuration": schedule.get("initial_configuration"),
            "stop_time_formula": "t_stop(n,cfg_index) = N_therm + cfg_index*N_sep",
            "cfg_realization_formula": schedule.get("cfg_realization_formula"),
            "seed_bytes_formula": schedule.get("seed_bytes_formula"),
            "seed_lookup_path": schedule.get("seed_lookup_path"),
            "source_lookup_path": schedule.get("source_lookup_path"),
            "source_set_formula": schedule.get("source_set_formula"),
            "ensemble_schedule": ensemble_schedule,
            "emission_formulas": dict(schedule.get("emission_formulas", {})),
            "writeback_targets": dict(schedule.get("writeback_targets", {})),
            "next_single_residual_object_after_execution": (
                schedule.get("conditional_execution_receipt", {}) or {}
            ).get("next_single_residual_object_after_execution", "oph_hadron_stable_channel_sequence_evaluator"),
        },
        "notes": [
            "This receipt is the final honest bridge between theorem-side hadron preparation and executed runtime/systematics work.",
            "It does not emit N_therm or N_sep itself; those remain external execution inputs on the present hadron surface.",
            "Once the receipt is filled and execution is performed, the next work is measure writeback, published statistical errors, and declared continuum/volume/chiral systematics rather than another theorem-only selector search.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the hadron runtime schedule receipt artifact.")
    parser.add_argument("--payload", default=str(DEFAULT_PAYLOAD))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--n-therm", type=int, default=None)
    parser.add_argument("--n-sep", type=int, default=None)
    parser.add_argument(
        "--schedule-provenance",
        default=None,
        help="Explicit provenance string for externally supplied schedule scalars.",
    )
    args = parser.parse_args()

    payload = json.loads(Path(args.payload).read_text(encoding="utf-8"))
    artifact = build_artifact(payload)
    artifact = fill_runtime_receipt(
        artifact,
        n_therm=args.n_therm,
        n_sep=args.n_sep,
        schedule_provenance=args.schedule_provenance,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
