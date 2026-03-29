#!/usr/bin/env python3
"""Realize the cfg/source support shell for the active hadron ensembles.

Chain role: instantiate deterministic config/source identifiers for the seeded
stable-channel ensembles before any correlator arrays are populated.

Mathematics: canonical JSON seed material plus SHA-256 hashing for reproducible
cfg/source labels on each ensemble.

OPH-derived inputs: the seeded ensemble family, stable-channel population
schema, and baryon contraction plan.

Output: the cfg/source payload artifact consumed by the stable-channel sequence
evaluator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.hadron.production_execution_support import ingest_dump_into_payload

DEFAULT_FULL_UNQUENCHED = ROOT / "particles" / "runs" / "hadron" / "full_unquenched_correlator.json"
DEFAULT_SEQUENCE_POPULATION = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_population.json"
DEFAULT_CONTRACTION_PLAN = ROOT / "particles" / "runs" / "hadron" / "proton_contraction_plan.json"
DEFAULT_RUNTIME_RECEIPT = ROOT / "particles" / "runs" / "hadron" / "runtime_schedule_receipt_N_therm_and_N_sep.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "stable_channel_cfg_source_measure_payload.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _cfg_seed_material(
    ensemble_id: str,
    beta: float,
    l_size: int,
    t_size: int,
    am_l: float,
    am_s: float,
    cfg_index: int,
) -> list[object]:
    return [
        ensemble_id,
        format(beta, ".17g"),
        l_size,
        t_size,
        format(am_l, ".17g"),
        format(am_s, ".17g"),
        cfg_index,
    ]


def _cfg_seed_hash(material: list[object]) -> str:
    serialized = json.dumps(material, ensure_ascii=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def build_artifact(
    full_unquenched: dict,
    sequence_population: dict,
    contraction_plan: dict,
    runtime_receipt: dict | None = None,
) -> dict:
    runtime_receipt = runtime_receipt or {}
    receipt_scalars = dict(runtime_receipt.get("required_schedule_scalars", {}))
    receipt_filled = (
        receipt_scalars.get("N_therm") is not None and receipt_scalars.get("N_sep") is not None
    )
    receipt_schedule = {
        str(entry["ensemble_id"]): entry
        for entry in (runtime_receipt.get("execution_contract") or {}).get("ensemble_schedule", [])
    }
    ensemble_payloads = []
    ensemble_schedule = []
    for ensemble in sequence_population.get("ensemble_sequences", []):
        ensemble_id = ensemble["ensemble_id"]
        l_size = int(ensemble["L"])
        t_size = int(ensemble["T"])
        cfg_ids = [f"{ensemble_id}__cfg0", f"{ensemble_id}__cfg1"]
        cfg_seed_materials = {
            cfg_id: _cfg_seed_material(
                ensemble_id,
                float(ensemble["beta"]),
                l_size,
                t_size,
                float(ensemble["am_l"]),
                float(ensemble["am_s"]),
                cfg_index,
            )
            for cfg_index, cfg_id in enumerate(cfg_ids)
        }
        cfg_seed_hashes = {
            cfg_id: _cfg_seed_hash(material)
            for cfg_id, material in cfg_seed_materials.items()
        }
        source_descriptors = [
            {"src_id": "s0", "kind": "point", "coords": [0, 0, 0, 0]},
            {"src_id": "s1", "kind": "point", "coords": [l_size // 2, l_size // 2, l_size // 2, t_size // 2]},
        ]
        ensemble_payloads.append(
            {
                "ensemble_id": ensemble_id,
                "family_index": ensemble["family_index"],
                "beta": ensemble["beta"],
                "L": l_size,
                "T": t_size,
                "aLambda_msbar3": ensemble["aLambda_msbar3"],
                "am_l": ensemble["am_l"],
                "am_s": ensemble["am_s"],
                "cfg_ids": cfg_ids,
                "n_cfg": len(cfg_ids),
                "source_descriptors_by_cfg": {
                    cfg_id: [dict(source) for source in source_descriptors]
                    for cfg_id in cfg_ids
                },
                "n_src_per_cfg": len(source_descriptors),
                "t_support": ensemble["t_support"],
                "t_lambda_msbar3": ensemble["t_lambda_msbar3"],
                "cfg_seed_status": "deterministic_sha256_seed_bridge_closed_cfg_arrays_unrealized",
                "cfg_realization_contract": {
                    "cfg_id_formula": f"{ensemble_id}__cfg{{j}} for j = 0,1",
                    "cfg_seed_hash_formula": "SHA256(JSON([ensemble_id, \"%.17g\" % beta, L, T, \"%.17g\" % am_l, \"%.17g\" % am_s, cfg_index]))",
                    "cfg_seed_hash_algorithm": "sha256",
                    "cfg_seed_hash_serialization": "json_array_with_fixed_17_digit_float_strings",
                    "cfg_seed_hash_inputs": ["ensemble_id", "beta", "L", "T", "am_l", "am_s", "cfg_index"],
                    "cfg_seed_materials": cfg_seed_materials,
                    "cfg_seed_hashes": cfg_seed_hashes,
                    "cfg_realization_kernel": "fixed_schedule_rhmc_hmc",
                    "cfg_realization_formula": "U_{n,c} = K_n^(N_therm + c*N_sep)(U_cold; seed_{n,c})",
                    "measure_invariant": "dmu_n(U)",
                    "initial_configuration": "U_cold = 1",
                    "source_set_formula": "[[0,0,0,0], [L//2, L//2, L//2, T//2]]",
                    "irreducible_stochastic_primitive": "realized_cfgs_on_seeded_2p1_family",
                },
                "cfg_source_shape_contract": {
                    "n_cfg": len(cfg_ids),
                    "n_src_per_cfg": len(source_descriptors),
                    "t_extent": len(ensemble["t_support"]),
                },
                "pi_iso": {
                    "cfg_source_corr_t": [],
                    "cfg_source_corr_shape": [len(cfg_ids), len(source_descriptors), len(ensemble["t_support"])],
                    "cfg_source_corr_formula": "sum_x Re tr_c,spin[gamma5 S_l(x;src) gamma5 S_l(src;x)]",
                },
                "N_iso": {
                    "cfg_source_corr_direct_t": [],
                    "cfg_source_corr_exchange_t": [],
                    "cfg_source_corr_t": [],
                    "cfg_source_corr_shape": [len(cfg_ids), len(source_descriptors), len(ensemble["t_support"])],
                    "cfg_source_corr_direct_formula": "sum_x G_d, with G_d = Sc . tr(Gamma Sa Gamma_B Sb^T)",
                    "cfg_source_corr_exchange_formula": "sum_x G_x, with G_x = Sca Gamma_B^T Sb^T Gamma^T Sac",
                    "cfg_source_corr_formula": "cfg_source_corr_direct_t - cfg_source_corr_exchange_t",
                },
            }
        )
        ensemble_schedule.append(
            {
                "ensemble_id": ensemble_id,
                "family_index": ensemble["family_index"],
                "beta": ensemble["beta"],
                "L": l_size,
                "T": t_size,
                "aLambda_msbar3": ensemble["aLambda_msbar3"],
                "am_l": ensemble["am_l"],
                "am_s": ensemble["am_s"],
                "cfg_ids": cfg_ids,
                "n_cfg": len(cfg_ids),
                "n_src_per_cfg": len(source_descriptors),
                "trajectory_stop_formula": "N_therm + cfg_index*N_sep",
                "trajectory_stop_by_cfg": dict(
                    (receipt_schedule.get(ensemble_id) or {}).get("trajectory_stop_by_cfg", {})
                ),
                "trajectory_stop_by_cfg_formula": dict(
                    (receipt_schedule.get(ensemble_id) or {}).get("trajectory_stop_by_cfg_formula", {})
                ),
                "pi_iso_cfg_source_corr_shape": [len(cfg_ids), len(source_descriptors), len(ensemble["t_support"])],
                "N_iso_cfg_source_corr_shape": [len(cfg_ids), len(source_descriptors), len(ensemble["t_support"])],
            }
        )
    return {
        "artifact": "oph_hadron_stable_channel_cfg_source_measure_payload",
        "generated_utc": _timestamp(),
        "payload_law_id": "oph_qcd_2p1_stable_channel_cfg_source_measure_payload",
        "status": "law_closed_waiting_measure_realization",
        "proof_status": "law_closed_waiting_measure_realization",
        "predictive_promotion_allowed": False,
        "payload_realization_status": (
            "receipt_filled_backend_arrays_missing"
            if receipt_filled
            else "external_schedule_receipt_required_arrays_missing"
        ),
        "realized_support_primitive": "runtime_schedule_receipt_N_therm_and_N_sep",
        "cfg_source_payload_kind": "source_resolved_contraction_map",
        "cfg_support_realization_status": (
            "deterministic_sha256_seed_bridge_closed_receipt_filled_arrays_missing"
            if receipt_filled
            else "deterministic_sha256_seed_bridge_closed_external_schedule_receipt_missing"
        ),
        "support_realization_contract": {
            "n_cfg_per_ensemble": 2,
            "n_src_per_cfg": 2,
            "cfg_seed_hash_formula": "SHA256(JSON([ensemble_id, \"%.17g\" % beta, L, T, \"%.17g\" % am_l, \"%.17g\" % am_s, cfg_index]))",
            "cfg_seed_hash_algorithm": "sha256",
            "cfg_seed_hash_serialization": "json_array_with_fixed_17_digit_float_strings",
            "cfg_realization_kernel": "fixed_schedule_rhmc_hmc",
            "cfg_realization_formula": "U_{n,c} = K_n^(N_therm + c*N_sep)(U_cold; seed_{n,c})",
            "measure_invariant": "dmu_n(U)",
            "source_set_formula": "[[0,0,0,0], [L//2, L//2, L//2, T//2]]",
            "irreducible_stochastic_primitive": "realized_cfgs_on_seeded_2p1_family",
        },
        "support_realization_schedule": {
            "artifact": "oph_hadron_fixed_schedule_rhmc_hmc_execution",
            "status": (
                "receipt_filled_waiting_backend_execution"
                if receipt_filled
                else "external_receipt_required_before_execution"
            ),
            "runtime_receipt_artifact": runtime_receipt.get("artifact"),
            "kernel_id": "fixed_schedule_rhmc_hmc",
            "source_artifact": "oph_hadron_stable_channel_cfg_source_measure_payload",
            "realized_support_primitive": "executed_fixed_schedule_rhmc_hmc_on_seeded_2p1_family",
            "cfg_realization_formula": "U_{n,c} = K_n^(N_therm + cfg_index*N_sep)(U_cold; seed_{n,c})",
            "measure_invariant": "dmu_n(U)",
            "initial_configuration": "U_cold = 1",
            "seed_lookup_path": "ensemble_payloads[*].cfg_realization_contract.cfg_seed_hashes[cfg_id]",
            "seed_bytes_formula": "seed_{n,c} = bytes.fromhex(cfg_seed_hash_{n,c})",
            "source_lookup_path": "ensemble_payloads[*].source_descriptors_by_cfg[cfg_id]",
            "source_set_formula": "[[0,0,0,0], [L//2, L//2, L//2, T//2]]",
            "required_schedule_scalars": {
                "N_therm": receipt_scalars.get("N_therm"),
                "N_sep": receipt_scalars.get("N_sep"),
            },
            "strict_externality_theorem": "no_predictive_selector_for_schedule_scalars_on_current_surface",
            "runtime_obstruction_status": (
                "external_runtime_receipt_filled_backend_arrays_pending"
                if receipt_filled
                else "schedule_scalars_are_external_runtime_inputs"
            ),
            "runtime_obstruction_reason": "the deterministic seed law and writeback map are fixed, but N_therm and N_sep are irreducibly external on the present emitted surface",
            "strongest_strictly_smaller_constructive_primitive": "runtime_schedule_receipt_N_therm_and_N_sep",
            "emission_formulas": {
                "pi_iso.cfg_source_corr_t": "p_pi^(n,c,s)(t) = sum_x Re tr_c,spin[gamma5 S_l(x;s) gamma5 S_l(s;x)]",
                "N_iso.cfg_source_corr_direct_t": "p_N,dir^(n,c,s)(t) = sum_x G_d",
                "N_iso.cfg_source_corr_exchange_t": "p_N,ex^(n,c,s)(t) = sum_x G_x",
                "N_iso.cfg_source_corr_t": "p_N^(n,c,s)(t) = p_N,dir^(n,c,s)(t) - p_N,ex^(n,c,s)(t)",
            },
            "ensemble_schedule": ensemble_schedule,
            "conditional_execution_receipt": {
                "stop_time_formula": "t_stop(n,cfg_index) = N_therm + cfg_index*N_sep",
                "cfg_source_array_write_shapes": {
                    "pi_iso.cfg_source_corr_t": "shape [n_cfg, n_src_per_cfg, t_extent]",
                    "N_iso.cfg_source_corr_direct_t": "shape [n_cfg, n_src_per_cfg, t_extent]",
                    "N_iso.cfg_source_corr_exchange_t": "shape [n_cfg, n_src_per_cfg, t_extent]",
                    "N_iso.cfg_source_corr_t": "shape [n_cfg, n_src_per_cfg, t_extent]",
                },
                "next_single_residual_object_after_execution": "oph_hadron_stable_channel_sequence_evaluator",
            },
            "writeback_targets": {
                "payload_artifact": "particles/runs/hadron/stable_channel_cfg_source_measure_payload.json",
                "payload_fields": [
                    "ensemble_payloads[*].pi_iso.cfg_source_corr_t",
                    "ensemble_payloads[*].N_iso.cfg_source_corr_direct_t",
                    "ensemble_payloads[*].N_iso.cfg_source_corr_exchange_t",
                    "ensemble_payloads[*].N_iso.cfg_source_corr_t",
                ],
                "evaluation_artifact": "particles/runs/hadron/stable_channel_sequence_evaluation.json",
                "evaluation_fields": [
                    "ensemble_evaluations[*].pi_iso.cfg_source_corr_t",
                    "ensemble_evaluations[*].N_iso.cfg_source_corr_direct_t",
                    "ensemble_evaluations[*].N_iso.cfg_source_corr_exchange_t",
                    "ensemble_evaluations[*].N_iso.cfg_source_corr_t",
                ],
            },
        },
        "strict_upstream_blocker": None,
        "source_artifacts": {
            "full_unquenched_correlator": full_unquenched.get("artifact"),
            "stable_channel_sequence_population": sequence_population.get("artifact"),
            "proton_contraction_plan": contraction_plan.get("artifact"),
        },
        "realization_formulas": {
            "cfg_realization": "U_{n,c} ~ dmu_n(U) with deterministic cfg seed derived from ensemble data and cfg index",
            "source_realization": "S_{n,c} = {[0,0,0,0], [L_n//2, L_n//2, L_n//2, T_n//2]}",
            "pi_cfg_source_corr": "p_pi^(n,c,s)(t) = sum_x Re tr_c,spin[gamma5 S_l(x;s) gamma5 S_l(s;x)]",
            "N_cfg_source_corr_direct": "p_N,dir^(n,c,s)(t) = sum_x G_d",
            "N_cfg_source_corr_exchange": "p_N,ex^(n,c,s)(t) = sum_x G_x",
            "N_cfg_source_corr": "p_N^(n,c,s)(t) = p_N,dir^(n,c,s)(t) - p_N,ex^(n,c,s)(t)",
            "cfg_average_pi": "bar_p_pi^(n,c)(t) = (1/|S_{n,c}|) * sum_{s in S_{n,c}} p_pi^(n)(c,s;t)",
            "cfg_average_N": "bar_p_N^(n,c)(t) = (1/|S_{n,c}|) * sum_{s in S_{n,c}} p_N^(n)(c,s;t)",
            "cfg_ensemble_average": "C_X^(n)(t) = (1/|C_n|) * sum_{c in C_n} bar_p_X^(n,c)(t)",
        },
        "ensemble_payloads": ensemble_payloads,
        "smallest_constructive_missing_object": (
            "backend_correlator_dump.production.json from real production RHMC/HMC execution on the theorem-emitted seeded family"
            if receipt_filled
            else "runtime_schedule_receipt_N_therm_and_N_sep"
        ),
        "next_theorem_after_measure_realization": "StableChannelForwardWindowConvergence",
        "notes": [
            "This artifact fixes the cfg/source payload law for the stable-channel measure arrays on each seeded ensemble.",
            "The deterministic point-source support contract is now populated on each seeded ensemble.",
            "The cfg seed bridge is fixed by canonical JSON serialization plus SHA-256, and the emitted next runtime object is the external schedule receipt `(N_therm, N_sep)` before execution of the fixed RHMC/HMC schedule on the seeded 2+1 family.",
            "After realized cfg arrays exist, the next hadron theorem remains forward-window convergence.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the stable-channel cfg/source measure payload artifact.")
    parser.add_argument("--full-unquenched", default=str(DEFAULT_FULL_UNQUENCHED))
    parser.add_argument("--sequence-population", default=str(DEFAULT_SEQUENCE_POPULATION))
    parser.add_argument("--contraction-plan", default=str(DEFAULT_CONTRACTION_PLAN))
    parser.add_argument("--runtime-receipt", default=str(DEFAULT_RUNTIME_RECEIPT))
    parser.add_argument(
        "--backend-dump",
        default=None,
        help="Optional validated production dump to write back into cfg/source arrays.",
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    full_unquenched = json.loads(Path(args.full_unquenched).read_text(encoding="utf-8"))
    sequence_population = json.loads(Path(args.sequence_population).read_text(encoding="utf-8"))
    contraction_plan = json.loads(Path(args.contraction_plan).read_text(encoding="utf-8"))
    runtime_receipt_path = Path(args.runtime_receipt)
    runtime_receipt = (
        json.loads(runtime_receipt_path.read_text(encoding="utf-8"))
        if runtime_receipt_path.exists()
        else None
    )
    artifact = build_artifact(full_unquenched, sequence_population, contraction_plan, runtime_receipt)
    if args.backend_dump:
        backend_dump = json.loads(Path(args.backend_dump).read_text(encoding="utf-8"))
        artifact = ingest_dump_into_payload(artifact, backend_dump, runtime_receipt)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
