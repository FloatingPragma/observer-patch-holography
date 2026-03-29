#!/usr/bin/env python3
"""Attach evaluation and convergence scaffolding to stable-channel sequences.

Chain role: combine the sequence shell with cfg/source payload metadata so the
hadron lane can accumulate jackknife averages and forward-window certificates.

Mathematics: cfg/source averaging contracts, jackknife axes, and placeholder
forward-window/convergence slots.

OPH-derived inputs: the stable-channel sequence population and deterministic
cfg/source payload shell.

Output: the stable-channel evaluation artifact that becomes numerical once
cfg/source arrays are realized.
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

from particles.hadron.production_execution_support import populate_evaluation_from_dump

DEFAULT_SEQUENCE_POPULATION = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_population.json"
DEFAULT_CFG_SOURCE_PAYLOAD = ROOT / "particles" / "runs" / "hadron" / "stable_channel_cfg_source_measure_payload.json"
DEFAULT_RUNTIME_RECEIPT = ROOT / "particles" / "runs" / "hadron" / "runtime_schedule_receipt_N_therm_and_N_sep.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_evaluation.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    sequence_population: dict,
    cfg_source_payload: dict | None = None,
    runtime_receipt: dict | None = None,
) -> dict:
    cfg_source_payload = cfg_source_payload or {}
    runtime_receipt = runtime_receipt or {}
    payload_index = {
        ensemble["ensemble_id"]: ensemble
        for ensemble in cfg_source_payload.get("ensemble_payloads", [])
    }
    ensemble_evaluations = []
    for ensemble in sequence_population.get("ensemble_sequences", []):
        payload_ensemble = payload_index.get(ensemble["ensemble_id"], {})
        payload_artifact = cfg_source_payload.get("artifact", sequence_population.get("artifact"))
        ensemble_ref = f"{payload_artifact}::{ensemble['ensemble_id']}"
        t_size = len(ensemble["t_support"])
        raw_forward_window_cardinality = max(int(ensemble["T"]) // 2 - 2, 0)
        ensemble_evaluations.append(
            {
                "ensemble_id": ensemble["ensemble_id"],
                "family_index": ensemble["family_index"],
                "T": ensemble["T"],
                "aLambda_msbar3": ensemble["aLambda_msbar3"],
                "raw_measure_payload_ref": ensemble_ref,
                "jk_scheme": {
                    "cfg_axis": "delete-1-config",
                    "source_average_inside_cfg": True,
                    "resampling": "jackknife",
                },
                "n_cfg": payload_ensemble.get("n_cfg", ensemble.get("n_cfg")),
                "n_src_per_cfg": payload_ensemble.get("n_src_per_cfg", ensemble.get("n_src_per_cfg")),
                "cfg_ids": payload_ensemble.get("cfg_ids", []),
                "source_descriptors_by_cfg": payload_ensemble.get("source_descriptors_by_cfg", {}),
                "cfg_support_realization_status": (
                    payload_ensemble.get("cfg_seed_status")
                    or payload_ensemble.get("cfg_realization_contract", {}).get("cfg_seed_hash_formula")
                    or "missing"
                ),
                "cfg_source_shape_contract": payload_ensemble.get("cfg_source_shape_contract", {}),
                "tau_int_cfg": None,
                "n_eff_cfg": None,
                "pi_iso": {
                    "cfg_source_corr_t": (payload_ensemble.get("pi_iso") or {}).get("cfg_source_corr_t", []),
                    "cfg_source_corr_shape": (payload_ensemble.get("pi_iso") or {}).get("cfg_source_corr_shape"),
                    "corr_t": ensemble["pi_iso"].get("corr_t", []),
                    "corr_t_stderr": ensemble["pi_iso"].get("corr_t_stderr", []),
                    "corr_t_jk": [],
                    "am_eff_t": ensemble["pi_iso"].get("am_eff_t", []),
                    "am_eff_t_stderr": ensemble["pi_iso"].get("am_eff_t_stderr", []),
                    "am_eff_t_jk": [],
                    "corr_t_candidate_length": t_size,
                    "am_eff_t_candidate_length": max(t_size - 1, 0),
                    "raw_forward_window_candidate_cardinality": raw_forward_window_cardinality,
                    "forward_window_t": [],
                    "forward_window_runs": [],
                    "selected_forward_window": [],
                    "selected_forward_window_cardinality": 0,
                    "forward_window_limit_exists": False,
                    "log_convexity_residual_t": [],
                    "log_convexity_residual_t_stderr": [],
                    "log_convexity_residual_t_jk": [],
                    "log_convexity_upper_t": [],
                    "log_convexity_certified_t": [],
                    "tail_drop_t": [],
                    "tail_drop_t_stderr": [],
                    "tail_drop_t_jk": [],
                    "tail_drop_lower_t": [],
                    "monotone_tail_t": [],
                    "mirror_tail_indicator_t": [],
                    "mirror_tail_indicator_t_stderr": [],
                    "mirror_tail_indicator_t_jk": [],
                    "mirror_to_noise_t": [],
                    "mirror_to_drift_t": [],
                    "mirror_suppressed_t": [],
                    "plateau_noise_floor_t": [],
                    "plateau_flat_t": [],
                    "forward_certificate_t": [],
                    "am_ground_candidate": None,
                    "am_ground_stat_err": None,
                    "am_ground_sys_err": None,
                    "am_ground_candidate_err": None,
                    "ratio_to_lambda_msbar3_candidate": None,
                    "mass_gev_candidate": None,
                    "published_statistical_error": None,
                    "published_systematics": {
                        "status": "pending",
                        "sigma_sys": None,
                        "delta_cont": None,
                        "delta_vol": None,
                        "delta_chi": None,
                    },
                    "convergence_status": "awaiting_measure_evaluation",
                    "sequence_status": "awaiting_measure_evaluation",
                },
                "N_iso": {
                    "cfg_source_corr_direct_t": (payload_ensemble.get("N_iso") or {}).get("cfg_source_corr_direct_t", []),
                    "cfg_source_corr_exchange_t": (payload_ensemble.get("N_iso") or {}).get("cfg_source_corr_exchange_t", []),
                    "cfg_source_corr_t": (payload_ensemble.get("N_iso") or {}).get("cfg_source_corr_t", []),
                    "cfg_source_corr_shape": (payload_ensemble.get("N_iso") or {}).get("cfg_source_corr_shape"),
                    "corr_direct_t": ensemble["N_iso"].get("corr_direct_t", []),
                    "corr_direct_t_stderr": ensemble["N_iso"].get("corr_direct_t_stderr", []),
                    "corr_direct_t_jk": [],
                    "corr_exchange_t": ensemble["N_iso"].get("corr_exchange_t", []),
                    "corr_exchange_t_stderr": ensemble["N_iso"].get("corr_exchange_t_stderr", []),
                    "corr_exchange_t_jk": [],
                    "corr_t": ensemble["N_iso"].get("corr_t", []),
                    "corr_t_stderr": ensemble["N_iso"].get("corr_t_stderr", []),
                    "corr_t_jk": [],
                    "corr_sign_t": [],
                    "corr_sign_t_jk": [],
                    "sign_stable_t": [],
                    "am_eff_t": ensemble["N_iso"].get("am_eff_t", []),
                    "am_eff_t_stderr": ensemble["N_iso"].get("am_eff_t_stderr", []),
                    "am_eff_t_jk": [],
                    "corr_t_candidate_length": t_size,
                    "am_eff_t_candidate_length": max(t_size - 1, 0),
                    "raw_forward_window_candidate_cardinality": raw_forward_window_cardinality,
                    "forward_window_t": [],
                    "forward_window_runs": [],
                    "selected_forward_window": [],
                    "selected_forward_window_cardinality": 0,
                    "forward_window_limit_exists": False,
                    "log_convexity_residual_t": [],
                    "log_convexity_residual_t_stderr": [],
                    "log_convexity_residual_t_jk": [],
                    "log_convexity_upper_t": [],
                    "log_convexity_certified_t": [],
                    "tail_drop_t": [],
                    "tail_drop_t_stderr": [],
                    "tail_drop_t_jk": [],
                    "tail_drop_lower_t": [],
                    "monotone_tail_t": [],
                    "mirror_tail_indicator_t": [],
                    "mirror_tail_indicator_t_stderr": [],
                    "mirror_tail_indicator_t_jk": [],
                    "mirror_to_noise_t": [],
                    "mirror_to_drift_t": [],
                    "mirror_suppressed_t": [],
                    "plateau_noise_floor_t": [],
                    "plateau_flat_t": [],
                    "forward_certificate_t": [],
                    "direct_minus_exchange_residual_t": [],
                    "direct_minus_exchange_residual_t_stderr": [],
                    "direct_minus_exchange_residual_t_jk": [],
                    "direct_minus_exchange_consistent_t": [],
                    "am_ground_candidate": None,
                    "am_ground_stat_err": None,
                    "am_ground_sys_err": None,
                    "am_ground_candidate_err": None,
                    "ratio_to_lambda_msbar3_candidate": None,
                    "mass_gev_candidate": None,
                    "published_statistical_error": None,
                    "published_systematics": {
                        "status": "pending",
                        "sigma_sys": None,
                        "delta_cont": None,
                        "delta_vol": None,
                        "delta_chi": None,
                    },
                    "convergence_status": "awaiting_measure_evaluation",
                    "sequence_status": "awaiting_measure_evaluation",
                },
            }
        )
    return {
        "artifact": "oph_hadron_stable_channel_sequence_evaluator",
        "generated_utc": _timestamp(),
        "status": "awaiting_measure_evaluation",
        "proof_status": "awaiting_measure_evaluation",
        "predictive_promotion_allowed": False,
        "source_artifact": cfg_source_payload.get("artifact", sequence_population.get("artifact")),
        "measure_evaluation_law_id": "oph_qcd_2p1_stable_channel_cfg_source_jackknife_evaluation",
        "theorem_candidate": "StableChannelForwardWindowConvergence",
        "forward_window_certificate_family": "oph_hadron_forward_window_certificate_family",
        "support_realization_schedule": cfg_source_payload.get("support_realization_schedule"),
        "runtime_receipt_artifact": runtime_receipt.get("artifact"),
        "next_single_residual_object_after_cfg_execution": "oph_hadron_stable_channel_sequence_evaluator",
        "pi_formulae": {
            "corr_estimator": "(1/(n_cfg*n_src_per_cfg)) * sum_{cfg,src} c_pi^{cfg,src}(t)",
            "cfg_average_inside_cfg": "(1/|C_n|) * sum_cfg ((1/|S_{n,cfg}|) * sum_src c_pi^{cfg,src}(t))",
            "am_eff": "log(corr_t[t] / corr_t[t+1])",
            "forward_window": "{t : 1 <= t+1 < floor(T/2)}",
            "log_convexity_residual": "corr_t[t]^2 - corr_t[t-1]*corr_t[t+1]",
            "tail_drop": "am_eff_t[t] - am_eff_t[t+1]",
            "mirror_tail_indicator": "exp(-am_eff_t[t] * (T - 2*t))"
        },
        "n_formulae": {
            "corr_estimator": "(1/(n_cfg*n_src_per_cfg)) * sum_{cfg,src} (c_dir^{cfg,src}(t) - c_ex^{cfg,src}(t))",
            "cfg_average_inside_cfg": "(1/|C_n|) * sum_cfg ((1/|S_{n,cfg}|) * sum_src (c_dir^{cfg,src}(t) - c_ex^{cfg,src}(t)))",
            "am_eff": "log(|corr_t[t]| / |corr_t[t+1]|)",
            "forward_window": "{t : 1 <= t+1 < floor(T/2)}",
            "log_convexity_residual": "|corr_t[t]|^2 - |corr_t[t-1]|*|corr_t[t+1]|",
            "tail_drop": "am_eff_t[t] - am_eff_t[t+1]",
            "mirror_tail_indicator": "exp(-am_eff_t[t] * (T - 2*t))"
        },
        "ensemble_evaluations": ensemble_evaluations,
        "smallest_constructive_missing_object": cfg_source_payload.get(
            "smallest_constructive_missing_object",
            "runtime_schedule_receipt_N_therm_and_N_sep",
        ),
        "notes": [
            "The stable-channel evaluation law is fixed at cfg/source jackknife level; this artifact tracks the per-ensemble arrays needed before forward-window convergence can certify masses.",
            "Once the external runtime schedule receipt is supplied and the executed schedule writes the cfg/source arrays, this evaluator is the next single residual object before forward-window convergence can certify masses.",
            "Masses remain unset until the realized cfg/source arrays, evaluation arrays, and convergence certificates are populated."
        ]
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the stable-channel sequence-evaluation artifact.")
    parser.add_argument("--sequence-population", default=str(DEFAULT_SEQUENCE_POPULATION))
    parser.add_argument("--cfg-source-payload", default=str(DEFAULT_CFG_SOURCE_PAYLOAD))
    parser.add_argument("--runtime-receipt", default=str(DEFAULT_RUNTIME_RECEIPT))
    parser.add_argument(
        "--backend-dump",
        default=None,
        help="Optional validated production dump to evaluate into jackknife/forward-window outputs.",
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    sequence_population = json.loads(Path(args.sequence_population).read_text(encoding="utf-8"))
    cfg_source_payload_path = Path(args.cfg_source_payload)
    cfg_source_payload = (
        json.loads(cfg_source_payload_path.read_text(encoding="utf-8"))
        if cfg_source_payload_path.exists()
        else None
    )
    runtime_receipt_path = Path(args.runtime_receipt)
    runtime_receipt = (
        json.loads(runtime_receipt_path.read_text(encoding="utf-8"))
        if runtime_receipt_path.exists()
        else None
    )
    artifact = build_artifact(sequence_population, cfg_source_payload, runtime_receipt)
    if args.backend_dump:
        backend_dump = json.loads(Path(args.backend_dump).read_text(encoding="utf-8"))
        artifact = populate_evaluation_from_dump(artifact, backend_dump)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
