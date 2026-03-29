#!/usr/bin/env python3
"""Validate production-side hadron closure inputs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def _load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _is_finite_number(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def _get_schedule_scalars(receipt: dict[str, Any]) -> tuple[Any, Any]:
    scalars = receipt.get("required_schedule_scalars", {})
    if not scalars and isinstance(receipt.get("execution_contract"), dict):
        scalars = receipt["execution_contract"].get("required_schedule_scalars", {})
    return scalars.get("N_therm"), scalars.get("N_sep")


def _channel_summary(channel: dict[str, Any]) -> dict[str, Any]:
    published = channel.get("published_systematics")
    sigma_sys = published.get("sigma_sys") if isinstance(published, dict) else None
    return {
        "has_cfg_source_corr_t": isinstance(channel.get("cfg_source_corr_t"), list) and len(channel.get("cfg_source_corr_t")) > 0,
        "am_ground_candidate_finite": _is_finite_number(channel.get("am_ground_candidate")),
        "mass_gev_candidate_finite": _is_finite_number(channel.get("mass_gev_candidate")),
        "published_statistical_error_finite": _is_finite_number(channel.get("published_statistical_error")),
        "published_systematics_sigma_sys_finite": _is_finite_number(sigma_sys),
        "forward_window_limit_exists": bool(channel.get("forward_window_limit_exists")),
        "selected_forward_window_nonempty": isinstance(channel.get("selected_forward_window"), list) and len(channel.get("selected_forward_window")) > 0,
        "systematics_status": published.get("status") if isinstance(published, dict) else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate current-boundary OPH hadron closure status.")
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--dump", required=False)
    parser.add_argument("--evaluation", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    receipt = _load_json(args.receipt)
    dump = _load_json(args.dump) if args.dump else None
    evaluation = _load_json(args.evaluation)
    n_therm, n_sep = _get_schedule_scalars(receipt)
    receipt_filled = _is_finite_number(n_therm) and _is_finite_number(n_sep)

    ensembles = {}
    channel_checks = []
    for entry in evaluation.get("ensemble_evaluations", []):
        ens_id = entry.get("ensemble_id", "<unknown>")
        ensembles[ens_id] = {}
        for channel_name in ("pi_iso", "N_iso"):
            summary = _channel_summary(entry.get(channel_name, {}))
            ensembles[ens_id][channel_name] = summary
            channel_checks.extend(
                [
                    summary["has_cfg_source_corr_t"],
                    summary["am_ground_candidate_finite"],
                    summary["mass_gev_candidate_finite"],
                    summary["published_statistical_error_finite"],
                    summary["published_systematics_sigma_sys_finite"],
                    summary["forward_window_limit_exists"],
                    summary["selected_forward_window_nonempty"],
                ]
            )

    all_channels_ok = bool(channel_checks) and all(channel_checks)
    production_dump_present = dump is not None and dump.get("production_execution") is True
    public_ready = all_channels_ok and production_dump_present and receipt_filled
    if public_ready:
        residual = None
    elif not receipt_filled:
        residual = "runtime_schedule_receipt_N_therm_and_N_sep with explicit external N_therm/N_sep inputs"
    elif not production_dump_present:
        residual = "backend_correlator_dump.production.json from real production RHMC/HMC execution on the theorem-emitted seeded family"
    else:
        residual = (
            "stable_channel_sequence_evaluation with populated forward-window and published statistical/systematic fields "
            "for pi_iso and N_iso"
        )
    report = {
        "artifact": "oph_hadron_closure_validation_report",
        "required_schedule_scalars": {"N_therm": n_therm, "N_sep": n_sep},
        "production_dump_present": production_dump_present,
        "receipt_filled": receipt_filled,
        "all_channels_ok": all_channels_ok,
        "closure_grade": "public_unsuppression_ready" if public_ready else "execution_incomplete",
        "public_unsuppression_ready": public_ready,
        "smallest_live_residual_object": residual,
        "ensembles": ensembles,
    }
    Path(args.output).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
