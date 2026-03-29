#!/usr/bin/env python3
"""Render the hadron runtime/systematics status Markdown surface."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RECEIPT = ROOT / "particles" / "runs" / "hadron" / "runtime_schedule_receipt_N_therm_and_N_sep.json"
DEFAULT_DUMP = ROOT / "particles" / "runs" / "hadron" / "backend_correlator_dump.production.json"
DEFAULT_EVALUATION = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_evaluation.json"
DEFAULT_CLOSURE = ROOT / "particles" / "runs" / "hadron" / "hadron_production_closure_validation_report.json"
DEFAULT_OUT = ROOT / "particles" / "HADRON_SYSTEMATICS_STATUS.md"


def _load_optional_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _channel_row(ensemble_id: str, channel_name: str, channel: dict[str, Any]) -> str:
    published = channel.get("published_systematics") or {}
    return (
        f"| `{ensemble_id}` | `{channel_name}` | `{channel.get('forward_window_limit_exists', False)}`"
        f" | `{channel.get('am_ground_candidate')}` | `{channel.get('published_statistical_error')}`"
        f" | `{published.get('sigma_sys')}` | `{channel.get('mass_gev_candidate')}` |"
    )


def build_markdown(
    receipt: dict[str, Any],
    evaluation: dict[str, Any],
    closure: dict[str, Any] | None,
    dump: dict[str, Any] | None,
) -> str:
    required = receipt.get("required_schedule_scalars") or {}
    dump_ensembles = (dump or {}).get("ensembles") or {}
    lines = [
        "# Hadron Pipeline Status",
        "",
        "This file tracks the execution and publication contract for the active hadron",
        "prediction path separately from the main particle table.",
        "",
        "## Execution Receipt",
        f"- artifact: `{receipt.get('artifact')}`",
        f"- receipt_kind: `{receipt.get('receipt_kind')}`",
        f"- status: `{receipt.get('status')}`",
        f"- `N_therm`: `{required.get('N_therm')}`",
        f"- `N_sep`: `{required.get('N_sep')}`",
        f"- kernel_id: `{(receipt.get('execution_contract') or {}).get('kernel_id')}`",
        f"- initial_configuration: `{(receipt.get('execution_contract') or {}).get('initial_configuration')}`",
        f"- next_single_residual_object_after_execution: `{(receipt.get('execution_contract') or {}).get('next_single_residual_object_after_execution')}`",
        "",
        "## Stable-Channel Execution State",
        "| ensemble_id | cfg_ids | n_cfg | n_src_per_cfg | t_extent | arrays_written |",
        "| --- | --- | ---: | ---: | ---: | --- |",
    ]
    for schedule in (receipt.get("execution_contract") or {}).get("ensemble_schedule", []):
        ensemble_id = str(schedule["ensemble_id"])
        dump_entry = dump_ensembles.get(ensemble_id) or {}
        arrays_written = bool(dump_entry.get("cfgs"))
        lines.append(
            f"| `{ensemble_id}` | `{','.join(schedule.get('cfg_ids', []))}` | `{len(schedule.get('cfg_ids', []))}`"
            f" | `{schedule.get('n_src_per_cfg')}` | `{schedule.get('t_extent')}` | `{'yes' if arrays_written else 'no'}` |"
        )
    lines.extend(
        [
            "",
            "## Stable-Channel Numerical State",
            "| ensemble_id | channel | forward_window_limit_exists | am_ground_candidate | stat_err | sys_err | mass_gev_candidate |",
            "| --- | --- | --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for ensemble in evaluation.get("ensemble_evaluations", []):
        ensemble_id = str(ensemble["ensemble_id"])
        lines.append(_channel_row(ensemble_id, "pi_iso", ensemble.get("pi_iso") or {}))
        lines.append(_channel_row(ensemble_id, "N_iso", ensemble.get("N_iso") or {}))
    lines.extend(
        [
            "",
            "## Published Systematics Budget",
            "| ensemble_id | channel | statistics | continuum | volume | chiral | publishable |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    public_ready = bool((closure or {}).get("public_unsuppression_ready"))
    for ensemble in evaluation.get("ensemble_evaluations", []):
        ensemble_id = str(ensemble["ensemble_id"])
        for channel_name in ("pi_iso", "N_iso"):
            channel = ensemble.get(channel_name) or {}
            published = channel.get("published_systematics") or {}
            publishable = (
                public_ready
                and channel.get("forward_window_limit_exists")
                and channel.get("published_statistical_error") is not None
                and published.get("sigma_sys") is not None
            )
            status = published.get("status") or "pending"
            lines.append(
                f"| `{ensemble_id}` | `{channel_name}` | `{('complete' if channel.get('published_statistical_error') is not None else 'pending')}`"
                f" | `{status}` | `{status}` | `{status}` | `{'yes' if publishable else 'no'}` |"
            )
    lines.extend(
        [
            "",
            "## Residual Open Object",
            f"- current: `{(closure or {}).get('smallest_live_residual_object')}`",
            "- after-execution: `StableChannelForwardWindowConvergence`",
            "- after-convergence: `oph_hadron_rho_scattering_readout`",
            "",
            "## Working Rule",
            "Keep hadron rows out of the public particle table until the active pipeline",
            "has a real production backend dump plus published statistical and",
            "continuum/volume/chiral systematics on the seeded `2+1`, QED-off stable-channel",
            "branch.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Render HADRON_SYSTEMATICS_STATUS.md from live artifacts.")
    parser.add_argument("--receipt", default=str(DEFAULT_RECEIPT))
    parser.add_argument("--dump", default=str(DEFAULT_DUMP))
    parser.add_argument("--evaluation", default=str(DEFAULT_EVALUATION))
    parser.add_argument("--closure", default=str(DEFAULT_CLOSURE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    receipt = json.loads(Path(args.receipt).read_text(encoding="utf-8"))
    evaluation = json.loads(Path(args.evaluation).read_text(encoding="utf-8"))
    closure = _load_optional_json(Path(args.closure))
    dump = _load_optional_json(Path(args.dump))

    markdown = build_markdown(receipt, evaluation, closure, dump)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(markdown, encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
