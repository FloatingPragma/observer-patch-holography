#!/usr/bin/env python3
"""Smoke-test the hadron runtime schedule receipt artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
FULL_SCRIPT = ROOT / "particles" / "hadron" / "derive_full_unquenched_correlator.py"
POP_SCRIPT = ROOT / "particles" / "hadron" / "derive_stable_channel_sequence_population.py"
PAYLOAD_SCRIPT = ROOT / "particles" / "hadron" / "derive_stable_channel_cfg_source_measure_payload.py"
SCRIPT = ROOT / "particles" / "hadron" / "derive_runtime_schedule_receipt_n_therm_and_n_sep.py"
def main() -> int:
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = pathlib.Path(tmp_dir)
        payload_path = tmp_path / "payload.json"
        receipt_path = tmp_path / "receipt.json"
        subprocess.run([sys.executable, str(FULL_SCRIPT)], check=True, cwd=ROOT)
        subprocess.run([sys.executable, str(POP_SCRIPT)], check=True, cwd=ROOT)
        subprocess.run(
            [
                sys.executable,
                str(PAYLOAD_SCRIPT),
                "--runtime-receipt",
                str(tmp_path / "missing_receipt.json"),
                "--output",
                str(payload_path),
            ],
            check=True,
            cwd=ROOT,
        )
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--payload",
                str(payload_path),
                "--output",
                str(receipt_path),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(receipt_path.read_text(encoding="utf-8"))
    if payload.get("artifact") != "runtime_schedule_receipt_N_therm_and_N_sep":
        print("wrong hadron runtime receipt artifact id", file=sys.stderr)
        return 1
    if payload.get("status") != "awaiting_external_runtime_inputs":
        print("hadron runtime receipt should wait on external runtime inputs", file=sys.stderr)
        return 1
    required = payload.get("required_schedule_scalars", {})
    if required.get("N_therm", "sentinel") is not None or required.get("N_sep", "sentinel") is not None:
        print("runtime receipt should keep schedule scalars unset until supplied externally", file=sys.stderr)
        return 1
    if payload.get("execution_contract", {}).get("next_single_residual_object_after_execution") != "oph_hadron_stable_channel_sequence_evaluator":
        print("runtime receipt should hand off to the stable-channel sequence evaluator after execution", file=sys.stderr)
        return 1
    if payload.get("completion_mode") != "execution_and_systematics_contract":
        print("runtime receipt should classify the hadron frontier as execution/systematics work", file=sys.stderr)
        return 1
    if payload.get("execution_and_systematics_contract", {}).get("systematics_contract", {}).get("sigma_sys_formula") != "sqrt(delta_cont^2 + delta_vol^2 + delta_chi^2)":
        print("runtime receipt should expose the published hadron systematics contract", file=sys.stderr)
        return 1
    publication_fields = payload.get("execution_and_systematics_contract", {}).get("provenance_contract", {}).get("required_publication_fields", [])
    expected_fields = {
        "stable_channel_sequence_evaluation.ensemble_evaluations[*].pi_iso.published_systematics",
        "stable_channel_sequence_evaluation.ensemble_evaluations[*].pi_iso.published_statistical_error",
        "stable_channel_sequence_evaluation.ensemble_evaluations[*].N_iso.published_systematics",
        "stable_channel_sequence_evaluation.ensemble_evaluations[*].N_iso.published_statistical_error",
    }
    if set(publication_fields) != expected_fields:
        print("runtime receipt should point publication requirements at channel-local stable-channel fields", file=sys.stderr)
        return 1
    for entry in payload.get("execution_contract", {}).get("ensemble_schedule", []):
        if not isinstance(entry.get("t_extent"), int) or entry.get("t_extent", 0) <= 0:
            print("runtime receipt should expose positive t_extent values for every scheduled ensemble", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
