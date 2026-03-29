#!/usr/bin/env python3
"""Smoke-test the stable-channel sequence-evaluation artifact."""

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
RECEIPT_SCRIPT = ROOT / "particles" / "hadron" / "derive_runtime_schedule_receipt_n_therm_and_n_sep.py"
EVAL_SCRIPT = ROOT / "particles" / "hadron" / "derive_stable_channel_sequence_evaluation.py"
def main() -> int:
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = pathlib.Path(tmp_dir)
        payload_path = tmp_path / "payload.json"
        receipt_path = tmp_path / "receipt.json"
        eval_path = tmp_path / "evaluation.json"
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
                str(RECEIPT_SCRIPT),
                "--payload",
                str(payload_path),
                "--output",
                str(receipt_path),
            ],
            check=True,
            cwd=ROOT,
        )
        subprocess.run(
            [
                sys.executable,
                str(PAYLOAD_SCRIPT),
                "--runtime-receipt",
                str(receipt_path),
                "--output",
                str(payload_path),
            ],
            check=True,
            cwd=ROOT,
        )
        subprocess.run(
            [
                sys.executable,
                str(EVAL_SCRIPT),
                "--cfg-source-payload",
                str(payload_path),
                "--runtime-receipt",
                str(receipt_path),
                "--output",
                str(eval_path),
            ],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(eval_path.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_hadron_stable_channel_sequence_evaluator":
        print("wrong stable-channel sequence-evaluation artifact id", file=sys.stderr)
        return 1
    if payload.get("status") != "awaiting_measure_evaluation":
        print("sequence-evaluation artifact should record measure evaluation as the next step", file=sys.stderr)
        return 1
    if payload.get("measure_evaluation_law_id") != "oph_qcd_2p1_stable_channel_cfg_source_jackknife_evaluation":
        print("sequence-evaluation artifact should expose the cfg/source jackknife evaluation law", file=sys.stderr)
        return 1
    if payload.get("theorem_candidate") != "StableChannelForwardWindowConvergence":
        print("sequence-evaluation artifact should expose the forward-window convergence theorem", file=sys.stderr)
        return 1
    evaluations = payload.get("ensemble_evaluations", [])
    if not evaluations:
        print("sequence-evaluation artifact should expose per-ensemble evaluation shells", file=sys.stderr)
        return 1
    if not str(evaluations[0].get("raw_measure_payload_ref", "")).startswith("oph_hadron_stable_channel_cfg_source_measure_payload::"):
        print("sequence-evaluation artifact should point raw payload refs at the cfg/source payload artifact", file=sys.stderr)
        return 1
    if "forward_window_t" not in evaluations[0]["pi_iso"] or "forward_window_t" not in evaluations[0]["N_iso"]:
        print("sequence-evaluation artifact should expose forward-window placeholders", file=sys.stderr)
        return 1
    if "corr_t_jk" not in evaluations[0]["pi_iso"] or "corr_t_jk" not in evaluations[0]["N_iso"]:
        print("sequence-evaluation artifact should expose jackknife correlator placeholders", file=sys.stderr)
        return 1
    if evaluations[0].get("cfg_support_realization_status") == "missing":
        print("sequence-evaluation artifact should expose the cfg support realization status", file=sys.stderr)
        return 1
    if "forward_certificate_t" not in evaluations[0]["pi_iso"] or "forward_certificate_t" not in evaluations[0]["N_iso"]:
        print("sequence-evaluation artifact should expose forward-window certificate placeholders", file=sys.stderr)
        return 1
    if evaluations[0]["pi_iso"].get("raw_forward_window_candidate_cardinality", 0) <= 0:
        print("sequence-evaluation artifact should expose a nonempty raw forward-window candidate size", file=sys.stderr)
        return 1
    if evaluations[0]["pi_iso"].get("corr_t_candidate_length", 0) <= 0:
        print("sequence-evaluation artifact should expose candidate correlator lengths", file=sys.stderr)
        return 1
    if payload.get("forward_window_certificate_family") != "oph_hadron_forward_window_certificate_family":
        print("sequence-evaluation artifact should name the forward-window certificate family", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "runtime_schedule_receipt_N_therm_and_N_sep":
        print("sequence-evaluation artifact should inherit the external runtime schedule receipt frontier", file=sys.stderr)
        return 1
    if payload.get("runtime_receipt_artifact") != "runtime_schedule_receipt_N_therm_and_N_sep":
        print("sequence-evaluation artifact should point at the runtime receipt artifact", file=sys.stderr)
        return 1
    for channel_name in ("pi_iso", "N_iso"):
        channel = evaluations[0][channel_name]
        if "published_statistical_error" not in channel or "published_systematics" not in channel:
            print("sequence-evaluation artifact should expose channel-local publication placeholders", file=sys.stderr)
            return 1
        if (channel.get("published_systematics") or {}).get("status") != "pending":
            print("sequence-evaluation artifact should mark channel-local published systematics as pending before execution", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
