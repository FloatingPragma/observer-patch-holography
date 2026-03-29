#!/usr/bin/env python3
"""Smoke-test the current hadron-lane audit artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FULL_SCRIPT = ROOT / "particles" / "hadron" / "derive_full_unquenched_correlator.py"
POP_SCRIPT = ROOT / "particles" / "hadron" / "derive_stable_channel_sequence_population.py"
PAYLOAD_SCRIPT = ROOT / "particles" / "hadron" / "derive_stable_channel_cfg_source_measure_payload.py"
RECEIPT_SCRIPT = ROOT / "particles" / "hadron" / "derive_runtime_schedule_receipt_n_therm_and_n_sep.py"
EVAL_SCRIPT = ROOT / "particles" / "hadron" / "derive_stable_channel_sequence_evaluation.py"
SURROGATE_SCRIPT = ROOT / "particles" / "hadron" / "derive_hadron_surrogate_execution_bridge_status.py"
GEOMETRY_SCRIPT = ROOT / "particles" / "hadron" / "derive_hadron_production_geometry_summary.py"
AUDIT_SCRIPT = ROOT / "particles" / "hadron" / "derive_current_hadron_lane_audit.py"
AUDIT = ROOT / "particles" / "runs" / "hadron" / "current_hadron_lane_audit.json"


def main() -> int:
    subprocess.run([sys.executable, str(FULL_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(POP_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PAYLOAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(RECEIPT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PAYLOAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EVAL_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SURROGATE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(GEOMETRY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(AUDIT_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(AUDIT.read_text(encoding="utf-8"))
    if payload.get("promotion_verdict") != "suppress_from_public_surface":
        print("hadron audit should suppress the current public hadron surface", file=sys.stderr)
        return 1
    blockers = set(payload.get("pipeline_classification", {}).get("blockers", []))
    required = {
        "stable_channel_cfg_source_measure_payload",
        "stable_channel_sequence_evaluation",
        "StableChannelForwardWindowConvergence",
        "finite_volume_resonance_and_spectrum_readout",
    }
    if not required.issubset(blockers):
        print(f"missing expected hadron blockers: {sorted(required - blockers)}", file=sys.stderr)
        return 1
    frontier = set(payload.get("minimal_closure_frontier", []))
    if {"full_unquenched_correlator", "stable_channel_sequence_population", "stable_channel_sequence_evaluation", "stable_channel_groundstate_readout", "finite_volume_resonance_and_spectrum_readout"} - frontier:
        print("hadron audit should expose the sharpened minimal closure frontier", file=sys.stderr)
        return 1
    next_artifact = payload.get("recommended_next_predictive_artifact", {}).get("name")
    if next_artifact != "backend_correlator_dump.production.json":
        print("hadron audit should point to the real production backend correlator dump after the runtime receipt exists", file=sys.stderr)
        return 1
    if payload.get("smallest_constructive_missing_object") != "backend_correlator_dump.production.json from real production RHMC/HMC execution on the theorem-emitted seeded family":
        print("hadron audit should reduce to the real production backend correlator dump once the receipt is explicit", file=sys.stderr)
        return 1
    next_missing = payload.get("pipeline_classification", {}).get("summary", {}).get("next_missing_object")
    if next_missing != "backend_correlator_dump.production.json from real production RHMC/HMC execution on the theorem-emitted seeded family":
        print("hadron audit summary should agree that the next missing object is the real production backend correlator dump", file=sys.stderr)
        return 1
    if payload.get("smallest_missing_theorem_after_full_unquenched") != "StableChannelForwardWindowConvergence":
        print("hadron audit should identify forward-window convergence as the next theorem after full unquenching", file=sys.stderr)
        return 1
    notes = " ".join(payload.get("notes", []))
    if "execution-and-systematics contract" not in notes or "backend correlator dump" not in notes:
        print("hadron audit should describe both the runtime contract and the production dump frontier", file=sys.stderr)
        return 1
    surrogate = payload.get("surrogate_execution_bridge") or {}
    if surrogate.get("status") != "surrogate_hmc_execution_bridge_complete":
        print("hadron audit should surface the surrogate bridge as a separate diagnostic layer", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
