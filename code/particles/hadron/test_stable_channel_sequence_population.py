#!/usr/bin/env python3
"""Smoke-test the stable-channel sequence-population artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FULL_SCRIPT = ROOT / "particles" / "hadron" / "derive_full_unquenched_correlator.py"
POP_SCRIPT = ROOT / "particles" / "hadron" / "derive_stable_channel_sequence_population.py"
OUTPUT = ROOT / "particles" / "runs" / "hadron" / "stable_channel_sequence_population.json"


def main() -> int:
    subprocess.run([sys.executable, str(FULL_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(POP_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_hadron_stable_channel_sequence_population":
        print("wrong stable-channel sequence-population artifact id", file=sys.stderr)
        return 1
    if payload.get("status") != "law_closed_waiting_measure_evaluation":
        print("sequence-population artifact should now record the closed law waiting on evaluation", file=sys.stderr)
        return 1
    if payload.get("predictive_promotion_allowed") is not False:
        print("sequence-population artifact must remain non-promoted", file=sys.stderr)
        return 1
    sequences = payload.get("ensemble_sequences", [])
    if not sequences:
        print("sequence-population artifact should expose per-ensemble sequence shells", file=sys.stderr)
        return 1
    if "t_support" not in sequences[0] or "t_lambda_msbar3" not in sequences[0]:
        print("sequence-population artifact should expose per-ensemble time support", file=sys.stderr)
        return 1
    if set(sequences[0]) & {"mass_gev", "ratio_to_lambda_msbar3", "am_ground"}:
        print("sequence-population artifact must not promote masses before convergence", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
