#!/usr/bin/env python3
"""Smoke-test the stable-channel hadron ground-state readout artifact."""

from __future__ import annotations

import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "particles" / "runs" / "hadron" / "stable_channel_groundstate_readout.json"


def main() -> int:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_hadron_stable_channel_groundstate_readout":
        print("wrong hadron stable-channel artifact id", file=sys.stderr)
        return 1
    if payload.get("proof_status") != "candidate_only":
        print("stable-channel readout should remain candidate_only", file=sys.stderr)
        return 1
    channels = payload.get("channels", {})
    if {"pi_iso", "N_iso"} - set(channels):
        print("stable-channel artifact should cover pi_iso and N_iso", file=sys.stderr)
        return 1
    if payload.get("theorem_candidate") != "StableChannelForwardWindowConvergence":
        print("stable-channel artifact should expose the forward-window convergence theorem candidate", file=sys.stderr)
        return 1
    if channels["pi_iso"].get("ratio_to_lambda_msbar3", "missing") is not None:
        print("pi_iso ratio should remain unset until convergence closes", file=sys.stderr)
        return 1
    if channels["N_iso"].get("ratio_to_lambda_msbar3", "missing") is not None:
        print("N_iso ratio should remain unset until convergence closes", file=sys.stderr)
        return 1
    if not channels["pi_iso"].get("per_ensemble") or not channels["N_iso"].get("per_ensemble"):
        print("stable-channel artifact should carry per-ensemble channel families", file=sys.stderr)
        return 1
    promoted = set(channels["N_iso"].get("promoted_channel_fields", []))
    if {"corr_direct_t", "corr_exchange_t", "ratio_to_lambda_msbar3"} - promoted:
        print("N_iso promoted field list should include the stable-channel payload from e152", file=sys.stderr)
        return 1
    frontier = set(payload.get("minimal_closure_frontier", []))
    if "stable_channel_groundstate_readout" not in frontier:
        print("stable-channel frontier should be tracked explicitly", file=sys.stderr)
        return 1
    availability = payload.get("data_availability", {})
    if availability.get("full_unquenched_correlator_status") != "predictive_ensemble_seeded_candidate":
        print("artifact should point to the seeded unquenched correlator producer", file=sys.stderr)
        return 1
    upstream = payload.get("upstream", {})
    if upstream.get("stable_channel_sequence_evaluation_status") != "awaiting_measure_evaluation":
        print("artifact should point to the new stable-channel sequence evaluator stage", file=sys.stderr)
        return 1
    if availability.get("raw_correlator_arrays_present") is not False:
        print("stable-channel arrays should remain empty until the unquenched producer is populated", file=sys.stderr)
        return 1
    if availability.get("effective_mass_sequences_present") is not False:
        print("effective-mass sequences should remain empty until the producer emits them", file=sys.stderr)
        return 1
    if availability.get("full_baryon_contractions_present") is not True:
        print("stable-channel artifact should treat the direct-minus-exchange baryon contraction law as closed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
