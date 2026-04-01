#!/usr/bin/env python3
"""Emit the eventual spectral-floor side condition beneath the UV/BW faithfulness term."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_eventual_spectral_floor_scaffold.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_payload() -> dict[str, object]:
    return {
        "artifact": "oph_bw_fixed_local_collar_eventual_spectral_floor_scaffold",
        "generated_utc": _timestamp(),
        "status": "minimal_faithfulness_side_extension",
        "public_promotion_allowed": False,
        "exact_missing_object": "eventual_fixed_local_collar_spectral_floor_for_transported_marginals",
        "role": (
            "Package the sole nonlatent faithfulness-side clause beneath the fixed-local-collar "
            "faithful modular-defect witness."
        ),
        "contract": {
            "for_fixed_models": "every fixed local collar model (m, delta)",
            "must_emit": (
                "exists lambda_bar_{m,delta} > 0 with lambda_{*,n,m,delta} >= "
                "lambda_bar_{m,delta} eventually"
            ),
        },
        "meaning": (
            "Along the realized refinement chain, each fixed transported local-collar marginal stays "
            "uniformly faithful with a collarwise eigenvalue floor bounded away from zero."
        ),
        "unlocks": [
            "fixed_local_collar_faithful_modular_defect_vanishing",
            "vanishing_carried_collar_schedule_on_fixed_local_collars",
            "canonical_scaling_cap_pair_realization_from_transported_cap_marginals",
        ],
        "notes": [
            "This clause does not by itself emit the carried-collar schedule or the scaling-limit cap pair.",
            "It is the live faithfulness-side clause still external to the emitted theorem chain.",
            "On the local-Gibbs plus exponential-mixing pullback branch, the recovery/Markov side is already latent once epsilon_{n,m,delta} -> 0.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the fixed-local-collar eventual spectral-floor scaffold."
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(build_payload(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
