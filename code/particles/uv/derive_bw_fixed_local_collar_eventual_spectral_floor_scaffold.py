#!/usr/bin/env python3
"""Emit the eventual spectral-floor side condition beneath the UV/BW faithfulness term."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from bw_collar_honesty import build_comparison_reference_floor_transfer


ROOT = Path(__file__).resolve().parents[2]
EXACT_MARKOV_MODULUS = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_exact_markov_modulus_scaffold.json"
)
COMMON_FLOOR = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_modular_transport_common_floor_scaffold.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_eventual_spectral_floor_scaffold.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _artifact_ref(path: Path) -> str:
    return f"code/{path.relative_to(ROOT).as_posix()}"


def build_payload() -> dict[str, object]:
    return {
        "artifact": "oph_bw_fixed_local_collar_eventual_spectral_floor_scaffold",
        "generated_utc": _timestamp(),
        "status": "legacy_coarse_wrapper",
        "public_promotion_allowed": False,
        "exact_missing_object": "eventual_fixed_local_collar_spectral_floor_for_transported_marginals",
        "role": (
            "Package the legacy coarse all-marginal spectral-floor wrapper above the smaller "
            "modular-transport common-floor clause actually consumed by the proof."
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
        "coarsens_live_artifact": _artifact_ref(COMMON_FLOOR),
        "live_exact_missing_object": "eventual_fixed_local_collar_common_floor_on_modular_transport_marginals",
        "comparison_reference_floor_transfer": build_comparison_reference_floor_transfer(
            exact_markov_artifact=_artifact_ref(EXACT_MARKOV_MODULUS),
            spectral_floor_artifact=_artifact_ref(COMMON_FLOOR),
        ),
        "unlocks": [
            "fixed_local_collar_faithful_modular_defect_vanishing",
            "vanishing_carried_collar_schedule_on_fixed_local_collars",
            "canonical_scaling_cap_pair_realization_from_transported_cap_marginals",
        ],
        "notes": [
            "This clause is intentionally larger than the live primitive blocker.",
            "The faithful modular-defect proof only consumes the finite modular-transport family Xi^{mod}_{m,delta}, so the smaller common-floor artifact is the live target and this all-marginal floor is only a legacy coarse wrapper.",
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
