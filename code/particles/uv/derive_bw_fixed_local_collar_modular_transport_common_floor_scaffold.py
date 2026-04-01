#!/usr/bin/env python3
"""Emit the exact modular-transport common-floor scaffold beneath UV/BW faithfulness."""

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
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_modular_transport_common_floor_scaffold.json"
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _artifact_ref(path: Path) -> str:
    return f"code/{path.relative_to(ROOT).as_posix()}"


def build_payload() -> dict[str, object]:
    return {
        "artifact": "oph_bw_fixed_local_collar_modular_transport_common_floor_scaffold",
        "generated_utc": _timestamp(),
        "status": "minimal_faithfulness_side_extension",
        "public_promotion_allowed": False,
        "exact_missing_object": "eventual_fixed_local_collar_common_floor_on_modular_transport_marginals",
        "role": (
            "Package the sole nonlatent faithfulness-side clause beneath the fixed-local-collar "
            "faithful modular-defect witness, restricted to the finite modular-transport family "
            "that the proof actually consumes."
        ),
        "contract": {
            "for_fixed_models": "every fixed local collar model (m, delta)",
            "relevant_family": "Xi^{mod}_{m,delta}",
            "must_emit": (
                "exists lambda_bar_{m,delta} > 0 and N_{m,delta} such that for all n >= N_{m,delta} "
                "and every X in Xi^{mod}_{m,delta}, rho_{n->m,X} >= lambda_bar_{m,delta} * 1"
            ),
        },
        "meaning": (
            "Along the realized refinement chain, the finitely many modular-transport marginals "
            "used by the faithful modular-defect estimate stay uniformly faithful with one common "
            "collarwise eigenvalue floor bounded away from zero."
        ),
        "comparison_reference_floor_transfer": build_comparison_reference_floor_transfer(
            exact_markov_artifact=_artifact_ref(EXACT_MARKOV_MODULUS),
            spectral_floor_artifact=_artifact_ref(DEFAULT_OUT),
        ),
        "unlocks": [
            "fixed_local_collar_faithful_modular_defect_vanishing",
            "vanishing_carried_collar_schedule_on_fixed_local_collars",
            "canonical_scaling_cap_pair_realization_from_transported_cap_marginals",
        ],
        "notes": [
            "This clause is smaller than an all-marginal transported spectral floor because the faithful modular-defect proof only consumes the finite family Xi^{mod}_{m,delta}.",
            "It does not by itself emit the carried-collar schedule or the scaling-limit cap pair.",
            "When paired with the exact-Markov modulus witness, this same floor transfers to the exact-Markov comparison marginals on sufficiently late stages, so no second comparison-state spectral input remains hidden below the faithful modular-defect step.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the fixed-local-collar modular-transport common-floor scaffold."
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
