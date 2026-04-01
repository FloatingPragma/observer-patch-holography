#!/usr/bin/env python3
"""Emit the faithful modular-defect scaffold beneath the carried-collar schedule."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bw_collar_honesty import (
    build_comparison_reference_floor_transfer,
    build_schedule_term_frontier,
)


ROOT = Path(__file__).resolve().parents[2]
EXTRACTION_SCAFFOLD = ROOT / "particles" / "runs" / "uv" / "bw_scaling_limit_cap_pair_extraction_scaffold.json"
RAW_DATUM = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_markov_faithfulness_datum.json"
CONSTRUCTIVE_RECOVERY = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_constructive_recovery_scaffold.json"
)
EXACT_MARKOV_MODULUS = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_exact_markov_modulus_scaffold.json"
)
COMMON_FLOOR = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_modular_transport_common_floor_scaffold.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_faithful_modular_defect_scaffold.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_ref(path: Path) -> str:
    return f"code/{path.relative_to(ROOT).as_posix()}"


def build_payload(
    extraction_scaffold: dict[str, Any],
    raw_datum: dict[str, Any],
    exact_markov_modulus: dict[str, Any],
) -> dict[str, Any]:
    faithful_component = raw_datum["contract"]["must_emit"][1]
    return {
        "artifact": "oph_bw_fixed_local_collar_faithful_modular_defect_scaffold",
        "generated_utc": _timestamp(),
        "status": "minimal_modular_defect_extension",
        "public_promotion_allowed": False,
        "exact_missing_object": "fixed_local_collar_faithful_modular_defect_vanishing",
        "parent_missing_witness": extraction_scaffold["remaining_missing_emitted_witness"],
        "parent_extraction_object": extraction_scaffold["precise_missing_object_name"],
        "smaller_comparison_witness": exact_markov_modulus["exact_missing_object"],
        "smaller_comparison_witness_artifact": _artifact_ref(EXACT_MARKOV_MODULUS),
        "blocking_side_condition_artifact": _artifact_ref(COMMON_FLOOR),
        "role": (
            "Package the faithful modular-defect term that sits directly beneath the full carried-collar "
            "schedule and above the exact-Markov comparison witness on each fixed local collar model."
        ),
        "contract": {
            "for_fixed_models": raw_datum["contract"]["for_fixed_models"],
            "must_emit": "4 * lambda_{*,n,m,delta}^{-1} * delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0",
            "components": [
                exact_markov_modulus["contract"]["must_emit"],
                faithful_component,
            ],
            "meaning": (
                "On each fixed faithful collar model, the modular-additivity defect carried by the "
                "exact-Markov comparison state vanishes along the realized refinement chain."
            ),
        },
        "reduction_from_smaller_inputs": {
            "comparison_witness_artifact": _artifact_ref(EXACT_MARKOV_MODULUS),
            "faithfulness_side_condition": faithful_component,
            "faithfulness_side_condition_artifact": _artifact_ref(COMMON_FLOOR),
            "comparison_reference_floor_transfer": build_comparison_reference_floor_transfer(
                exact_markov_artifact=_artifact_ref(EXACT_MARKOV_MODULUS),
                spectral_floor_artifact=_artifact_ref(COMMON_FLOOR),
            ),
            "theorem": (
                "Pair the latent comparison-reference floor transfer with the modular-transport "
                "proposition: after the exact-Markov comparison marginals inherit the same eventual "
                "floor up to lambda_bar_{m,delta} / 2, the modular-additivity defect differs from the "
                "exact Markov value by at most 4 * lambda_*^{-1} * delta^M with "
                "lambda_* = lambda_bar_{m,delta} / 2 for all sufficiently large n."
            ),
            "eventual_common_floor": "lambda_* = lambda_bar_{m,delta} / 2 on all relevant marginals for large n",
            "status_on_fill": "faithful_modular_defect_closed",
        },
        "position_inside_carried_schedule": {
            "this_term": "4 * lambda_{*,n,m,delta}^{-1} * delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0",
            "other_term_still_needed": "r_FR(epsilon_{n,m,delta}) -> 0",
            "combined_parent_formula": extraction_scaffold["remaining_missing_emitted_witness_formula"],
        },
        "joint_schedule_term_frontier": build_schedule_term_frontier(
            constructive_recovery_artifact=_artifact_ref(CONSTRUCTIVE_RECOVERY),
            faithful_modular_defect_artifact=_artifact_ref(DEFAULT_OUT),
            carried_schedule_artifact=_artifact_ref(
                ROOT / "particles" / "runs" / "uv" / "bw_carried_collar_schedule_scaffold.json"
            ),
        ),
        "why_this_is_intermediate": [
            "This isolates the faithful modular-additivity burden from the full carried-collar schedule.",
            "It is smaller than the carried schedule because it omits the separate Fawzi-Renner recovery remainder.",
            "It is larger than the exact-Markov comparison witness because it still needs the eventual modular-transport common floor, which is now surfaced as its own artifact.",
        ],
        "notes": [
            "This scaffold does not claim the faithful modular-defect term is already emitted on the live corpus.",
            "It captures the second carried-error term singled out in the technical supplement, not the full carried schedule.",
            "The only nonlatent lower side condition on the live branch is the eventual modular-transport common floor recorded by the linked faithfulness-side artifact.",
            "No second spectral-floor artifact for the exact-Markov comparison family is missing: that common floor is inherited from the transported marginals once the exact-Markov modulus goes to zero on the fixed collar model.",
            "Closing this term together with the constructive recovery remainder closes the carried-collar schedule itself.",
            "The actual solver frontier above this witness is now recorded as the two-term pair, not as a separately targeted schedule object.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the fixed-local-collar faithful modular-defect scaffold."
    )
    parser.add_argument("--extraction-scaffold", default=str(EXTRACTION_SCAFFOLD))
    parser.add_argument("--raw-datum", default=str(RAW_DATUM))
    parser.add_argument("--exact-markov-modulus", default=str(EXACT_MARKOV_MODULUS))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(
        _load_json(Path(args.extraction_scaffold)),
        _load_json(Path(args.raw_datum)),
        _load_json(Path(args.exact_markov_modulus)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
