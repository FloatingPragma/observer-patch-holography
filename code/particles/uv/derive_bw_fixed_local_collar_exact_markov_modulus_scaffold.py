#!/usr/bin/env python3
"""Emit the fixed-local-collar exact-Markov modulus scaffold."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
RAW_DATUM = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_markov_faithfulness_datum.json"
FAITHFUL_MODULAR_DEFECT = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_faithful_modular_defect_scaffold.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_exact_markov_modulus_scaffold.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_ref(path: Path) -> str:
    return f"code/{path.relative_to(ROOT).as_posix()}"


def build_payload(raw_datum: dict[str, Any]) -> dict[str, Any]:
    cmi_component = raw_datum["contract"]["must_emit"][0]
    faithful_component = raw_datum["contract"]["must_emit"][1]
    return {
        "artifact": "oph_bw_fixed_local_collar_exact_markov_modulus_scaffold",
        "generated_utc": _timestamp(),
        "status": "minimal_local_comparison_extension",
        "public_promotion_allowed": False,
        "exact_missing_object": "fixed_local_collar_exact_markov_modulus_vanishing",
        "parent_raw_datum": raw_datum["exact_missing_object"],
        "parent_missing_witness": raw_datum["parent_missing_witness"],
        "parent_extraction_object": raw_datum["parent_extraction_object"],
        "role": (
            "Package the smaller exact-Markov comparison witness on each fixed local collar model: "
            "the exact-Markov distance modulus evaluated along the realized refinement chain tends to zero."
        ),
        "contract": {
            "for_fixed_models": raw_datum["contract"]["for_fixed_models"],
            "must_emit": "delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0",
            "definition": (
                "delta^M_{m,delta}(epsilon) = sup{ inf_{sigma in M_{m,delta}} ||rho - sigma||_1 : "
                "I(A_{m,delta}:D_{m,delta}|B_{m,delta})_rho <= epsilon }"
            ),
            "derived_from_component": cmi_component,
            "theorem_basis": (
                "On one fixed finite-dimensional collar model, compactness of state space and continuity "
                "of conditional mutual information force the exact-Markov distance modulus to vanish as "
                "epsilon -> 0."
            ),
            "must_not_assume": raw_datum["contract"]["must_not_assume"],
        },
        "feeds_follow_on_modular_defect": {
            "artifact": _artifact_ref(FAITHFUL_MODULAR_DEFECT),
            "formula": "4 * lambda_{*,n,m,delta}^{-1} * delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0",
            "still_needs_side_condition": faithful_component,
        },
        "already_packaged_below_this_witness": raw_datum["already_packaged_below_this_datum"],
        "why_this_is_smaller": [
            "This witness strips the fixed-local-collar datum down to the exact-Markov comparison convergence alone.",
            "It isolates the compactness-and-continuity step before any faithful spectral lower bound is used.",
            "The current raw datum remains larger because it bundles this exact-Markov convergence with eventual collarwise faithfulness.",
        ],
        "notes": [
            "This scaffold does not claim the exact-Markov modulus witness is already emitted on the live corpus.",
            "It is strictly smaller than the fixed-local-collar Markov/faithfulness datum because it removes the spectral-floor side condition.",
            "By itself it does not control the modular-additivity defect; that requires the faithful lower spectral bound as a separate side condition.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the fixed-local-collar exact-Markov modulus scaffold.")
    parser.add_argument("--raw-datum", default=str(RAW_DATUM))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(_load_json(Path(args.raw_datum)))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
