#!/usr/bin/env python3
"""Emit the carried-collar vanishing-schedule scaffold beneath cap-pair extraction."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bw_collar_honesty import (
    CARRIED_SCHEDULE_FORMULA,
    CONSTRUCTIVE_RECOVERY_FORMULA,
    build_local_honesty_gate,
    build_local_obligation_ledger,
    build_schedule_term_frontier,
)


ROOT = Path(__file__).resolve().parents[2]
EXTRACTION_SCAFFOLD = ROOT / "particles" / "runs" / "uv" / "bw_scaling_limit_cap_pair_extraction_scaffold.json"
RAW_DATUM = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_markov_faithfulness_datum.json"
CONSTRUCTIVE_RECOVERY = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_constructive_recovery_scaffold.json"
)
EXACT_MARKOV_MODULUS = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_exact_markov_modulus_scaffold.json"
FAITHFUL_MODULAR_DEFECT = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_faithful_modular_defect_scaffold.json"
)
COMMON_FLOOR = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_modular_transport_common_floor_scaffold.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "uv" / "bw_carried_collar_schedule_scaffold.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_ref(path: Path) -> str:
    return f"code/{path.relative_to(ROOT).as_posix()}"


def build_payload(
    extraction_scaffold: dict[str, Any],
    raw_datum: dict[str, Any],
    constructive_recovery: dict[str, Any],
    exact_markov_modulus: dict[str, Any],
    faithful_modular_defect: dict[str, Any],
) -> dict[str, Any]:
    return {
        "artifact": "oph_bw_carried_collar_schedule_scaffold",
        "generated_utc": _timestamp(),
        "status": "minimal_emitted_witness_extension",
        "public_promotion_allowed": False,
        "exact_missing_object": extraction_scaffold["remaining_missing_emitted_witness"],
        "theorem_object": extraction_scaffold["precise_missing_object_name"],
        "derived_missing_witness": extraction_scaffold["remaining_missing_emitted_witness"],
        "derived_missing_witness_artifact": _artifact_ref(DEFAULT_OUT),
        "raw_input_frontier": raw_datum["exact_missing_object"],
        "raw_input_frontier_artifact": _artifact_ref(RAW_DATUM),
        "intermediate_substep": exact_markov_modulus["exact_missing_object"],
        "parent_extraction_object": extraction_scaffold["precise_missing_object_name"],
        "smaller_raw_datum": raw_datum["exact_missing_object"],
        "smaller_raw_datum_artifact": _artifact_ref(RAW_DATUM),
        "single_live_missing_clause_artifact": _artifact_ref(COMMON_FLOOR),
        "role": (
            "Package the exact carried-collar vanishing schedule that sits directly above the fixed-local-collar "
            "Markov/faithfulness datum and directly below scaling-limit cap-pair extraction."
        ),
        "schedule_contract": {
            "formula": CARRIED_SCHEDULE_FORMULA,
            "components": extraction_scaffold["remaining_missing_emitted_witness_components"],
            "domain": extraction_scaffold["remaining_missing_emitted_witness_domain"],
            "meaning": (
                "On each fixed local collar model, the carried-collar remainder is explicitly bounded by one "
                "transport-fragility term and one faithfulness-weighted collar defect term, and that combined "
                "schedule vanishes along the realized refinement chain."
            ),
        },
        "reduction_from_raw_datum": {
            "reduction_theorem": extraction_scaffold["schedule_reduction_theorem"],
            "required_raw_components": raw_datum["contract"]["must_emit"],
            "constructive_recovery_witness_artifact": _artifact_ref(CONSTRUCTIVE_RECOVERY),
            "derived_exact_markov_comparison_witness_artifact": _artifact_ref(EXACT_MARKOV_MODULUS),
            "faithful_modular_defect_witness_artifact": _artifact_ref(FAITHFUL_MODULAR_DEFECT),
            "status_on_fill": "carried_collar_schedule_closed",
        },
        "schedule_term_witnesses": [
            {
                "id": constructive_recovery["exact_missing_object"],
                "artifact": _artifact_ref(CONSTRUCTIVE_RECOVERY),
                "role": "markov_side_recovery_term",
            },
            {
                "id": faithful_modular_defect["exact_missing_object"],
                "artifact": _artifact_ref(FAITHFUL_MODULAR_DEFECT),
                "role": "faithfulness_weighted_modular_term",
            },
        ],
        "termwise_closure_frontier": build_schedule_term_frontier(
            constructive_recovery_artifact=_artifact_ref(CONSTRUCTIVE_RECOVERY),
            faithful_modular_defect_artifact=_artifact_ref(FAITHFUL_MODULAR_DEFECT),
            carried_schedule_artifact=_artifact_ref(DEFAULT_OUT),
        ),
        "decomposed_error_terms": {
            "constructive_recovery_remainder": {
                "id": constructive_recovery["exact_missing_object"],
                "artifact": _artifact_ref(CONSTRUCTIVE_RECOVERY),
                "formula": constructive_recovery["contract"]["must_emit"],
                "source_component": constructive_recovery["contract"]["derived_from_component"],
                "meaning": constructive_recovery["role"],
            },
            "faithful_modular_defect_remainder": {
                "id": faithful_modular_defect["exact_missing_object"],
                "artifact": _artifact_ref(FAITHFUL_MODULAR_DEFECT),
                "formula": faithful_modular_defect["contract"]["must_emit"],
                "smaller_comparison_witness": exact_markov_modulus["exact_missing_object"],
            },
        },
        "obligation_ledger": build_local_obligation_ledger(
            constructive_recovery_artifact=_artifact_ref(CONSTRUCTIVE_RECOVERY),
            exact_markov_artifact=_artifact_ref(EXACT_MARKOV_MODULUS),
            faithful_modular_defect_artifact=_artifact_ref(FAITHFUL_MODULAR_DEFECT),
            carried_schedule_artifact=_artifact_ref(DEFAULT_OUT),
        ),
        "honesty_gate": build_local_honesty_gate(
            carried_schedule_artifact=_artifact_ref(DEFAULT_OUT),
            constructive_recovery_artifact=_artifact_ref(CONSTRUCTIVE_RECOVERY),
            exact_markov_artifact=_artifact_ref(EXACT_MARKOV_MODULUS),
            faithful_modular_defect_artifact=_artifact_ref(FAITHFUL_MODULAR_DEFECT),
        ),
        "why_this_witness_is_sharp": [
            "This witness is the exact emitted schedule the extraction theorem consumes; it is smaller than cap-pair extraction but larger than the raw collarwise datum.",
            "The fixed-local-collar Markov/faithfulness datum remains the sharpest lower object because it isolates the two raw controls before they are assembled into one vanishing schedule.",
            "The exact-Markov comparison witness and faithful modular-defect term now expose the two lower local substeps inside that reduction.",
            "The only nonlatent lower input beneath this derived schedule is the eventual modular-transport common floor feeding the faithful modular-defect side; the recovery/Markov side is already latent once epsilon_{n,m,delta} -> 0.",
            "Once this schedule is emitted, the first UV exact object left is cap-pair extraction itself.",
        ],
        "promotion_boundary": {
            "closed_below_this_witness": [
                "reference_cap_local_test_system",
                "projectively_compatible_transported_cap_marginal_family",
                "asymptotic_transport_equivalence_certificate",
                raw_datum["exact_missing_object"],
            ],
            "not_closed_here": [
                extraction_scaffold["precise_missing_object_name"],
                "ordered_null_cut_pair_rigidity",
            ],
        },
        "notes": [
            "This scaffold does not claim the schedule is already emitted on the live corpus.",
            "It records the exact witness contract consumed by scaling-limit cap-pair extraction so solver work no longer jumps directly from the raw datum to the extraction theorem.",
            "The carried-collar witness is theorem-generated from the two term witnesses recorded here, so it should no longer be treated as an independent primitive solver target.",
            "On the current branch, the only nonlatent lower side condition still outside the emitted chain is the eventual fixed-local-collar modular-transport common floor recorded by the linked faithfulness-side artifact.",
            "The honesty gate records exactly which lower witnesses are still insufficient on their own, preventing premature promotion from one-sided local control.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the BW carried-collar vanishing-schedule scaffold.")
    parser.add_argument("--extraction-scaffold", default=str(EXTRACTION_SCAFFOLD))
    parser.add_argument("--raw-datum", default=str(RAW_DATUM))
    parser.add_argument("--constructive-recovery", default=str(CONSTRUCTIVE_RECOVERY))
    parser.add_argument("--exact-markov-modulus", default=str(EXACT_MARKOV_MODULUS))
    parser.add_argument("--faithful-modular-defect", default=str(FAITHFUL_MODULAR_DEFECT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(
        _load_json(Path(args.extraction_scaffold)),
        _load_json(Path(args.raw_datum)),
        _load_json(Path(args.constructive_recovery)),
        _load_json(Path(args.exact_markov_modulus)),
        _load_json(Path(args.faithful_modular_defect)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
