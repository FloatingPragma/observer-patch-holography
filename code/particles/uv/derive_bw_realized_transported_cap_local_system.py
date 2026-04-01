#!/usr/bin/env python3
"""Emit the quotient/local-level transported cap-local UV system."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from bw_collar_honesty import (
    build_local_honesty_gate,
    build_local_obligation_ledger,
    build_schedule_term_frontier,
)


ROOT = Path(__file__).resolve().parents[2]
EXTRACTION_SCAFFOLD = ROOT / "particles" / "runs" / "uv" / "bw_scaling_limit_cap_pair_extraction_scaffold.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "uv" / "bw_realized_transported_cap_local_system.json"
RAW_DATUM = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_markov_faithfulness_datum.json"
CARRIED_SCHEDULE = ROOT / "particles" / "runs" / "uv" / "bw_carried_collar_schedule_scaffold.json"
CONSTRUCTIVE_RECOVERY = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_constructive_recovery_scaffold.json"
)
EXACT_MARKOV_MODULUS = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_exact_markov_modulus_scaffold.json"
FAITHFUL_MODULAR_DEFECT = (
    ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_faithful_modular_defect_scaffold.json"
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_ref(path: Path) -> str:
    return f"code/{path.relative_to(ROOT).as_posix()}"


def build_payload(extraction_scaffold: dict[str, object]) -> dict[str, object]:
    filled = list(extraction_scaffold["fills_contract_witnesses"])
    remaining = str(extraction_scaffold["remaining_missing_emitted_witness"])
    remaining_formula = str(extraction_scaffold.get("remaining_missing_emitted_witness_formula", ""))
    smaller_raw_datum = str(extraction_scaffold.get("smaller_remaining_raw_datum", ""))
    smaller_raw_components = list(extraction_scaffold.get("smaller_remaining_raw_datum_components", []))
    intermediate_witness_chain = list(extraction_scaffold.get("intermediate_witness_chain", []))
    term_frontier = build_schedule_term_frontier(
        constructive_recovery_artifact=_artifact_ref(CONSTRUCTIVE_RECOVERY),
        faithful_modular_defect_artifact=_artifact_ref(FAITHFUL_MODULAR_DEFECT),
        carried_schedule_artifact=_artifact_ref(CARRIED_SCHEDULE),
    )
    return {
        "artifact": "oph_realized_transported_cap_local_system",
        "generated_utc": _timestamp(),
        "status": "constructed_prelimit_system_two_lower_emitted_witnesses_still_missing",
        "public_promotion_allowed": False,
        "role": (
            "Package the quotient/local *-isomorphism-level transported cap-local system "
            "that the scaling-limit cap-pair extraction theorem uses as input."
        ),
        "fills_contract_witnesses": filled,
        "remaining_missing_emitted_witness": remaining,
        "source_scaffold_artifact": extraction_scaffold["artifact"],
        "reference_cap_local_test_system": {
            "status": "packaged",
            "kind": "countable_directed_reference_test_family",
            "construction": (
                "Fix one reference cut presentation/chart along the cap boundary and use "
                "rational matrix-coefficient cap-local subalgebras on supported local regions "
                "s ⊂ C as the directed test family T_{m,s} ⊂ A_m(s)."
            ),
            "equivalence_level": "quotient_or_local_star_isomorphism_class",
        },
        "projectively_compatible_transported_cap_marginal_family": {
            "status": "packaged",
            "construction": (
                "For n ≥ m, define ω_{n→m}^{(χ)}(X) = ω_{ℓ_n}(τ_{n→m}^{(χ)}(X)), where τ_{n→m}^{(χ)} "
                "is refinement inclusion plus admissible recharting into the fixed reference chart."
            ),
            "compatibility": [
                "isotony of the directed local family",
                "projective restriction compatibility",
                "channel composition along the realized refinement chain",
            ],
        },
        "asymptotic_transport_equivalence_certificate": {
            "status": "packaged",
            "equivalence_statement": (
                "Admissible transport choices differ only by the derived boundary recharting/gauge action, "
                "so they define one OPH-stable local *-isomorphism class on the quotient test system."
            ),
            "scope": "reference_chart_transport_class",
        },
        "remaining_missing_witness_contract": {
            "id": remaining,
            "artifact": _artifact_ref(CARRIED_SCHEDULE),
            "formula": remaining_formula,
            "for_fixed_models": "every fixed local collar model (m, delta)",
            "meaning": (
                "Carried-collar remainder schedules vanish along the realized refinement chain "
                "after transport to one fixed collar model."
            ),
        },
        "smaller_remaining_raw_datum": {
            "id": smaller_raw_datum,
            "artifact": _artifact_ref(RAW_DATUM),
            "components": smaller_raw_components,
            "role": "raw fixed-local-collar datum that implies the eta schedule once emitted",
        },
        "remaining_witness_decomposition": intermediate_witness_chain,
        "schedule_term_witnesses": [
            {
                "id": "constructive_recovery_remainder_vanishing",
                "artifact": _artifact_ref(CONSTRUCTIVE_RECOVERY),
                "role": "markov_side_recovery_term",
            },
            {
                "id": "fixed_local_collar_faithful_modular_defect_vanishing",
                "artifact": _artifact_ref(FAITHFUL_MODULAR_DEFECT),
                "role": "faithfulness_weighted_modular_term",
            },
        ],
        "actual_solver_missing_emitted_witnesses": term_frontier["missing_emitted_witnesses"],
        "derived_remaining_input_witness": term_frontier["derived_parent_witness"],
        "derived_remaining_input_witness_closure_theorem": term_frontier["closure_theorem"],
        "remaining_witness_obligation_ledger": build_local_obligation_ledger(
            constructive_recovery_artifact=_artifact_ref(CONSTRUCTIVE_RECOVERY),
            exact_markov_artifact=_artifact_ref(EXACT_MARKOV_MODULUS),
            faithful_modular_defect_artifact=_artifact_ref(FAITHFUL_MODULAR_DEFECT),
            carried_schedule_artifact=_artifact_ref(CARRIED_SCHEDULE),
        ),
        "remaining_witness_honesty_gate": build_local_honesty_gate(
            carried_schedule_artifact=_artifact_ref(CARRIED_SCHEDULE),
            constructive_recovery_artifact=_artifact_ref(CONSTRUCTIVE_RECOVERY),
            exact_markov_artifact=_artifact_ref(EXACT_MARKOV_MODULUS),
            faithful_modular_defect_artifact=_artifact_ref(FAITHFUL_MODULAR_DEFECT),
            include_prelimit_system_artifact=_artifact_ref(DEFAULT_OUT),
        ),
        "remaining_witness_term_frontier": term_frontier,
        "promotion_boundary": {
            "closed_here": [
                "reference_cap_local_test_system",
                "projectively_compatible_transported_cap_marginal_family",
                "asymptotic_transport_equivalence_certificate",
            ],
            "not_closed_here": [
                "canonical_scaling_cap_pair_realization_from_transported_cap_marginals",
                "ordered_null_cut_pair_rigidity",
            ],
        },
        "next_exact_object": {
            "id": "canonical_scaling_cap_pair_realization_from_transported_cap_marginals",
            "remaining_witness": remaining,
        },
        "notes": [
            "This artifact is constructive but not yet the realized scaling-limit cap pair.",
            "It packages the prelimit transported cap-local system at the quotient/local *-isomorphism level only.",
            "The carried-collar witness now comes with a finer lower local family: constructive recovery, exact-Markov comparison convergence, faithful modular-defect vanishing, then the full eta schedule.",
            "The actual emitted solver frontier is the two-term pair beneath the derived eta schedule, not the schedule viewed as a separate primitive target.",
            "The honesty gate makes explicit that this prelimit package is still insufficient on its own for cap-pair promotion.",
            "The remaining emitted witnesses for cap-pair promotion are the constructive-recovery and faithful modular-defect terms on fixed local collar models.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the realized transported cap-local system artifact.")
    parser.add_argument("--extraction-scaffold", default=str(EXTRACTION_SCAFFOLD))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(_load_json(Path(args.extraction_scaffold)))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
