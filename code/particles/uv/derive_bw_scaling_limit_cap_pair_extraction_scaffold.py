#!/usr/bin/env python3
"""Emit the canonical scaling-limit cap-pair extraction scaffold."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from bw_collar_honesty import (
    CARRIED_SCHEDULE_FORMULA,
    CMI_COMPONENT,
    FAITHFUL_COMPONENT,
    build_comparison_reference_floor_transfer,
    build_local_honesty_gate,
    build_local_obligation_ledger,
    build_schedule_term_frontier,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "uv" / "bw_scaling_limit_cap_pair_extraction_scaffold.json"
PRELIMIT_SYSTEM = ROOT / "particles" / "runs" / "uv" / "bw_realized_transported_cap_local_system.json"
RAW_DATUM = ROOT / "particles" / "runs" / "uv" / "bw_fixed_local_collar_markov_faithfulness_datum.json"
CARRIED_SCHEDULE = ROOT / "particles" / "runs" / "uv" / "bw_carried_collar_schedule_scaffold.json"
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


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _artifact_ref(path: Path) -> str:
    return f"code/{path.relative_to(ROOT).as_posix()}"


def build_payload() -> dict[str, object]:
    filled_witnesses = [
        "reference_cap_local_test_system",
        "projectively_compatible_transported_cap_marginal_family",
        "asymptotic_transport_equivalence_certificate",
    ]
    remaining_witness = "vanishing_carried_collar_schedule_on_fixed_local_collars"
    smaller_raw_datum = "fixed_local_collar_markov_faithfulness_datum"
    remaining_formula = CARRIED_SCHEDULE_FORMULA
    schedule_term_frontier = build_schedule_term_frontier(
        constructive_recovery_artifact=_artifact_ref(CONSTRUCTIVE_RECOVERY),
        faithful_modular_defect_artifact=_artifact_ref(FAITHFUL_MODULAR_DEFECT),
        carried_schedule_artifact=_artifact_ref(CARRIED_SCHEDULE),
    )
    return {
        "artifact": "oph_bw_scaling_limit_cap_pair_extraction_scaffold",
        "generated_utc": _timestamp(),
        "status": "constructive_prelimit_system_two_lower_emitted_witnesses_still_missing",
        "public_promotion_allowed": False,
        "exact_missing_object": "scaling_limit_cap_pair_extraction",
        "precise_missing_object_name": "canonical_scaling_cap_pair_realization_from_transported_cap_marginals",
        "theorem_contract_name": "conditional_scaling_limit_cap_pair_extraction_theorem",
        "goal": (
            "Extract a canonical scaling-limit cap pair (A_infty(C), omega_infty^C) "
            "from transported cap marginals along the realized refinement-stable branch."
        ),
        "input_contract": {
            "must_use": [
                "transported_cap_marginals along a realized refinement chain",
                "inherited strip data with oriented null generator and ordered cut pair",
                "support-map control on the quasi-local cap net",
            ],
            "must_not_assume": [
                "type-I survival in the scaling limit",
                "a pre-existing BW automorphism theorem",
                "a fixed-cutoff density-matrix identification of the limit",
            ],
        },
        "fills_contract_witnesses": filled_witnesses,
        "remaining_missing_emitted_witness": remaining_witness,
        "remaining_missing_emitted_witness_artifact": _artifact_ref(CARRIED_SCHEDULE),
        "remaining_missing_emitted_witness_formula": remaining_formula,
        "remaining_missing_emitted_witness_components": [
            "r_FR(epsilon_{n,m,delta})",
            "lambda_{*,n,m,delta}",
            "delta^M_{m,delta}(epsilon_{n,m,delta})",
        ],
        "remaining_missing_emitted_witness_domain": "every fixed local collar model",
        "smaller_remaining_raw_datum": smaller_raw_datum,
        "smaller_remaining_raw_datum_artifact": _artifact_ref(RAW_DATUM),
        "smaller_remaining_raw_datum_components": [
            CMI_COMPONENT,
            FAITHFUL_COMPONENT,
        ],
        "single_live_missing_clause_artifact": _artifact_ref(COMMON_FLOOR),
        "single_live_missing_clause_closure_lemma": build_comparison_reference_floor_transfer(
            exact_markov_artifact=_artifact_ref(EXACT_MARKOV_MODULUS),
            spectral_floor_artifact=_artifact_ref(COMMON_FLOOR),
        ),
        "markov_side_status": "latent_from_epsilon_to_zero",
        "faithfulness_side_status": "open",
        "intermediate_witness_chain": [
            {
                "id": "constructive_recovery_remainder_vanishing",
                "artifact": _artifact_ref(CONSTRUCTIVE_RECOVERY),
            },
            {
                "id": "fixed_local_collar_exact_markov_modulus_vanishing",
                "artifact": _artifact_ref(EXACT_MARKOV_MODULUS),
            },
            {
                "id": "fixed_local_collar_faithful_modular_defect_vanishing",
                "artifact": _artifact_ref(FAITHFUL_MODULAR_DEFECT),
            },
            {
                "id": remaining_witness,
                "artifact": _artifact_ref(CARRIED_SCHEDULE),
            },
        ],
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
        "actual_solver_missing_emitted_witnesses": schedule_term_frontier["missing_emitted_witnesses"],
        "derived_remaining_input_witness": schedule_term_frontier["derived_parent_witness"],
        "derived_remaining_input_witness_closure_theorem": schedule_term_frontier["closure_theorem"],
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
            include_prelimit_system_artifact=_artifact_ref(PRELIMIT_SYSTEM),
        ),
        "remaining_witness_term_frontier": schedule_term_frontier,
        "schedule_reduction_theorem": (
            "For each fixed local collar model, the eta_{n,m,delta} vanishing schedule follows once the collarwise conditional mutual information vanishes and the faithful lower spectral bound is eventually uniform."
        ),
        "missing_input_witnesses": [
            "constructive_recovery_remainder_vanishing",
            "fixed_local_collar_faithful_modular_defect_vanishing",
        ],
        "prelimit_realized_system_artifact": _artifact_ref(PRELIMIT_SYSTEM),
        "theorem_assumptions": [
            "countable directed cap-local test system in one fixed reference chart",
            "transported pullback states on each fixed local cap algebra",
            "projective restriction compatibility across the directed local system",
            "asymptotic transport coherence across admissible transport choices",
            "vanishing carried collar schedules on every fixed local collar model",
        ],
        "checks": [
            "local_weak_star_extraction",
            "diagonal_subnet_extraction_on_countable_local_test_family",
            "projective_gluing_of_limit_functionals",
            "GNS_completion_of_the_local_quotient_family",
            "vanishing_carried_collar_errors",
            "refinement_stable_transport_compatibility",
            "uniqueness_of_extracted_limit_action",
        ],
        "output_certificate": {
            "extracted_pair": "(A_infty(C), omega_infty^C)",
            "canonicality": "unique local *-isomorphism class inside the asymptotic transport-equivalence class",
            "typeI_required": False,
            "residual_modular_class": "q_BW(C) remains open until ordered-cut rigidity closes",
            "status_on_fill": "realized_limit_pair_closed",
        },
        "collar_corollary": (
            "Exact fixed-cutoff collar splice and additivity identities pass to the extracted pair whenever they are used only through carried-collar remainders that vanish on fixed local collar models."
        ),
        "follow_on_object": {
            "id": "ordered_null_cut_pair_rigidity",
            "role": (
                "Collapse the residual cap-preserving conformal freedom on the realized "
                "limit pair to the unique BW hyperbolic subgroup."
            ),
        },
        "final_target": "sigma_t^{omega_infty^C} = alpha_{lambda_C(2 pi t)}",
        "notes": [
            "This is the first half of the UV/BW internalization route.",
            "The local quotient/prelimit transported cap-local system can already be packaged from the current corpus at the reference-chart and asymptotic transport-equivalence level.",
            "The compactness/extraction theorem itself is not the missing proof; the current remaining blocker is the two-term carried-collar emitted frontier beneath the derived eta schedule.",
            "The carried schedule itself is now recorded as theorem-generated from those two term witnesses: r_FR(epsilon_{n,m,delta}) -> 0 and 4 * lambda_{*,n,m,delta}^{-1} * delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0.",
            "That combined schedule reduces one level lower to a fixed-local-collar Markov/faithfulness datum: collarwise CMI vanishing plus an eventual common floor on the finite modular-transport family.",
            "Inside that raw datum, the constructive recovery witness, the exact-Markov comparison modulus, and the faithful modular-defect term are now split as separate lower local scaffolds.",
            "No separate exact-Markov-reference spectral input remains hidden below the faithful modular-defect term: on one fixed collar model, the eventual floor for the transported marginals transfers to the exact-Markov comparison marginals once delta^M -> 0.",
            "On the current branch the only nonlatent lower input still outside the emitted chain is the eventual fixed-local-collar modular-transport common floor feeding the faithful modular-defect term; the recovery/Markov side is already latent once epsilon_{n,m,delta} -> 0.",
            "The artifact now also carries a machine-readable honesty gate so prelimit transport packaging cannot be mistaken for the missing emitted collar schedule.",
            "It does not by itself prove the BW automorphism law.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the BW scaling-limit cap-pair extraction scaffold.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload()
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
