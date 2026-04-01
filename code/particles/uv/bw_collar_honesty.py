#!/usr/bin/env python3
"""Shared honest-closure metadata for the UV/BW fixed-collar frontier."""

from __future__ import annotations


CMI_COMPONENT = "epsilon_{n,m,delta} = I(A_{m,delta}:D_{m,delta}|B_{m,delta})_{omega_{n->m}} -> 0"
FAITHFUL_COMPONENT = (
    "exists lambda_bar_{m,delta} > 0 and N_{m,delta} such that for all n >= N_{m,delta} "
    "and every X in Xi^{mod}_{m,delta}, rho_{n->m,X} >= lambda_bar_{m,delta} * 1"
)
EXACT_MARKOV_FORMULA = "delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0"
CONSTRUCTIVE_RECOVERY_FORMULA = "r_FR(epsilon_{n,m,delta}) -> 0"
FAITHFUL_MODULAR_DEFECT_FORMULA = (
    "4 * lambda_{*,n,m,delta}^{-1} * delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0"
)
CARRIED_SCHEDULE_FORMULA = (
    "eta_{n,m,delta} = r_FR(epsilon_{n,m,delta}) + "
    "4 * lambda_{*,n,m,delta}^{-1} * delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0"
)

EXACT_MARKOV_WITNESS_ID = "fixed_local_collar_exact_markov_modulus_vanishing"
CONSTRUCTIVE_RECOVERY_ID = "constructive_recovery_remainder_vanishing"
FAITHFUL_MODULAR_DEFECT_ID = "fixed_local_collar_faithful_modular_defect_vanishing"
CARRIED_SCHEDULE_ID = "vanishing_carried_collar_schedule_on_fixed_local_collars"
COMPARISON_REFERENCE_FLOOR_TRANSFER_ID = "exact_markov_reference_eventual_common_floor_transfer"


def build_comparison_reference_floor_transfer(
    *,
    exact_markov_artifact: str,
    spectral_floor_artifact: str,
) -> dict[str, object]:
    return {
        "id": COMPARISON_REFERENCE_FLOOR_TRANSFER_ID,
        "requires_exact_markov_artifact": exact_markov_artifact,
        "requires_spectral_floor_artifact": spectral_floor_artifact,
        "statement": (
            "On one fixed local collar model, if the modular-transport marginals satisfy "
            "rho_{n->m,X} >= lambda_bar_{m,delta} * 1 eventually for every X in Xi^{mod}_{m,delta} "
            "and the exact-Markov replacements satisfy "
            "||rho_X - sigma_{n,m,delta,X}||_1 <= delta^M_{m,delta}(epsilon_{n,m,delta}) -> 0, then "
            "for all sufficiently large n the comparison marginals also satisfy "
            "sigma_{n,m,delta,X} >= (lambda_bar_{m,delta} / 2) * 1."
        ),
        "why_this_is_enough": (
            "No separate spectral-floor hypothesis on the exact-Markov comparison family is needed "
            "before applying the modular-transport estimate."
        ),
        "status_on_fill": "exact_markov_reference_common_floor_closed",
    }


def build_local_obligation_ledger(
    *,
    constructive_recovery_artifact: str,
    exact_markov_artifact: str,
    faithful_modular_defect_artifact: str,
    carried_schedule_artifact: str,
) -> list[dict[str, object]]:
    return [
        {
            "id": "collarwise_markov_input",
            "kind": "raw_input",
            "formula": CMI_COMPONENT,
            "status": "open_until_emitted_on_every_fixed_model",
            "unlocks": [
                EXACT_MARKOV_WITNESS_ID,
                CONSTRUCTIVE_RECOVERY_ID,
            ],
        },
        {
            "id": "collarwise_faithfulness_input",
            "kind": "raw_input",
            "formula": FAITHFUL_COMPONENT,
            "status": "open_until_emitted_on_every_fixed_model",
            "unlocks": [FAITHFUL_MODULAR_DEFECT_ID],
        },
        {
            "id": EXACT_MARKOV_WITNESS_ID,
            "kind": "derived_local_witness",
            "formula": EXACT_MARKOV_FORMULA,
            "artifact": exact_markov_artifact,
            "depends_on": ["collarwise_markov_input"],
            "status": "open_until_emitted",
        },
        {
            "id": CONSTRUCTIVE_RECOVERY_ID,
            "kind": "derived_error_term",
            "formula": CONSTRUCTIVE_RECOVERY_FORMULA,
            "artifact": constructive_recovery_artifact,
            "depends_on": ["collarwise_markov_input"],
            "status": "open_until_emitted",
        },
        {
            "id": FAITHFUL_MODULAR_DEFECT_ID,
            "kind": "derived_error_term",
            "formula": FAITHFUL_MODULAR_DEFECT_FORMULA,
            "artifact": faithful_modular_defect_artifact,
            "depends_on": [
                EXACT_MARKOV_WITNESS_ID,
                "collarwise_faithfulness_input",
            ],
            "status": "open_until_emitted",
        },
        {
            "id": CARRIED_SCHEDULE_ID,
            "kind": "emitted_parent_witness",
            "formula": CARRIED_SCHEDULE_FORMULA,
            "artifact": carried_schedule_artifact,
            "depends_on": [
                CONSTRUCTIVE_RECOVERY_ID,
                FAITHFUL_MODULAR_DEFECT_ID,
            ],
            "status": "open_until_emitted",
        },
    ]


def build_local_honesty_gate(
    *,
    carried_schedule_artifact: str,
    constructive_recovery_artifact: str,
    exact_markov_artifact: str,
    faithful_modular_defect_artifact: str,
    include_prelimit_system_artifact: str | None = None,
) -> dict[str, object]:
    insufficient = [
        {
            "id": CONSTRUCTIVE_RECOVERY_ID,
            "artifact": constructive_recovery_artifact,
            "why_not": (
                "This removes only the Fawzi-Renner recovery term; it does not control the separate "
                "faithfulness-weighted modular defect remainder."
            ),
        },
        {
            "id": EXACT_MARKOV_WITNESS_ID,
            "artifact": exact_markov_artifact,
            "why_not": (
                "This only controls distance to the exact-Markov set; it does not by itself bound "
                "the spectral-weighted modular defect term."
            ),
        },
        {
            "id": FAITHFUL_MODULAR_DEFECT_ID,
            "artifact": faithful_modular_defect_artifact,
            "why_not": (
                "This still omits the separate Fawzi-Renner recovery remainder, so the full eta "
                "schedule is not yet emitted."
            ),
        },
    ]
    if include_prelimit_system_artifact:
        insufficient.append(
            {
                "id": "realized_transported_cap_local_system",
                "artifact": include_prelimit_system_artifact,
                "why_not": (
                    "The prelimit transport package fixes the local system but does not itself emit "
                    "the vanishing carried-collar schedule consumed by cap-pair extraction."
                ),
            }
        )
    return {
        "status": "open",
        "promotion_rule": (
            "No theorem promotion is honest until the carried-collar schedule itself is emitted on "
            "every fixed local collar model."
        ),
        "required_raw_inputs": [
            CMI_COMPONENT,
            FAITHFUL_COMPONENT,
        ],
        "derived_witness_chain": [
            CONSTRUCTIVE_RECOVERY_ID,
            EXACT_MARKOV_WITNESS_ID,
            FAITHFUL_MODULAR_DEFECT_ID,
            CARRIED_SCHEDULE_ID,
        ],
        "insufficient_on_their_own": insufficient,
        "closure_artifact": carried_schedule_artifact,
    }


def build_schedule_term_frontier(
    *,
    constructive_recovery_artifact: str,
    faithful_modular_defect_artifact: str,
    carried_schedule_artifact: str,
) -> dict[str, object]:
    return {
        "status": "open_two_term_emitted_frontier",
        "missing_emitted_witnesses": [
            {
                "id": CONSTRUCTIVE_RECOVERY_ID,
                "artifact": constructive_recovery_artifact,
                "formula": CONSTRUCTIVE_RECOVERY_FORMULA,
                "role": "markov_side_recovery_term",
            },
            {
                "id": FAITHFUL_MODULAR_DEFECT_ID,
                "artifact": faithful_modular_defect_artifact,
                "formula": FAITHFUL_MODULAR_DEFECT_FORMULA,
                "role": "faithfulness_weighted_modular_term",
            },
        ],
        "derived_parent_witness": {
            "id": CARRIED_SCHEDULE_ID,
            "artifact": carried_schedule_artifact,
            "formula": CARRIED_SCHEDULE_FORMULA,
        },
        "closure_theorem": (
            "If the constructive recovery remainder and the faithful modular-defect term both vanish "
            "on every fixed local collar model, then the carried-collar schedule follows by termwise "
            "addition of the two nonnegative remainders."
        ),
        "status_on_fill": "carried_collar_schedule_closed",
        "why_this_sharpens_the_frontier": (
            "The carried-collar schedule is no longer treated as an independent primitive solver target; "
            "the actual emitted frontier is the two-term pair that generates it."
        ),
    }
