#!/usr/bin/env python3
"""Deterministic candidate generator for the issue-647 mining campaign.

The generator runs only under the enablement policy: candidate generation
and evaluation enabled, public-data access and target registration sealed.
It reads the frozen registries, executes the implemented slot producers,
scores every candidate with the frozen integer weights, orders them with the
frozen tie breaker, and writes ``outputs/candidate_registry.json``. It reads
no public measurement, no comparison payload, and no target value, and it
uses no randomness and no clock.

The first implemented producer is ``z6_charge_line_congruences``. Its
candidates are exact discrete relations of the pinned diagonal
:math:`\\mathbb{Z}_6` descent: the representation descent congruence, its
color-singlet and color-triplet charge corollaries, the six-fold flux-class
count, and the fractional-charge kill rule. Every relation is certified by
finite exhaustion over the residue classes of the kernel action recorded in
the registry itself, so the independent verifier can recompute the entire
content. Slots without an implemented producer are recorded as awaiting
generation, and scoring stays sealed until the registry is complete.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

PACKAGE_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = PACKAGE_ROOT.parent.parent
OUTPUT_PATH = PACKAGE_ROOT / "outputs" / "candidate_registry.json"

SCHEMA = "oph.invariant_mining.candidate_registry.v1"
STATUS = "GENERATED_PARTIAL__SCORING_SEALED"
IMPLEMENTED_SLOTS = ("z6_charge_line_congruences",)


class GenerationError(ValueError):
    """The candidate generator refused to run or emit."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise GenerationError(message)


def canonical_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    ).encode("utf-8")


def sha256_bytes(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"JSON document must be an object: {path}")
    return value


# ---------------------------------------------------------------------------
# Z6 descent arithmetic: exact finite exhaustion
# ---------------------------------------------------------------------------


def z6_kernel_phase_class(triality: int, duality: int, six_y: int) -> int:
    """Exact class of the pinned kernel action on a representation.

    The diagonal kernel generator acts on a representation with
    :math:`\\mathrm{SU}(3)` triality ``t``, :math:`\\mathrm{SU}(2)` duality
    ``d``, and hypercharge ``Y`` by the phase
    :math:`\\exp(2\\pi i\\,(t/3 + d/2 + Y))`. With ``six_y = 6Y`` the phase
    class in :math:`\\mathbb{Z}_6` is ``(2t + 3d + six_y) mod 6``; the
    representation descends to the quotient exactly when the class is zero.
    """

    return (2 * (triality % 3) + 3 * (duality % 2) + six_y) % 6


def descent_congruence_table() -> list[dict[str, int | bool]]:
    """Exhaustive descent table over all residue classes of the kernel action."""

    table = []
    for triality in range(3):
        for duality in range(2):
            for six_y in range(6):
                phase = z6_kernel_phase_class(triality, duality, six_y)
                table.append(
                    {
                        "triality": triality,
                        "duality": duality,
                        "six_y_mod_6": six_y,
                        "kernel_phase_class": phase,
                        "descends": phase == 0,
                    }
                )
    return table


def _descent_corollaries() -> dict[str, Any]:
    table = descent_congruence_table()
    color_singlet_weak_singlet = [
        row for row in table if row["triality"] == 0 and row["duality"] == 0
    ]
    singlet_descending = [
        row["six_y_mod_6"] for row in color_singlet_weak_singlet if row["descends"]
    ]
    color_triplet = [row for row in table if row["triality"] == 1]
    triplet_descending = sorted(
        {
            (row["duality"], row["six_y_mod_6"])
            for row in color_triplet
            if row["descends"]
        }
    )
    return {
        "descending_rows": sum(1 for row in table if row["descends"]),
        "total_rows": len(table),
        "color_singlet_weak_singlet_descending_six_y": singlet_descending,
        "color_triplet_descending_duality_six_y": [
            list(pair) for pair in triplet_descending
        ],
    }


# ---------------------------------------------------------------------------
# Z6 slot producer
# ---------------------------------------------------------------------------


def produce_z6_candidates(
    slots: dict[str, Any],
    grammar: dict[str, Any],
    nuisances: dict[str, Any],
    exposure: dict[str, Any],
) -> list[dict[str, Any]]:
    slot = next(
        row
        for row in slots["slots"]
        if row["slot_id"] == "z6_charge_line_congruences"
    )
    allowed_classes = set(slot["allowed_grammar_class_ids"])
    corollaries = _descent_corollaries()
    require(
        corollaries["descending_rows"] == 6
        and corollaries["color_singlet_weak_singlet_descending_six_y"] == [0]
        and corollaries["color_triplet_descending_duality_six_y"]
        == [[0, 4], [1, 1]],
        "descent exhaustion drift",
    )

    baseline_contract = {
        "baseline_class": "minimal Standard Model plus general relativity",
        "global_form_policy": (
            "the gauge group global form is a free discrete choice among "
            "SU(3)xSU(2)xU(1) and its Z2, Z3, and Z6 quotients; every "
            "Standard Model matter field satisfies the descent congruence, "
            "so the choice is unconstrained by the established matter "
            "content"
        ),
        "field_content": "established three-generation matter, no extra light states assumed",
        "counterexample": (
            "a color-singlet weak-singlet representation with 6Y = 3 is a "
            "consistent representation of the unquotiented global form and "
            "violates the descent congruence, so the baseline image is not "
            "contained in the relation"
        ),
        "baseline_freedom_class": "PERMITS_AS_FREE_DISCRETE_CHOICE",
    }
    shared_nuisance_matrix = {
        "basis_and_coordinate_choice": {
            "status": "QUOTIENTED",
            "certificate": (
                "the relation is stated on character data (triality, "
                "duality, hypercharge class), which is invariant under "
                "basis and coordinate changes"
            ),
        },
        "physical_sector_choice": {
            "status": "NAMED_CONDITIONAL_BRANCH",
            "branch": "z6-global-form-source-selection",
            "note": (
                "the diagonal Z6 selection is source-derived on the "
                "declared carrier; its physical attachment is open, so "
                "every candidate carries the named branch label"
            ),
        },
        "source_admissible_completion": {
            "status": "NAMED_COMPLETION_PREMISE",
            "premise": "completions preserving the pinned diagonal kernel",
        },
        "detector_transfer_calibration": {
            "status": "DEFERRED_UNTIL_PHYSICALIZATION",
        },
        "common_dimensionful_scale": {
            "status": "NOT_APPLICABLE_DIMENSIONLESS_DISCRETE",
        },
    }
    exposure_surfaces = sorted(
        row["exposure_id"]
        for row in exposure["surfaces"]
        if "z6_charge_line_congruences" in row["applies_to_slot_ids"]
    )

    def expression(kind: str, terminals: list[str], operator_units: int) -> dict[str, Any]:
        budget = grammar["complexity_budget"]
        depth = 3
        nodes = 1 + 1 + len(terminals)
        require(depth <= budget["maximum_ast_depth"], "depth budget exceeded")
        require(nodes <= budget["maximum_ast_nodes"], "node budget exceeded")
        require(
            len(set(terminals)) <= budget["maximum_distinct_registered_terminals"],
            "terminal budget exceeded",
        )
        return {
            "form": kind,
            "registered_terminals": terminals,
            "ast_depth": depth,
            "ast_nodes": nodes,
            "complexity_units": operator_units,
        }

    candidates = [
        {
            "candidate_id": "z6-descent-congruence",
            "slot_id": slot["slot_id"],
            "grammar_class": "discrete_congruence_parity_ordering_multiplicity",
            "statement": (
                "every physical representation satisfies the descent "
                "congruence 2t + 3d + 6Y = 0 mod 6 through the pinned "
                "diagonal kernel"
            ),
            "expression": expression(
                "declared_modulus_class[6] of inner_product(integer "
                "coefficients, character vector)",
                ["registered_integer", "registered_character"],
                1 * 2 + 3 + 4,
            ),
            "relation_certificate": {
                "kind": "finite_exhaustion_over_kernel_residue_classes",
                "rows": 36,
                "descending_rows": 6,
                "recomputable_from": "descent_congruence_table",
            },
            "kill_rule": (
                "one confirmed physical state violating the congruence "
                "falsifies the selected global form branch"
            ),
        },
        {
            "candidate_id": "z6-color-singlet-integer-charge",
            "slot_id": slot["slot_id"],
            "grammar_class": "exact_zero_forbidden_transition",
            "statement": (
                "color-singlet weak-singlet states carry integer electric "
                "charge; fractionally charged color singlets are forbidden"
            ),
            "expression": expression(
                "exact_zero of declared_modulus_class[6] restricted to "
                "trivial triality and duality",
                ["registered_character", "registered_integer"],
                1 * 2 + 3 + 4,
            ),
            "relation_certificate": {
                "kind": "finite_exhaustion_over_kernel_residue_classes",
                "restriction": "triality 0, duality 0",
                "descending_six_y_classes": [0],
                "recomputable_from": "descent_congruence_table",
            },
            "kill_rule": (
                "one confirmed free color-singlet particle with fractional "
                "electric charge falsifies the selected global form branch"
            ),
        },
        {
            "candidate_id": "z6-color-triplet-charge-classes",
            "slot_id": slot["slot_id"],
            "grammar_class": "discrete_congruence_parity_ordering_multiplicity",
            "statement": (
                "color-triplet states carry electric charge in the "
                "one-third classes fixed by the congruence: weak singlets "
                "sit at 6Y = 4 mod 6 and weak doublets at 6Y = 1 mod 6"
            ),
            "expression": expression(
                "declared_modulus_class[6] restricted to unit triality",
                ["registered_character", "registered_integer"],
                1 * 2 + 3 + 4,
            ),
            "relation_certificate": {
                "kind": "finite_exhaustion_over_kernel_residue_classes",
                "restriction": "triality 1",
                "descending_duality_six_y": [[0, 4], [1, 1]],
                "recomputable_from": "descent_congruence_table",
            },
            "kill_rule": (
                "one confirmed color-triplet state outside the descending "
                "classes falsifies the selected global form branch"
            ),
        },
        {
            "candidate_id": "z6-flux-class-count",
            "slot_id": slot["slot_id"],
            "grammar_class": "trace_determinant_character_index",
            "statement": (
                "the quotient fundamental group is cyclic of order six, so "
                "magnetic flux classes take exactly six values"
            ),
            "expression": expression(
                "index of the pinned quotient kernel",
                ["registered_index"],
                1 + 2,
            ),
            "relation_certificate": {
                "kind": "pinned_kernel_order",
                "kernel_order": 6,
                "source": "exact diagonal Z6 quotient",
            },
            "kill_rule": (
                "a physical flux classification with a different finite "
                "order falsifies the selected global form branch"
            ),
        },
        {
            "candidate_id": "z6-line-lattice-exclusion",
            "slot_id": slot["slot_id"],
            "grammar_class": "exact_zero_forbidden_transition",
            "statement": (
                "electric line classes outside the descended "
                "representation lattice are forbidden in the quotient "
                "branch"
            ),
            "expression": expression(
                "exact_zero of non-descending character classes",
                ["registered_character", "registered_index"],
                1 * 2 + 3 + 4,
            ),
            "relation_certificate": {
                "kind": "finite_exhaustion_over_kernel_residue_classes",
                "forbidden_rows": 30,
                "recomputable_from": "descent_congruence_table",
            },
            "kill_rule": (
                "one confirmed electric line class outside the descended "
                "lattice falsifies the selected global form branch"
            ),
        },
    ]

    for candidate in candidates:
        require(
            candidate["grammar_class"] in allowed_classes,
            f"grammar class not allowed for slot: {candidate['candidate_id']}",
        )
        candidate["source_ancestry"] = {
            "feature_ids": list(slot["required_feature_ids"]),
            "target_ancestry": "SOURCE_ONLY_NO_PUBLIC_COMPARISON_INPUT",
        }
        candidate["nuisance_matrix"] = shared_nuisance_matrix
        candidate["baseline_contract"] = baseline_contract
        candidate["physicalization_status"] = (
            "OPEN_MAP__CONTINUUM_GAUGE_AND_LABORATORY_ATTACHMENT_REQUIRED"
        )
        candidate["exposure_surfaces"] = exposure_surfaces
        candidate["weight_claims"] = {
            "exact_global_certificate": True,
            "independent_recomputation": True,
            "physicalization_complete": False,
            "baseline_freedom_counterexample": True,
            "all_registered_completions_covered": False,
            "all_continuous_parameters_covered": False,
            "conditional_branch": True,
            "open_physical_map": True,
        }
    return candidates


# ---------------------------------------------------------------------------
# Scoring with the frozen weights
# ---------------------------------------------------------------------------


def score_candidate(candidate: dict[str, Any], ranking: dict[str, Any]) -> int:
    weights = ranking["weights"]
    claims = candidate["weight_claims"]
    score = 0
    if claims["exact_global_certificate"]:
        score += weights["exact_global_certificate"]
    if claims["independent_recomputation"]:
        score += weights["independent_recomputation"]
    if claims["physicalization_complete"]:
        score += weights["physicalization_complete"]
    if claims["baseline_freedom_counterexample"]:
        score += weights["baseline_freedom_counterexample"]
    if claims["all_registered_completions_covered"]:
        score += weights["all_registered_completions_covered"]
    if claims["all_continuous_parameters_covered"]:
        score += weights["all_continuous_parameters_covered"]
    if claims["conditional_branch"]:
        score += weights["conditional_branch_penalty"]
    if claims["open_physical_map"]:
        score += weights["open_physical_map_penalty"]
    score += (
        weights["expression_complexity_unit_penalty"]
        * candidate["expression"]["complexity_units"]
    )
    return score


def rank_candidates(
    candidates: list[dict[str, Any]], ranking: dict[str, Any]
) -> list[dict[str, Any]]:
    for candidate in candidates:
        candidate["score"] = score_candidate(candidate, ranking)
    ordered = sorted(
        candidates,
        key=lambda candidate: (
            -candidate["score"],
            candidate["expression"]["complexity_units"],
            candidate["expression"]["form"],
            candidate["candidate_id"],
        ),
    )
    for position, candidate in enumerate(ordered, start=1):
        candidate["rank"] = position
    return ordered


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_registry() -> dict[str, Any]:
    policy = load_json(PACKAGE_ROOT / "policy" / "pregeneration_policy.json")
    require(
        policy.get("candidate_generator_enabled") is True
        and policy.get("candidate_evaluator_enabled") is True,
        "candidate generation is not enabled",
    )
    require(
        policy.get("public_data_access_enabled") is False
        and policy.get("target_payloads_registered") is False,
        "comparison surfaces must stay sealed during generation",
    )
    slots = load_json(PACKAGE_ROOT / "data" / "producer_slots.json")
    grammar = load_json(PACKAGE_ROOT / "data" / "observable_grammar.json")
    nuisances = load_json(PACKAGE_ROOT / "data" / "nuisance_registry.json")
    ranking = load_json(PACKAGE_ROOT / "data" / "ranking_policy.json")
    exposure = load_json(PACKAGE_ROOT / "data" / "exposed_data_registry.json")

    candidates = rank_candidates(
        produce_z6_candidates(slots, grammar, nuisances, exposure), ranking
    )
    budget = grammar["complexity_budget"]
    require(
        len(candidates) <= budget["maximum_generated_candidates"],
        "candidate budget exceeded",
    )
    forbidden_terminals = set(grammar["forbidden_terminals"])
    for candidate in candidates:
        require(
            not (
                set(candidate["expression"]["registered_terminals"])
                & forbidden_terminals
            ),
            f"forbidden terminal in {candidate['candidate_id']}",
        )

    direct_slot = policy["direct_n_contract"]["slot_id"]
    slot_states = {}
    for row in slots["slots"]:
        slot_id = row["slot_id"]
        if slot_id == direct_slot:
            slot_states[slot_id] = "EXCLUDED_EXTERNALLY_ORDERED_DIRECT_N"
        elif slot_id in IMPLEMENTED_SLOTS:
            slot_states[slot_id] = "GENERATED"
        else:
            slot_states[slot_id] = "REGISTERED_AWAITING_GENERATION"

    registry = {
        "schema": SCHEMA,
        "issue": 647,
        "status": STATUS,
        "generation_inputs": {
            "policy_sha256": sha256_bytes(
                (PACKAGE_ROOT / "policy" / "pregeneration_policy.json").read_bytes()
            ),
            "producer_slots_sha256": sha256_bytes(
                (PACKAGE_ROOT / "data" / "producer_slots.json").read_bytes()
            ),
            "observable_grammar_sha256": sha256_bytes(
                (PACKAGE_ROOT / "data" / "observable_grammar.json").read_bytes()
            ),
            "nuisance_registry_sha256": sha256_bytes(
                (PACKAGE_ROOT / "data" / "nuisance_registry.json").read_bytes()
            ),
            "ranking_policy_sha256": sha256_bytes(
                (PACKAGE_ROOT / "data" / "ranking_policy.json").read_bytes()
            ),
            "exposed_data_registry_sha256": sha256_bytes(
                (PACKAGE_ROOT / "data" / "exposed_data_registry.json").read_bytes()
            ),
        },
        "slot_generation_states": slot_states,
        "descent_congruence_table": descent_congruence_table(),
        "candidates": candidates,
        "candidate_count": len(candidates),
        "scoring_boundary": {
            "physical_scoring_permitted": False,
            "comparison_access_permitted": False,
            "reason": (
                "scoring waits for the complete registry: every slot must "
                "reach GENERATED or a typed negative state, and issue 639 "
                "owns the single sealed comparison"
            ),
        },
        "target_cleanliness": {
            "public_measurement_read": False,
            "comparison_payload_read": False,
            "target_value_read": False,
        },
    }
    registry["registry_sha256"] = sha256_bytes(canonical_bytes(
        {key: value for key, value in registry.items() if key != "registry_sha256"}
    ))
    return registry


def write_registry() -> Path:
    OUTPUT_PATH.write_bytes(canonical_bytes(build_registry()))
    return OUTPUT_PATH


def verify_registry() -> None:
    committed = OUTPUT_PATH.read_bytes()
    if committed != canonical_bytes(build_registry()):
        raise SystemExit("candidate registry differs from deterministic rebuild")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        verify_registry()
        print("CANDIDATE_REGISTRY_VALID")
        return 0
    print(write_registry())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
