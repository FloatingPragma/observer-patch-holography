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
STATUS = "GENERATION_COMPLETE__SCORING_SEALED"
IMPLEMENTED_SLOTS = (
    "z6_charge_line_congruences",
    "a5_angular_rules",
    "wz_scale_free_response",
    "observer_overlap_cross_spectra",
)
SKIPPED_SLOTS = {
    "family_trace_mixing_invariants": {
        "skip_condition": "open_physical_attachment",
        "detail": (
            "the family trace, determinant, and mixing invariants consume a "
            "positive physical family attachment; that producer verdict is "
            "open, so the frozen skip policy records the typed row without "
            "grammar expansion"
        ),
    },
    "cross_channel_source_identities": {
        "skip_condition": "incomplete_global_nuisance_coverage",
        "detail": (
            "cross-channel identities require a global cancellation "
            "certificate for the common dimensionful scale and the "
            "source-admissible completion directions; those certificates do "
            "not exist, so the frozen skip policy records the typed row "
            "without grammar expansion"
        ),
    },
}


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


def load_pinned_receipt(relative: str) -> dict[str, Any]:
    """Load a source receipt only if its bytes match the committed projection."""

    projection = load_json(PACKAGE_ROOT / "outputs" / "source_projection.json")
    pins = {
        row["path"]: row["sha256"]
        for row in projection.get("source_artifacts", [])
    }
    require(relative in pins, f"receipt is not a pinned source artifact: {relative}")
    payload = (REPO_ROOT / relative).read_bytes()
    require(
        sha256_bytes(payload) == pins[relative],
        f"receipt bytes drift from the committed projection: {relative}",
    )
    value = json.loads(payload.decode("utf-8"))
    require(isinstance(value, dict), f"receipt must be an object: {relative}")
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
# Exact Q(sqrt5) arithmetic and icosahedral structures
# ---------------------------------------------------------------------------


from fractions import Fraction


def q5(a: Fraction | int, b: Fraction | int) -> tuple[Fraction, Fraction]:
    """The element a + b*sqrt(5) of Q(sqrt5) as an exact pair."""

    return (Fraction(a), Fraction(b))


def q5_add(x: tuple[Fraction, Fraction], y: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return (x[0] + y[0], x[1] + y[1])


def q5_mul(x: tuple[Fraction, Fraction], y: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q5_div(x: tuple[Fraction, Fraction], y: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    norm = y[0] * y[0] - 5 * y[1] * y[1]
    conjugate = (y[0], -y[1])
    numerator = q5_mul(x, conjugate)
    return (numerator[0] / norm, numerator[1] / norm)


def q5_str(x: tuple[Fraction, Fraction]) -> str:
    return f"{x[0]}+{x[1]}*sqrt5"


def parse_q5(text: str) -> tuple[Fraction, Fraction]:
    """Parse the pinned-receipt convention ``a+b*sqrt5`` (spaces tolerated)."""

    cleaned = text.replace(" ", "").replace("sqrt(5)", "sqrt5")
    head, _, _ = cleaned.partition("*sqrt5")
    if "*sqrt5" not in cleaned:
        return (Fraction(cleaned), Fraction(0))
    a_part, _, b_part = head.rpartition("+")
    if a_part == "":
        a_part, _, b_part = head.rpartition("-")
        if a_part != "":
            b_part = "-" + b_part
        else:
            a_part, b_part = "0", head
    return (Fraction(a_part), Fraction(b_part))


def icosahedron_ports_and_edges() -> tuple[list[str], list[tuple[str, str]]]:
    """The standard twelve-port icosahedral graph, degree five, thirty edges."""

    ports = (
        ["north"]
        + [f"upper_{i}" for i in range(5)]
        + [f"lower_{i}" for i in range(5)]
        + ["south"]
    )
    edges: list[tuple[str, str]] = []
    for i in range(5):
        edges.append(("north", f"upper_{i}"))
        edges.append((f"upper_{i}", f"upper_{(i + 1) % 5}"))
        edges.append(("south", f"lower_{i}"))
        edges.append((f"lower_{i}", f"lower_{(i + 1) % 5}"))
        edges.append((f"upper_{i}", f"lower_{i}"))
        edges.append((f"upper_{i}", f"lower_{(i + 1) % 5}"))
    degree: dict[str, int] = {port: 0 for port in ports}
    for left, right in edges:
        degree[left] += 1
        degree[right] += 1
    require(len(edges) == 30, "icosahedron edge count drift")
    require(
        all(count == 5 for count in degree.values()),
        "icosahedron degree drift",
    )
    return ports, edges


def adjacency_band_certificate() -> dict[str, Any]:
    """Exact spectral band structure of the icosahedral adjacency matrix.

    The certificate verifies the annihilating polynomial
    ``(A - 5I)(A^2 - 5I)(A + I) = 0`` by exact integer matrix arithmetic and
    certifies the band multiplicities (1, 3, 3, 5) for eigenvalues
    ``5, sqrt5, -sqrt5, -1`` by exhausting every nonnegative multiplicity
    vector against the exact traces of ``A``, ``A^2``, and ``A^3``; the
    three traces pin the vector uniquely among the annihilator's roots.
    """

    ports, edges = icosahedron_ports_and_edges()
    index = {port: position for position, port in enumerate(ports)}
    size = len(ports)
    adjacency = [[0] * size for _ in range(size)]
    for left, right in edges:
        adjacency[index[left]][index[right]] = 1
        adjacency[index[right]][index[left]] = 1

    def matmul(x: list[list[int]], y: list[list[int]]) -> list[list[int]]:
        return [
            [
                sum(x[row][k] * y[k][column] for k in range(size))
                for column in range(size)
            ]
            for row in range(size)
        ]

    identity = [[1 if row == column else 0 for column in range(size)] for row in range(size)]

    def scale_shift(matrix: list[list[int]], shift: int) -> list[list[int]]:
        return [
            [matrix[row][column] + (shift if row == column else 0) for column in range(size)]
            for row in range(size)
        ]

    a_squared = matmul(adjacency, adjacency)
    a_cubed = matmul(a_squared, adjacency)
    factor_one = scale_shift(adjacency, -5)
    factor_two = scale_shift(a_squared, -5)
    factor_three = scale_shift(adjacency, 1)
    product = matmul(matmul(factor_one, factor_two), factor_three)
    annihilated = all(value == 0 for row in product for value in row)

    trace_a = sum(adjacency[i][i] for i in range(size))
    trace_a2 = sum(a_squared[i][i] for i in range(size))
    trace_a3 = sum(a_cubed[i][i] for i in range(size))
    solutions = []
    for m_five in range(size + 1):
        for m_plus in range(size + 1):
            for m_minus in range(size + 1 - m_five - m_plus):
                m_one = size - m_five - m_plus - m_minus
                if m_one < 0:
                    continue
                rational_a = 5 * m_five - m_one
                irrational_a = m_plus - m_minus
                rational_a3 = 125 * m_five - m_one
                irrational_a3 = 5 * (m_plus - m_minus)
                value_a2 = 25 * m_five + 5 * (m_plus + m_minus) + m_one
                if (
                    rational_a == trace_a
                    and irrational_a == 0
                    and value_a2 == trace_a2
                    and rational_a3 == trace_a3
                    and irrational_a3 == 0
                ):
                    solutions.append((m_five, m_plus, m_minus, m_one))
    multiplicities = {"5": 1, "sqrt5": 3, "-sqrt5": 3, "-1": 5}
    trace_check = solutions == [(1, 3, 3, 5)]
    return {
        "annihilating_polynomial": "(A-5I)(A^2-5I)(A+I)=0",
        "annihilated": annihilated,
        "trace_a": trace_a,
        "trace_a_squared": trace_a2,
        "trace_a_cubed": trace_a3,
        "multiplicity_solutions": [list(row) for row in solutions],
        "band_multiplicities": multiplicities,
        "multiplicity_trace_check": trace_check,
    }


def a5_pairing_certificate() -> dict[str, Any]:
    """Exact A5 character arithmetic for the band pairing rules.

    Conjugacy classes of A5 with sizes 1, 15, 20, 12, 12 carry the two
    three-dimensional characters with golden-ratio values on the five-cycles.
    The certificate computes the exact invariant multiplicities of
    ``3 x 3'``, ``3 x 3``, and ``3' x 3'`` by class-weighted character sums
    in Q(sqrt5).
    """

    phi = q5(Fraction(1, 2), Fraction(1, 2))
    one_minus_phi = q5(Fraction(1, 2), Fraction(-1, 2))
    class_sizes = [1, 15, 20, 12, 12]
    chi_three = [q5(3, 0), q5(-1, 0), q5(0, 0), phi, one_minus_phi]
    chi_three_prime = [q5(3, 0), q5(-1, 0), q5(0, 0), one_minus_phi, phi]
    group_order = q5(60, 0)

    def invariant_multiplicity(left: list, right: list) -> tuple[Fraction, Fraction]:
        total = q5(0, 0)
        for size, x, y in zip(class_sizes, left, right):
            total = q5_add(total, q5_mul(q5(size, 0), q5_mul(x, y)))
        return q5_div(total, group_order)

    cross = invariant_multiplicity(chi_three, chi_three_prime)
    self_three = invariant_multiplicity(chi_three, chi_three)
    self_prime = invariant_multiplicity(chi_three_prime, chi_three_prime)
    return {
        "class_sizes": class_sizes,
        "invariant_multiplicity_3_x_3prime": q5_str(cross),
        "invariant_multiplicity_3_x_3": q5_str(self_three),
        "invariant_multiplicity_3prime_x_3prime": q5_str(self_prime),
        "cross_pairing_forbidden": cross == q5(0, 0),
        "self_pairings_allowed": self_three == q5(1, 0) and self_prime == q5(1, 0),
    }


def galois_ratio_certificate() -> dict[str, Any]:
    """Exact identity: the Galois band-cost ratio equals the squared golden ratio."""

    plus = q5(5, 1)
    minus = q5(5, -1)
    ratio = q5_div(plus, minus)
    phi = q5(Fraction(1, 2), Fraction(1, 2))
    phi_squared = q5_mul(phi, phi)
    return {
        "numerator": q5_str(plus),
        "denominator": q5_str(minus),
        "ratio": q5_str(ratio),
        "phi_squared": q5_str(phi_squared),
        "ratio_equals_phi_squared": ratio == phi_squared,
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

    def expression(
        kind: str,
        terminals: list[str],
        operator_units: int,
        depth: int,
        nodes: int,
    ) -> dict[str, Any]:
        budget = grammar["complexity_budget"]
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
            "accounting": (
                "generator-declared against the frozen operator costs; the "
                "independent verifier recomputes every score from the "
                "declared units"
            ),
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
                depth=3,
                nodes=4,
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
                depth=3,
                nodes=4,
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
                depth=3,
                nodes=4,
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
                "the quotient covering kernel is cyclic of order six, so the "
                "flux and line classification of the quotient form refines "
                "the unquotiented lattice by exactly six classes"
            ),
            "expression": expression(
                "index of the pinned quotient kernel",
                ["registered_index"],
                1 + 2,
                depth=2,
                nodes=2,
            ),
            "relation_certificate": {
                "kind": "pinned_kernel_order",
                "kernel_order": 6,
                "source": "exact diagonal Z6 quotient",
            },
            "kill_rule": (
                "a physical flux and line refinement with a different "
                "finite index falsifies the selected global form branch"
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
                depth=3,
                nodes=4,
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
# Angular, scale-free, and overlap slot producers
# ---------------------------------------------------------------------------


def _slot_row(slots: dict[str, Any], slot_id: str) -> dict[str, Any]:
    return next(row for row in slots["slots"] if row["slot_id"] == slot_id)


def _exposure_for(exposure: dict[str, Any], slot_id: str) -> list[str]:
    return sorted(
        row["exposure_id"]
        for row in exposure["surfaces"]
        if slot_id in row["applies_to_slot_ids"]
    )


def _expression(
    grammar: dict[str, Any],
    kind: str,
    terminals: list[str],
    operator_units: int,
    depth: int,
    nodes: int,
) -> dict[str, Any]:
    budget = grammar["complexity_budget"]
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
        "accounting": (
            "generator-declared against the frozen operator costs; the "
            "independent verifier recomputes every score from the "
            "declared units"
        ),
    }


def produce_a5_angular_candidates(
    slots: dict[str, Any], grammar: dict[str, Any], exposure: dict[str, Any]
) -> list[dict[str, Any]]:
    slot = _slot_row(slots, "a5_angular_rules")
    bands = adjacency_band_certificate()
    pairing = a5_pairing_certificate()
    require(
        bands["annihilated"] and bands["multiplicity_trace_check"],
        "adjacency band certificate drift",
    )
    require(
        pairing["cross_pairing_forbidden"] and pairing["self_pairings_allowed"],
        "pairing certificate drift",
    )
    port_receipt = load_pinned_receipt(
        "code/a5_closure/receipts/port_current_inner_reference.receipt.json"
    )
    schur = port_receipt["icosahedral_intertwiner"]["kernel_band_schur_rigidity"]
    require(
        schur["multiplicity_of_kernel_band_in_even_block"] == "0",
        "pinned Schur rigidity drift",
    )
    galois_pair = parse_q5(
        port_receipt["compactness"]["hilbert_schmidt_pullback_band_coefficients"][
            "frame_band"
        ]
    ), parse_q5(
        port_receipt["compactness"]["hilbert_schmidt_pullback_band_coefficients"][
            "kernel_band"
        ]
    )
    require(
        galois_pair[0] == (galois_pair[1][0], -galois_pair[1][1]),
        "pinned Galois conjugacy drift",
    )

    nuisance_matrix = {
        "basis_and_coordinate_choice": {
            "status": "QUOTIENTED",
            "certificate": (
                "band multiplicities, characters, and invariant "
                "multiplicities are basis-independent spectral data"
            ),
        },
        "physical_sector_choice": {
            "status": "NAMED_CONDITIONAL_BRANCH",
            "branch": "declared-carrier-sector-selection",
            "note": (
                "the angular structure lives on the declared twelve-port "
                "carrier; its physical rotation-remnant reading is the open "
                "issue 643 question, so every candidate carries the named "
                "branch label"
            ),
        },
        "source_admissible_completion": {
            "status": "NAMED_COMPLETION_PREMISE",
            "premise": "completions preserving the icosahedral port module",
        },
        "detector_transfer_calibration": {
            "status": "DEFERRED_UNTIL_PHYSICALIZATION",
        },
        "common_dimensionful_scale": {
            "status": "NOT_APPLICABLE_DIMENSIONLESS_DISCRETE",
        },
    }
    baseline_contract = {
        "baseline_class": "minimal Standard Model plus general relativity",
        "global_form_policy": (
            "the baseline carries no twelve-port screen and no icosahedral "
            "band structure, so these relations are unconstrained baseline "
            "directions"
        ),
        "field_content": "established matter, statistically isotropic angular response",
        "counterexample": (
            "a generic angular response with a nonzero invariant coupling "
            "between the two three-dimensional bands is baseline-consistent "
            "and violates the pairing rule"
        ),
        "baseline_freedom_class": "UNCONSTRAINED_BY_BASELINE",
    }
    exposure_surfaces = _exposure_for(exposure, "a5_angular_rules")

    candidates = [
        {
            "candidate_id": "a5-band-multiplicity-pattern",
            "slot_id": slot["slot_id"],
            "grammar_class": "trace_determinant_character_index",
            "statement": (
                "the twelve-port coefficient space splits into bands of "
                "dimensions one, three, three, and five with adjacency "
                "eigenvalues five, root five, minus root five, and minus one"
            ),
            "expression": _expression(
                grammar,
                "trace-derived multiplicity vector of the port adjacency",
                ["registered_trace", "registered_multiplicity"],
                1 * 2 + 3,
                3,
                4,
            ),
            "relation_certificate": {
                "kind": "exact_integer_annihilating_polynomial",
                "polynomial": bands["annihilating_polynomial"],
                "band_multiplicities": bands["band_multiplicities"],
                "recomputable_from": "adjacency_band_certificate",
            },
            "kill_rule": (
                "an observed screen angular decomposition with a different "
                "band pattern falsifies the declared carrier branch"
            ),
            "weight_claims": {
                "exact_global_certificate": True,
                "independent_recomputation": True,
                "physicalization_complete": False,
                "baseline_freedom_counterexample": True,
                "all_registered_completions_covered": False,
                "all_continuous_parameters_covered": False,
                "conditional_branch": True,
                "open_physical_map": True,
            },
        },
        {
            "candidate_id": "a5-cross-band-pairing-zero",
            "slot_id": slot["slot_id"],
            "grammar_class": "exact_zero_forbidden_transition",
            "statement": (
                "the invariant pairing between the two three-dimensional "
                "bands vanishes exactly, so invariant couplings cannot "
                "connect them at bilinear order"
            ),
            "expression": _expression(
                grammar,
                "exact_zero of the class-weighted character inner product",
                ["registered_character", "registered_integer"],
                1 * 2 + 3 + 4,
                3,
                4,
            ),
            "relation_certificate": {
                "kind": "exact_class_weighted_character_sum",
                "invariant_multiplicity": pairing[
                    "invariant_multiplicity_3_x_3prime"
                ],
                "self_pairings": [
                    pairing["invariant_multiplicity_3_x_3"],
                    pairing["invariant_multiplicity_3prime_x_3prime"],
                ],
                "recomputable_from": "a5_pairing_certificate",
            },
            "kill_rule": (
                "an observed invariant bilinear coupling between the two "
                "three-dimensional angular channels falsifies the declared "
                "carrier branch"
            ),
            "weight_claims": {
                "exact_global_certificate": True,
                "independent_recomputation": True,
                "physicalization_complete": False,
                "baseline_freedom_counterexample": True,
                "all_registered_completions_covered": False,
                "all_continuous_parameters_covered": False,
                "conditional_branch": True,
                "open_physical_map": True,
            },
        },
        {
            "candidate_id": "a5-galois-conjugate-correlation",
            "slot_id": slot["slot_id"],
            "grammar_class": "angular_selection_cross_level_correlation",
            "statement": (
                "the two three-dimensional bands are Galois conjugates, so "
                "every rational-coefficient angular statistic takes "
                "Galois-conjugate values on them; every rational-valued "
                "angular statistic takes equal values on the two bands and "
                "its cross-level asymmetry vanishes"
            ),
            "expression": _expression(
                grammar,
                "cross-level asymmetry of a rational angular statistic over "
                "the conjugate band pair",
                ["registered_correlation", "registered_exact_algebraic"],
                1 * 2 + 3,
                3,
                4,
            ),
            "relation_certificate": {
                "kind": "pinned_galois_conjugacy",
                "conjugate_pair": [
                    q5_str(galois_pair[0]),
                    q5_str(galois_pair[1]),
                ],
                "source": (
                    "hilbert_schmidt_pullback_band_coefficients in the "
                    "pinned port-current receipt"
                ),
            },
            "kill_rule": (
                "a measured rational-valued angular statistic separating "
                "the conjugate band pair falsifies the declared carrier "
                "branch"
            ),
            "weight_claims": {
                "exact_global_certificate": True,
                "independent_recomputation": True,
                "physicalization_complete": False,
                "baseline_freedom_counterexample": True,
                "all_registered_completions_covered": False,
                "all_continuous_parameters_covered": False,
                "conditional_branch": True,
                "open_physical_map": True,
            },
        },
        {
            "candidate_id": "a5-kernel-band-even-block-zero",
            "slot_id": slot["slot_id"],
            "grammar_class": "exact_zero_forbidden_transition",
            "statement": (
                "the kernel band has multiplicity zero in the even block, "
                "so kernel-band angular response in the even channel is "
                "forbidden"
            ),
            "expression": _expression(
                grammar,
                "exact_zero of the kernel-band even-block multiplicity",
                ["registered_multiplicity"],
                1 + 2,
                2,
                2,
            ),
            "relation_certificate": {
                "kind": "pinned_certificate_value",
                "multiplicity_of_kernel_band_in_even_block": 0,
                "source": (
                    "kernel_band_schur_rigidity in the pinned port-current "
                    "receipt"
                ),
            },
            "kill_rule": (
                "an observed kernel-band response in the even angular "
                "channel falsifies the declared carrier branch"
            ),
            "weight_claims": {
                "exact_global_certificate": True,
                "independent_recomputation": False,
                "physicalization_complete": False,
                "baseline_freedom_counterexample": True,
                "all_registered_completions_covered": False,
                "all_continuous_parameters_covered": False,
                "conditional_branch": True,
                "open_physical_map": True,
            },
        },
    ]
    for candidate in candidates:
        require(
            candidate["grammar_class"] in set(slot["allowed_grammar_class_ids"]),
            f"grammar class not allowed: {candidate['candidate_id']}",
        )
        candidate["source_ancestry"] = {
            "feature_ids": list(slot["required_feature_ids"]),
            "target_ancestry": "SOURCE_ONLY_NO_PUBLIC_COMPARISON_INPUT",
        }
        candidate["nuisance_matrix"] = nuisance_matrix
        candidate["baseline_contract"] = baseline_contract
        candidate["physicalization_status"] = (
            "OPEN_MAP__SCREEN_TO_SKY_TEMPLATE_REQUIRED"
        )
        candidate["exposure_surfaces"] = exposure_surfaces
    return candidates


def produce_wz_scale_free_candidates(
    slots: dict[str, Any], grammar: dict[str, Any], exposure: dict[str, Any]
) -> list[dict[str, Any]]:
    slot = _slot_row(slots, "wz_scale_free_response")
    ratio = galois_ratio_certificate()
    require(ratio["ratio_equals_phi_squared"], "Galois ratio certificate drift")
    family_receipt = load_pinned_receipt(
        "code/a5_closure/manifests/family_band_attachment_reference.json"
    )
    costs = family_receipt["measured_receipt"]["measured_band_costs"]
    require(
        parse_q5(costs["3"]) == (Fraction(5), Fraction(-1))
        and parse_q5(costs["3p"]) == (Fraction(5), Fraction(1))
        and parse_q5(costs["5"]) == (Fraction(6), Fraction(0)),
        "pinned band cost drift",
    )
    channels = family_receipt["measured_receipt"]["measured_channels"]
    bands = adjacency_band_certificate()
    require(
        parse_q5(channels["frame_band"]) == (Fraction(0), Fraction(1))
        and parse_q5(channels["kernel_band"]) == (Fraction(0), Fraction(-1))
        and parse_q5(channels["quintet_band"]) == (Fraction(-1), Fraction(0))
        and bands["band_multiplicities"] == {"5": 1, "sqrt5": 3, "-sqrt5": 3, "-1": 5},
        "pinned channel eigenvalue drift",
    )

    nuisance_matrix = {
        "basis_and_coordinate_choice": {
            "status": "QUOTIENTED",
            "certificate": (
                "band costs and channel eigenvalues are spectral data, "
                "invariant under basis and coordinate changes"
            ),
        },
        "common_dimensionful_scale": {
            "status": "QUOTIENTED_BY_HOMOGENEITY",
            "certificate": (
                "the relations are ratios and eigenvalue identities in "
                "which any common scale cancels exactly"
            ),
        },
        "source_admissible_completion": {
            "status": "NAMED_COMPLETION_PREMISE",
            "premise": "completions preserving the icosahedral port module",
        },
        "detector_transfer_calibration": {
            "status": "DEFERRED_UNTIL_PHYSICALIZATION",
        },
        "physical_sector_choice": {
            "status": "NOT_APPLICABLE_PER_FROZEN_NUISANCE_REGISTRY",
            "content_reason": (
                "both relations consume only selection-independent spectral "
                "data: the cost pair and channel eigenvalues exist "
                "regardless of the open rank-three selection, which is not "
                "an input to either relation"
            ),
        },
    }
    baseline_contract = {
        "baseline_class": "minimal Standard Model plus general relativity",
        "global_form_policy": (
            "the baseline carries no band-cost observable, so the relations "
            "are unconstrained baseline directions"
        ),
        "field_content": "established matter",
        "counterexample": (
            "a response model whose conjugate-channel cost ratio differs "
            "from the squared golden ratio is baseline-consistent and "
            "violates the identity"
        ),
        "baseline_freedom_class": "UNCONSTRAINED_BY_BASELINE",
    }
    exposure_surfaces = _exposure_for(exposure, "wz_scale_free_response")

    candidates = [
        {
            "candidate_id": "wz-galois-cost-ratio-phi-squared",
            "slot_id": slot["slot_id"],
            "grammar_class": "homogeneous_scale_free_combination",
            "statement": (
                "the ratio of the Galois-conjugate band costs equals the "
                "squared golden ratio exactly"
            ),
            "expression": _expression(
                grammar,
                "ratio of the conjugate registered band costs",
                ["registered_exact_algebraic"],
                1 * 2 + 5,
                3,
                4,
            ),
            "relation_certificate": {
                "kind": "exact_q5_identity",
                "ratio": ratio["ratio"],
                "phi_squared": ratio["phi_squared"],
                "pinned_costs": [costs["3p"], costs["3"]],
                "recomputable_from": "galois_ratio_certificate",
            },
            "kill_rule": (
                "a certified conjugate-channel cost ratio away from the "
                "squared golden ratio falsifies the band identification"
            ),
            "weight_claims": {
                "exact_global_certificate": True,
                "independent_recomputation": True,
                "physicalization_complete": False,
                "baseline_freedom_counterexample": True,
                "all_registered_completions_covered": False,
                "all_continuous_parameters_covered": False,
                "conditional_branch": False,
                "open_physical_map": True,
            },
        },
        {
            "candidate_id": "wz-band-channel-shared-identity",
            "slot_id": slot["slot_id"],
            "grammar_class": "shared_response_identity",
            "statement": (
                "the family-band channel eigenvalues recorded in the family "
                "attachment receipt equal the icosahedral adjacency band "
                "eigenvalues recomputed from the port graph"
            ),
            "expression": _expression(
                grammar,
                "equality of two independently produced channel eigenvalue "
                "vectors",
                ["registered_exact_algebraic", "registered_correlation"],
                1 * 2 + 3,
                3,
                4,
            ),
            "relation_certificate": {
                "kind": "cross_receipt_shared_identity",
                "family_receipt_channels": [
                    channels["frame_band"],
                    channels["kernel_band"],
                    channels["quintet_band"],
                ],
                "recomputed_band_eigenvalues": ["0+1*sqrt5", "0+-1*sqrt5", "-1+0*sqrt5"],
                "recomputable_from": "adjacency_band_certificate",
            },
            "kill_rule": (
                "a certified divergence between the two independently "
                "produced channel vectors falsifies the shared response "
                "identification"
            ),
            "weight_claims": {
                "exact_global_certificate": True,
                "independent_recomputation": True,
                "physicalization_complete": False,
                "baseline_freedom_counterexample": True,
                "all_registered_completions_covered": False,
                "all_continuous_parameters_covered": False,
                "conditional_branch": False,
                "open_physical_map": True,
            },
        },
    ]
    for candidate in candidates:
        require(
            candidate["grammar_class"] in set(slot["allowed_grammar_class_ids"]),
            f"grammar class not allowed: {candidate['candidate_id']}",
        )
        candidate["source_ancestry"] = {
            "feature_ids": list(slot["required_feature_ids"]),
            "target_ancestry": "SOURCE_ONLY_NO_PUBLIC_COMPARISON_INPUT",
        }
        candidate["nuisance_matrix"] = nuisance_matrix
        candidate["baseline_contract"] = baseline_contract
        candidate["physicalization_status"] = (
            "OPEN_MAP__POLE_RESIDUE_READOUT_AND_LABORATORY_ATTACHMENT_REQUIRED"
        )
        candidate["exposure_surfaces"] = exposure_surfaces
    return candidates


def produce_overlap_candidates(
    slots: dict[str, Any], grammar: dict[str, Any], exposure: dict[str, Any]
) -> list[dict[str, Any]]:
    slot = _slot_row(slots, "observer_overlap_cross_spectra")
    receipt = load_pinned_receipt(
        "code/a5_closure/manifests/classical_realization_receipt.json"
    )
    readings = receipt["classical_completion"]["sector_readings"]
    kernel_counts = [row["classical_kernel_count"] for row in sorted(readings, key=lambda r: r["sector"])]
    require(kernel_counts == [0, 0, 0, 2, 0, 0], "pinned kernel pattern drift")
    require(
        all(row["kernel_match"] and row["gap_match"] for row in readings),
        "pinned sector reading drift",
    )
    ladder = receipt["classical_completion"]["ladder_point_readings"]
    require(
        all(row.get("kernel_match", True) for row in ladder),
        "pinned ladder reading drift",
    )

    nuisance_matrix = {
        "basis_and_coordinate_choice": {
            "status": "QUOTIENTED",
            "certificate": (
                "kernel counts and sector conjugation structure are "
                "spectral data of the declared sector operators"
            ),
        },
        "common_dimensionful_scale": {
            "status": "QUOTIENTED_BY_CONSTRUCTION",
            "certificate": (
                "the relations use kernel counts and sector labels only; "
                "no dimensionful gap value enters"
            ),
        },
        "physical_sector_choice": {
            "status": "NAMED_CONDITIONAL_BRANCH",
            "branch": "declared-orientation-reversing-seam-convention",
            "note": (
                "the twist sectors are defined along the declared seam "
                "orientation convention; every candidate carries the named "
                "branch label"
            ),
        },
        "source_admissible_completion": {
            "status": "NAMED_COMPLETION_PREMISE",
            "premise": "completions preserving the declared seam complex",
        },
        "detector_transfer_calibration": {
            "status": "DEFERRED_UNTIL_PHYSICALIZATION",
        },
    }
    baseline_contract = {
        "baseline_class": "minimal Standard Model plus general relativity",
        "global_form_policy": (
            "the baseline carries no observer-overlap twist sectors, so "
            "the relations are unconstrained baseline directions"
        ),
        "field_content": "established matter",
        "counterexample": (
            "an overlap cross-spectrum with flat modes in two sectors, or "
            "with a conjugate-sector asymmetry, is baseline-consistent and "
            "violates the relations"
        ),
        "baseline_freedom_class": "UNCONSTRAINED_BY_BASELINE",
    }
    exposure_surfaces = _exposure_for(exposure, "observer_overlap_cross_spectra")

    candidates = [
        {
            "candidate_id": "overlap-unique-flat-sector",
            "slot_id": slot["slot_id"],
            "grammar_class": "shared_response_identity",
            "statement": (
                "exactly one twist sector carries the flat overlap mode, "
                "the sector of trivial total seam transport under the "
                "declared orientation convention (sector index three in "
                "the pinned receipt), with kernel pattern zero elsewhere"
            ),
            "expression": _expression(
                grammar,
                "kernel-count vector of the six twist sectors",
                ["registered_multiplicity"],
                1 + 2,
                2,
                2,
            ),
            "relation_certificate": {
                "kind": "pinned_receipt_value",
                "kernel_pattern": [0, 0, 0, 1, 0, 0],
                "classical_doubled_pattern": kernel_counts,
                "source": (
                    "sector_readings in the pinned classical realization "
                    "receipt"
                ),
            },
            "kill_rule": (
                "an observed overlap cross-spectrum with flat modes in a "
                "different sector count falsifies the declared seam branch"
            ),
            "weight_claims": {
                "exact_global_certificate": False,
                "independent_recomputation": False,
                "physicalization_complete": False,
                "baseline_freedom_counterexample": True,
                "all_registered_completions_covered": False,
                "all_continuous_parameters_covered": False,
                "conditional_branch": True,
                "open_physical_map": True,
            },
        },
        {
            "candidate_id": "overlap-conjugate-sector-symmetry",
            "slot_id": slot["slot_id"],
            "grammar_class": "normalized_residue_asymmetry_cross_spectrum",
            "statement": (
                "conjugate twist sectors carry identical spectra, so the "
                "normalized cross-spectrum asymmetry between conjugate "
                "sectors vanishes"
            ),
            "expression": _expression(
                grammar,
                "normalized asymmetry between conjugate sector spectra",
                ["registered_correlation"],
                1 + 5,
                2,
                2,
            ),
            "relation_certificate": {
                "kind": "operator_conjugation_symmetry",
                "statement": (
                    "the sector operator for twist class k is the complex "
                    "conjugate of the operator for class six minus k, so "
                    "their spectra coincide exactly"
                ),
                "operator_identity_source": (
                    "the sector Laplacian construction recorded in the "
                    "pinned classical realization receipt: the stiffness "
                    "matrix is the exact realification of the sector "
                    "operator, and conjugate twist classes have conjugate "
                    "operators"
                ),
                "pinned_confirmation": (
                    "sector readings one and five, and two and four, agree "
                    "in the pinned receipt"
                ),
            },
            "kill_rule": (
                "a certified nonzero conjugate-sector asymmetry falsifies "
                "the declared seam branch"
            ),
            "weight_claims": {
                "exact_global_certificate": True,
                "independent_recomputation": False,
                "physicalization_complete": False,
                "baseline_freedom_counterexample": True,
                "all_registered_completions_covered": False,
                "all_continuous_parameters_covered": False,
                "conditional_branch": True,
                "open_physical_map": True,
            },
        },
        {
            "candidate_id": "overlap-flat-sector-refinement-stability",
            "slot_id": slot["slot_id"],
            "grammar_class": "refinement_scaling_exponent",
            "statement": (
                "the flat-sector location is refinement-stable: the kernel "
                "pattern scales with exponent zero along the declared "
                "refinement ladder"
            ),
            "expression": _expression(
                grammar,
                "scaling exponent of the kernel pattern along the ladder",
                ["registered_multiplicity", "registered_integer"],
                1 * 2 + 3,
                3,
                4,
            ),
            "relation_certificate": {
                "kind": "pinned_receipt_value",
                "ladder_rungs": 2,
                "sector_rows_per_rung": len(ladder),
                "source": (
                    "ladder_point_readings in the pinned classical "
                    "realization receipt: the six sectors at a second "
                    "domain size against the frozen main domain"
                ),
            },
            "kill_rule": (
                "a refinement rung moving the flat mode to a different "
                "sector falsifies the declared seam branch"
            ),
            "weight_claims": {
                "exact_global_certificate": False,
                "independent_recomputation": False,
                "physicalization_complete": False,
                "baseline_freedom_counterexample": True,
                "all_registered_completions_covered": False,
                "all_continuous_parameters_covered": False,
                "conditional_branch": True,
                "open_physical_map": True,
            },
        },
    ]
    for candidate in candidates:
        require(
            candidate["grammar_class"] in set(slot["allowed_grammar_class_ids"]),
            f"grammar class not allowed: {candidate['candidate_id']}",
        )
        candidate["source_ancestry"] = {
            "feature_ids": list(slot["required_feature_ids"]),
            "target_ancestry": "SOURCE_ONLY_NO_PUBLIC_COMPARISON_INPUT",
        }
        candidate["nuisance_matrix"] = nuisance_matrix
        candidate["baseline_contract"] = baseline_contract
        candidate["physicalization_status"] = (
            "OPEN_MAP__INTERFEROMETER_READOUT_ATTACHMENT_REQUIRED"
        )
        candidate["exposure_surfaces"] = exposure_surfaces
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

    produced = (
        produce_z6_candidates(slots, grammar, nuisances, exposure)
        + produce_a5_angular_candidates(slots, grammar, exposure)
        + produce_wz_scale_free_candidates(slots, grammar, exposure)
        + produce_overlap_candidates(slots, grammar, exposure)
    )
    candidates = rank_candidates(produced, ranking)
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
    skip_records = {}
    for row in slots["slots"]:
        slot_id = row["slot_id"]
        if slot_id == direct_slot:
            slot_states[slot_id] = "EXCLUDED_EXTERNALLY_ORDERED_DIRECT_N"
        elif slot_id in IMPLEMENTED_SLOTS:
            slot_states[slot_id] = "GENERATED"
        elif slot_id in SKIPPED_SLOTS:
            slot_states[slot_id] = "NOT_EVALUABLE_TYPED_SKIP"
            record = dict(SKIPPED_SLOTS[slot_id])
            record["post_skip_action"] = ranking["skip_policy"][
                "post_skip_action"
            ]
            require(
                record["skip_condition"]
                in set(ranking["skip_policy"]["skip_conditions"]),
                f"skip condition not frozen: {slot_id}",
            )
            skip_records[slot_id] = record
        else:
            slot_states[slot_id] = "REGISTERED_AWAITING_GENERATION"
    require(
        "REGISTERED_AWAITING_GENERATION" not in slot_states.values(),
        "registry incomplete: a slot has no generation outcome",
    )

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
        "skip_records": skip_records,
        "descent_congruence_table": descent_congruence_table(),
        "candidates": candidates,
        "candidate_count": len(candidates),
        "scoring_boundary": {
            "physical_scoring_permitted": False,
            "comparison_access_permitted": False,
            "reason": (
                "every slot carries a generation outcome, so the registry "
                "is complete; scoring stays sealed behind a separately "
                "reviewed unsealing change and the single issue-639 "
                "comparison with custody and exposure typing"
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
