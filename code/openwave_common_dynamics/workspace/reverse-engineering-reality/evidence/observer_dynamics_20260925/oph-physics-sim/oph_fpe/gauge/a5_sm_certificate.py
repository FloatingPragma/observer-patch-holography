from __future__ import annotations

import copy
import hashlib
import itertools
import json
import math
import re
from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Sequence

from jsonschema import Draft202012Validator

from oph_fpe.claims import CONTINUATION, with_claim_metadata


A5_TWELVE_PORT_STRUCTURAL_RECEIPT = "A5_TWELVE_PORT_STRUCTURAL_RECEIPT"
SM_ADJOINT_CHARACTER_MATCH_RECEIPT = "SM_ADJOINT_CHARACTER_MATCH_RECEIPT"
CONDITIONAL_EXTERIOR_ONE_GENERATION_WITNESS_RECEIPT = (
    "CONDITIONAL_EXTERIOR_ONE_GENERATION_WITNESS_RECEIPT"
)
NO_INVARIANT_PORT_PARTITION_8_3_1_RECEIPT = (
    "NO_INVARIANT_PORT_PARTITION_8_3_1_RECEIPT"
)
PHYSICAL_STANDARD_MODEL_PROMOTION_RECEIPT = (
    "PHYSICAL_STANDARD_MODEL_PROMOTION_RECEIPT"
)

REPORT_SCHEMA = "oph_a5_twelve_port_structural_certificate_v2"
VERIFICATION_SCHEMA = "oph_a5_twelve_port_structural_verification_v2"
DEFAULT_SCHEMA_PATH = (
    Path(__file__).resolve().parents[2]
    / "schemas"
    / "gauge"
    / "a5_sm_certificate.schema.json"
)

Permutation = tuple[int, int, int, int, int]
ActionRow = tuple[int, ...]


@dataclass(frozen=True)
class _Q5:
    """Exact element of Q(sqrt(5)), represented as a + b sqrt(5)."""

    rational: Fraction = Fraction(0)
    sqrt5: Fraction = Fraction(0)

    def __add__(self, other: _Q5 | int | Fraction) -> _Q5:
        rhs = _as_q5(other)
        return _Q5(self.rational + rhs.rational, self.sqrt5 + rhs.sqrt5)

    def __radd__(self, other: _Q5 | int | Fraction) -> _Q5:
        return self + other

    def __sub__(self, other: _Q5 | int | Fraction) -> _Q5:
        rhs = _as_q5(other)
        return _Q5(self.rational - rhs.rational, self.sqrt5 - rhs.sqrt5)

    def __rsub__(self, other: _Q5 | int | Fraction) -> _Q5:
        return _as_q5(other) - self

    def __mul__(self, other: _Q5 | int | Fraction) -> _Q5:
        rhs = _as_q5(other)
        return _Q5(
            self.rational * rhs.rational + 5 * self.sqrt5 * rhs.sqrt5,
            self.rational * rhs.sqrt5 + self.sqrt5 * rhs.rational,
        )

    def __rmul__(self, other: _Q5 | int | Fraction) -> _Q5:
        return self * other

    def __truediv__(self, other: int | Fraction) -> _Q5:
        divisor = Fraction(other)
        return _Q5(self.rational / divisor, self.sqrt5 / divisor)


_ZERO = _Q5()
_ONE = _Q5(Fraction(1))
_PHI = _Q5(Fraction(1, 2), Fraction(1, 2))
_PHI_BAR = _Q5(Fraction(1, 2), Fraction(-1, 2))
_CLASS_ORDER = ("1A", "2A", "3A", "5A", "5B")
_CLASS_SIZES = (1, 15, 20, 12, 12)
_IRREPS: dict[str, tuple[_Q5, ...]] = {
    "1": (_Q5(1), _Q5(1), _Q5(1), _Q5(1), _Q5(1)),
    "3": (_Q5(3), _Q5(-1), _ZERO, _PHI, _PHI_BAR),
    "3prime": (_Q5(3), _Q5(-1), _ZERO, _PHI_BAR, _PHI),
    "4": (_Q5(4), _ZERO, _Q5(1), _Q5(-1), _Q5(-1)),
    "5": (_Q5(5), _Q5(1), _Q5(-1), _ZERO, _ZERO),
}

_PHYSICAL_GATES = {
    "UD12_SOURCE_RECEIPT": "source-derived feasible twelve-unit defect domain and strict unit-splitting cost",
    "RP_A5_SOURCE_SELECTOR_RECEIPT": "source-derived icosahedral selector without an icosahedral template",
    "PORT_CURRENT_INNER_RECEIPT": "full-rank compact skew-adjoint commutator-closed current algebra with inner A5 action",
    "PORT_REFINEMENT_INTERTWINER_RECEIPT": "refinement-stable port-to-current intertwiner",
    "BLOCK_DETERMINANT_BALANCE_RECEIPT": "physical trace-balanced block descent",
    "PORT_SPIN_LIFT_RECEIPT": "physical compact-group spin lift",
    "AXIS_CENTER_DESCENT_RECEIPT": "physical identification of the six-axis residue with the gauge-group center quotient",
    "MATTER_REALIZATION_RECEIPT": "matter, hypercharge, and faithful tensor-kernel realization inside a source-complete menu",
    "EXTERIOR_PACKAGE_SELECTION_RECEIPT": "source-derived selection of the non-vacuum even exterior package",
    "HIGGS_SCALAR_SELECTION_RECEIPT": "source-derived selection of W as the physical scalar doublet",
    "A5_FAMILY_ATTACHMENT_RECEIPT": "attachment of the face-phase A5 carrier to physical families",
    "A5_FAMILY_DESCENT_RECEIPT": "derived breaking, hiding, or forgetting of A5 before general Yukawa matrices",
    "PORT_WEAK_INTERTWINER_RECEIPT": "physical intertwiner from screen invariants to weak-doublet copies",
    "PORT_LOAD_TRACE_RECEIPT": "physical identification of the normalized four-copy weak load",
    "CONTINUUM_LIMIT_RECEIPT": "controlled continuum limit of the finite carrier/current system",
    "SPIN_QFT_REALIZATION_RECEIPT": "local spin-statistics and quantum-field-theory realization",
    "NO_EXTRA_LIGHT_SECTOR_RECEIPT": "observer-visible exclusion of the vacuum singlet and other anomaly-free light sectors",
}

_EXTERIOR_CHECK_NAMES = (
    "exterior_carrier_trace_is_balanced",
    "exterior_package_dimension_is_15",
    "exterior_field_dimensions_and_hypercharges_match",
    "exterior_package_is_chiral_at_representation_level",
    "three_one_higgs_invariant_lines",
    "listed_perturbative_anomalies_cancel",
    "su2_witten_parity_is_even",
    "weak_doublet_count_is_four",
)

_EXPECTED_ONE_HIGGS_CANDIDATE_SPECS = (
    ("Q_H_Q", "Q", "H", "Q"),
    ("Q_Hdag_Q", "Q", "Hdag", "Q"),
    ("Q_H_L", "Q", "H", "L"),
    ("Q_Hdag_L", "Q", "Hdag", "L"),
    ("Q_H_u_c", "Q", "H", "u_c"),
    ("Q_Hdag_u_c", "Q", "Hdag", "u_c"),
    ("Q_H_d_c", "Q", "H", "d_c"),
    ("Q_Hdag_d_c", "Q", "Hdag", "d_c"),
    ("Q_H_e_c", "Q", "H", "e_c"),
    ("Q_Hdag_e_c", "Q", "Hdag", "e_c"),
    ("L_H_L", "L", "H", "L"),
    ("L_Hdag_L", "L", "Hdag", "L"),
    ("L_H_u_c", "L", "H", "u_c"),
    ("L_Hdag_u_c", "L", "Hdag", "u_c"),
    ("L_H_d_c", "L", "H", "d_c"),
    ("L_Hdag_d_c", "L", "Hdag", "d_c"),
    ("L_H_e_c", "L", "H", "e_c"),
    ("L_Hdag_e_c", "L", "Hdag", "e_c"),
    ("u_c_H_u_c", "u_c", "H", "u_c"),
    ("u_c_Hdag_u_c", "u_c", "Hdag", "u_c"),
    ("u_c_H_d_c", "u_c", "H", "d_c"),
    ("u_c_Hdag_d_c", "u_c", "Hdag", "d_c"),
    ("u_c_H_e_c", "u_c", "H", "e_c"),
    ("u_c_Hdag_e_c", "u_c", "Hdag", "e_c"),
    ("d_c_H_d_c", "d_c", "H", "d_c"),
    ("d_c_Hdag_d_c", "d_c", "Hdag", "d_c"),
    ("d_c_H_e_c", "d_c", "H", "e_c"),
    ("d_c_Hdag_e_c", "d_c", "Hdag", "e_c"),
    ("e_c_H_e_c", "e_c", "H", "e_c"),
    ("e_c_Hdag_e_c", "e_c", "Hdag", "e_c"),
)


def a5_sm_structural_certificate(
    imported_theorem_evidence: Sequence[Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Recompute the exact finite A5/twelve-port structural certificate.

    Imported theorem records are retained as non-promoting provenance. They are
    never used to compute a native receipt. The exact module-level match is not
    promoted to a physical port/current identification.
    """

    native = copy.deepcopy(_native_computation())
    checks = dict(native.pop("checks"))
    structural_receipt = all(checks.values())
    sm_match = bool(checks["sm_adjoint_character_matches_port_module"])
    partition_no_go = bool(checks["no_invariant_8_3_1_port_partition"])
    exterior_witness = all(checks[name] for name in _EXTERIOR_CHECK_NAMES)
    imports = _normalize_imported_evidence(imported_theorem_evidence)
    physical = _physical_promotion(
        structural_receipt=structural_receipt,
        sm_match=sm_match,
        partition_no_go=partition_no_go,
        exterior_witness=exterior_witness,
    )

    report: dict[str, Any] = {
        "schema": REPORT_SCHEMA,
        "mode": REPORT_SCHEMA,
        A5_TWELVE_PORT_STRUCTURAL_RECEIPT: structural_receipt,
        SM_ADJOINT_CHARACTER_MATCH_RECEIPT: sm_match,
        CONDITIONAL_EXTERIOR_ONE_GENERATION_WITNESS_RECEIPT: exterior_witness,
        NO_INVARIANT_PORT_PARTITION_8_3_1_RECEIPT: partition_no_go,
        PHYSICAL_STANDARD_MODEL_PROMOTION_RECEIPT: False,
        "receipt": structural_receipt,
        "checks": checks,
        "native_computation": native,
        "evidence_layers": {
            "native_exact_computation": {
                "producer": "oph_fpe.gauge.a5_sm_certificate",
                "method": "exact_integer_permutation_coset_character_matrix_and_fractional_exterior_arithmetic",
                "depends_on_imported_theorem_evidence": False,
                "simulation_receipt_eligible": True,
                "receipt": structural_receipt,
            },
            "imported_theorem_evidence": imports,
        },
        "physical_promotion": physical,
        "blockers": physical["blockers"],
        "claim_boundary": (
            "This certificate recomputes exact finite A5, A5/C5 twelve-port, character, "
            "icosahedral-adjacency, Standard-Model-adjoint restriction, and invariant-set "
            "anti-bridge arithmetic, together with the conditional exterior one-generation "
            "field and hypercharge rows, an exhaustive 15-unordered-fermion-pair by two-scalar "
            "ledger whose only three gauge-invariant lines are Q-H-u_c, Q-Hdag-d_c, and "
            "L-Hdag-e_c, and the anomaly and weak-parity witness. The "
            "character match is an isomorphism of real A5 "
            "modules, not an A5-invariant partition of twelve named ports into 8+3+1 "
            "generators and not a physical port/current map. Imported theorem evidence is "
            "provenance-only. Physical Standard Model promotion remains false until every "
            "listed source, current, refinement, global-form, and matter receipt is produced "
            "and independently verified."
        ),
    }
    report = with_claim_metadata(
        report,
        claim_level=CONTINUATION,
        receipt=A5_TWELVE_PORT_STRUCTURAL_RECEIPT,
        physical_claim=False,
        observable_id="a5_c5_twelve_port_structural_closure",
        fit_objective="exact_finite_structural_identity_without_physical_promotion",
    )
    report["certificate_payload_sha256"] = _payload_sha256(report)
    return report


def verify_a5_sm_structural_certificate(
    report: Mapping[str, Any],
    *,
    schema_path: Path | None = None,
) -> dict[str, Any]:
    """Independently recompute critical fields and reject producer tampering."""

    blockers: list[str] = []
    candidate = dict(report) if isinstance(report, Mapping) else {}
    path = Path(schema_path) if schema_path is not None else DEFAULT_SCHEMA_PATH
    try:
        schema = json.loads(path.read_text(encoding="utf-8"))
        errors = sorted(
            Draft202012Validator(schema).iter_errors(candidate),
            key=lambda error: tuple(str(part) for part in error.absolute_path),
        )
        blockers.extend(
            f"schema:{'/'.join(str(part) for part in error.absolute_path) or '<root>'}:{error.message}"
            for error in errors
        )
    except (OSError, json.JSONDecodeError) as exc:
        blockers.append(f"schema_unavailable_or_invalid:{type(exc).__name__}")

    expected_native = copy.deepcopy(_native_computation())
    expected_checks = dict(expected_native.pop("checks"))
    expected_structural = all(expected_checks.values())
    expected_sm = expected_checks["sm_adjoint_character_matches_port_module"]
    expected_no_go = expected_checks["no_invariant_8_3_1_port_partition"]
    expected_exterior = all(expected_checks[name] for name in _EXTERIOR_CHECK_NAMES)
    expected_physical = _physical_promotion(
        structural_receipt=expected_structural,
        sm_match=expected_sm,
        partition_no_go=expected_no_go,
        exterior_witness=expected_exterior,
    )
    expected_native_evidence = {
        "producer": "oph_fpe.gauge.a5_sm_certificate",
        "method": "exact_integer_permutation_coset_character_matrix_and_fractional_exterior_arithmetic",
        "depends_on_imported_theorem_evidence": False,
        "simulation_receipt_eligible": True,
        "receipt": expected_structural,
    }
    candidate_evidence = candidate.get("evidence_layers")
    candidate_native_evidence = (
        candidate_evidence.get("native_exact_computation")
        if isinstance(candidate_evidence, Mapping)
        else None
    )

    comparisons = {
        "native_computation_matches_recomputation": candidate.get("native_computation") == expected_native,
        "checks_match_recomputation": candidate.get("checks") == expected_checks,
        "structural_receipt_matches_recomputation": candidate.get(
            A5_TWELVE_PORT_STRUCTURAL_RECEIPT
        )
        is expected_structural,
        "sm_match_receipt_matches_recomputation": candidate.get(
            SM_ADJOINT_CHARACTER_MATCH_RECEIPT
        )
        is expected_sm,
        "exterior_witness_receipt_matches_recomputation": candidate.get(
            CONDITIONAL_EXTERIOR_ONE_GENERATION_WITNESS_RECEIPT
        )
        is expected_exterior,
        "partition_no_go_receipt_matches_recomputation": candidate.get(
            NO_INVARIANT_PORT_PARTITION_8_3_1_RECEIPT
        )
        is expected_no_go,
        "receipt_alias_matches_recomputation": candidate.get("receipt") is expected_structural,
        "native_evidence_layer_matches_recomputation": candidate_native_evidence
        == expected_native_evidence,
        "physical_promotion_matches_fail_closed_contract": candidate.get("physical_promotion")
        == expected_physical,
        "top_level_blockers_match_fail_closed_contract": candidate.get("blockers")
        == expected_physical["blockers"],
        "physical_promotion_receipt_is_false": candidate.get(
            PHYSICAL_STANDARD_MODEL_PROMOTION_RECEIPT
        )
        is False,
        "physical_claim_is_false": candidate.get("physical_claim") is False,
        "payload_hash_matches": candidate.get("certificate_payload_sha256")
        == _payload_sha256(candidate),
        "imports_do_not_promote": _imports_are_nonpromoting(candidate),
        "import_summary_is_internally_consistent": _import_summary_is_consistent(candidate),
    }
    blockers.extend(name for name, accepted in comparisons.items() if not accepted)
    return {
        "schema": VERIFICATION_SCHEMA,
        "verified": not blockers,
        "checks": comparisons,
        "blockers": blockers,
        "claim_boundary": (
            "Verification recomputes the finite structural evidence and enforces non-promotion "
            "of imported theorem records and physical Standard Model claims."
        ),
    }


def write_a5_sm_structural_certificate(
    out: Path,
    imported_theorem_evidence: Sequence[Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    report = a5_sm_structural_certificate(imported_theorem_evidence)
    verification = verify_a5_sm_structural_certificate(report)
    if verification["verified"] is not True:
        raise ValueError(f"A5 structural certificate failed verification: {verification['blockers']}")
    destination = Path(out)
    if destination.suffix.lower() != ".json":
        destination = destination / "a5_sm_structural_certificate.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


@lru_cache(maxsize=1)
def _native_computation() -> dict[str, Any]:
    group = tuple(permutation for permutation in itertools.permutations(range(5)) if _is_even(permutation))
    group_set = frozenset(group)
    identity: Permutation = (0, 1, 2, 3, 4)
    inverses = {element: _inverse(element) for element in group}
    order_profile = Counter(_permutation_order(element) for element in group)
    conjugacy = _label_conjugacy_classes(group, inverses)

    rot5: Permutation = (1, 2, 3, 4, 0)
    h5 = frozenset(_power(rot5, exponent) for exponent in range(5))
    cosets = _left_cosets(group, h5)
    element_to_coset = {
        element: index for index, coset in enumerate(cosets) for element in coset
    }
    representatives = tuple(min(coset) for coset in cosets)
    actions = {
        element: tuple(
            element_to_coset[_compose(element, representative)]
            for representative in representatives
        )
        for element in group
    }
    base_vertex = element_to_coset[identity]
    stabilizer = frozenset(
        element for element, action in actions.items() if action[base_vertex] == base_vertex
    )
    kernel = frozenset(
        element
        for element, action in actions.items()
        if action == tuple(range(len(cosets)))
    )
    orbit = frozenset(action[base_vertex] for action in actions.values())
    character = {
        label: _class_fixed_point_character(elements, actions)
        for label, elements in conjugacy.items()
    }
    character_constant = all(row["constant_on_class"] for row in character.values())
    permutation_character = tuple(
        _Q5(character[label]["fixed_cosets"]) for label in _CLASS_ORDER
    )

    character_certificate = _character_certificate(permutation_character)
    adjacency = _adjacency_certificate(
        group=group,
        h5=h5,
        actions=actions,
        base_vertex=base_vertex,
        vertex_count=len(cosets),
    )
    exterior = _exterior_one_generation_witness()
    invariant_sizes = _invariant_subset_size_counts(tuple(actions.values()), len(cosets))

    closure = all(_compose(left, right) in group_set for left in group for right in group)
    action_homomorphism = all(
        actions[_compose(left, right)] == _compose_action(actions[left], actions[right])
        for left in group
        for right in group
    )
    checks = {
        "a5_has_60_even_permutations": len(group) == 60 and all(_is_even(element) for element in group),
        "a5_is_closed_with_identity_and_inverses": closure
        and identity in group_set
        and all(inverses[element] in group_set for element in group),
        "a5_element_order_profile_is_1_15_20_24": order_profile
        == Counter({1: 1, 2: 15, 3: 20, 5: 24}),
        "a5_conjugacy_class_profile_is_1_15_20_12_12": tuple(
            len(conjugacy[label]) for label in _CLASS_ORDER
        )
        == _CLASS_SIZES,
        "h5_is_order_five_subgroup": len(h5) == 5
        and identity in h5
        and all(_compose(left, right) in h5 for left in h5 for right in h5),
        "a5_over_h5_has_twelve_cosets": len(cosets) == 12
        and sum(len(coset) for coset in cosets) == 60,
        "coset_action_is_a_homomorphism": action_homomorphism
        and all(sorted(action) == list(range(12)) for action in actions.values()),
        "coset_action_is_transitive": len(orbit) == 12,
        "base_stabilizer_is_h5_of_order_five": stabilizer == h5 and len(stabilizer) == 5,
        "orbit_stabilizer_is_60_equals_12_times_5": len(group) == len(orbit) * len(stabilizer),
        "coset_action_is_faithful": kernel == frozenset({identity}),
        "permutation_character_is_12_0_0_2_2": tuple(
            character[label]["fixed_cosets"] for label in _CLASS_ORDER
        )
        == (12, 0, 0, 2, 2)
        and character_constant,
        "a5_character_table_is_orthonormal_and_complete": character_certificate[
            "character_table_checks"
        ]["orthonormal"]
        and character_certificate["character_table_checks"]["dimension_squares_sum_to_60"],
        "port_character_decomposes_as_1_3_3prime_5": character_certificate[
            "port_module"
        ]["decomposition"]
        == {"1": 1, "3": 1, "3prime": 1, "4": 0, "5": 1},
        "port_character_is_multiplicity_free": character_certificate["port_module"][
            "multiplicity_free"
        ],
        "icosahedral_orbital_graph_is_5_regular_connected": adjacency["edge_count"] == 30
        and adjacency["degree_sequence"] == [5] * 12
        and adjacency["connected"]
        and adjacency["a5_invariant"],
        "icosahedral_adjacency_characteristic_polynomial_matches": adjacency[
            "characteristic_polynomial"
        ]["factorization_verified"],
        "icosahedral_adjacency_minimal_polynomial_annihilates": adjacency[
            "minimal_polynomial"
        ]["annihilator_is_zero"],
        "icosahedral_spectrum_has_ranks_1_3_3_5": adjacency["spectrum"]
        == {"5": 1, "sqrt(5)": 3, "-sqrt(5)": 3, "-1": 5},
        "sm_adjoint_character_matches_port_module": character_certificate[
            "sm_adjoint_match"
        ]["both_opposite_triplet_assignments_match"],
        "no_invariant_8_3_1_port_partition": sorted(invariant_sizes) == [0, 12]
        and not any(size in invariant_sizes for size in (1, 3, 8)),
        "exterior_carrier_trace_is_balanced": exterior["carrier"]["trace_balance"] == "0"
        and exterior["carrier"]["trace_balanced"],
        "exterior_package_dimension_is_15": exterior["matter_package"]["complex_dimension"] == 15
        and exterior["matter_package"]["lambda2_dimension"] == 10
        and exterior["matter_package"]["lambda4_dimension"] == 5,
        "exterior_field_dimensions_and_hypercharges_match": exterior["matter_package"][
            "field_signature_matches_canonical_generation"
        ],
        "exterior_package_is_chiral_at_representation_level": exterior["matter_package"][
            "disjoint_from_dual"
        ],
        "three_one_higgs_invariant_lines": exterior["one_higgs_invariant_lines"][
            "invariant_line_count"
        ]
        == 3
        and exterior["one_higgs_invariant_lines"]["all_channels_exact_singlets"]
        and exterior["one_higgs_invariant_lines"][
            "exhaustive_candidate_ledger_complete"
        ]
        and exterior["one_higgs_invariant_lines"][
            "nonzero_set_equals_three_canonical_channels"
        ],
        "listed_perturbative_anomalies_cancel": exterior["anomalies"]["all_coefficients_zero"],
        "su2_witten_parity_is_even": exterior["anomalies"]["su2_witten_parity_even"],
        "weak_doublet_count_is_four": exterior["weak_doublet_count"][
            "multiplicity_per_generation"
        ]
        == 4,
    }

    class_rows = []
    for label in _CLASS_ORDER:
        elements = conjugacy[label]
        representative = min(elements)
        class_rows.append(
            {
                "label": label,
                "element_order": _permutation_order(representative),
                "size": len(elements),
                "representative": list(representative),
                "fixed_cosets": character[label]["fixed_cosets"],
                "fixed_count_constant_on_class": character[label]["constant_on_class"],
            }
        )

    return {
        "checks": checks,
        "a5_group": {
            "presentation": "even_permutations_of_five_letters",
            "order": len(group),
            "identity": list(identity),
            "element_order_counts": {str(order): count for order, count in sorted(order_profile.items())},
            "conjugacy_classes": class_rows,
            "canonical_group_sha256": _sha256_json([list(element) for element in group]),
        },
        "h5_coset_action": {
            "generator": list(rot5),
            "subgroup_order": len(h5),
            "coset_count": len(cosets),
            "action_transitive": len(orbit) == len(cosets),
            "base_vertex": base_vertex,
            "base_stabilizer_order": len(stabilizer),
            "base_stabilizer_equals_h5": stabilizer == h5,
            "kernel_order": len(kernel),
            "faithful": kernel == frozenset({identity}),
            "orbit_stabilizer_identity": f"{len(group)}={len(orbit)}*{len(stabilizer)}",
            "permutation_character": [character[label]["fixed_cosets"] for label in _CLASS_ORDER],
            "canonical_action_sha256": _sha256_json(
                [list(actions[element]) for element in group]
            ),
        },
        "character_theory": character_certificate,
        "icosahedral_adjacency": adjacency,
        "exterior_one_generation_witness": exterior,
        "invariant_partition_antibridge": {
            "enumerated_subset_count": 1 << len(cosets),
            "invariant_subset_size_counts": {
                str(size): count for size, count in sorted(invariant_sizes.items())
            },
            "proper_nonempty_invariant_subset_exists": any(
                0 < size < len(cosets) for size in invariant_sizes
            ),
            "invariant_8_3_1_set_partition_exists": False,
            "no_go_verified": sorted(invariant_sizes) == [0, len(cosets)],
            "linear_module_8_3_1_decomposition_exists": character_certificate[
                "sm_adjoint_match"
            ]["both_opposite_triplet_assignments_match"],
            "distinction": (
                "Transitivity forbids a partition of the twelve named ports into separately "
                "A5-invariant subsets of sizes 8, 3, and 1. It does not forbid invariant "
                "linear subspaces: the real permutation module is isomorphic to an SM-adjoint "
                "restriction after choosing opposite triplet embeddings."
            ),
        },
    }


def _character_certificate(permutation_character: tuple[_Q5, ...]) -> dict[str, Any]:
    gram = {
        left: {
            right: _q5_text(_inner_product(row, _IRREPS[right]))
            for right in _IRREPS
        }
        for left, row in _IRREPS.items()
    }
    orthonormal = all(
        _inner_product(left_row, right_row) == (_ONE if left == right else _ZERO)
        for left, left_row in _IRREPS.items()
        for right, right_row in _IRREPS.items()
    )
    dimensions = {name: _q5_integer(row[0]) for name, row in _IRREPS.items()}
    decomposition = _decompose(permutation_character)

    assignments: list[dict[str, Any]] = []
    for color_triplet, weak_triplet in (("3prime", "3"), ("3", "3prime")):
        defining_character = _IRREPS[color_triplet]
        su3_adjoint = tuple(value * value - _ONE for value in defining_character)
        su3_decomposition = _decompose(su3_adjoint)
        sm_character = tuple(
            color + weak + hypercharge
            for color, weak, hypercharge in zip(
                su3_adjoint,
                _IRREPS[weak_triplet],
                _IRREPS["1"],
                strict=True,
            )
        )
        assignments.append(
            {
                "su3_defining_triplet": color_triplet,
                "su3_adjoint_character": [_q5_text(value) for value in su3_adjoint],
                "su3_adjoint_decomposition": su3_decomposition,
                "su3_adjoint_dimension": _q5_integer(su3_adjoint[0]),
                "su2_adjoint_triplet": weak_triplet,
                "u1_adjoint_irrep": "1",
                "sm_adjoint_character": [_q5_text(value) for value in sm_character],
                "sm_adjoint_dimension": _q5_integer(sm_character[0]),
                "matches_port_character": sm_character == permutation_character,
            }
        )

    return {
        "class_order": list(_CLASS_ORDER),
        "class_sizes": list(_CLASS_SIZES),
        "irreducible_character_table": {
            name: [_q5_text(value) for value in row] for name, row in _IRREPS.items()
        },
        "character_table_gram_matrix": gram,
        "character_table_checks": {
            "orthonormal": orthonormal,
            "dimension_squares_sum": sum(dimension * dimension for dimension in dimensions.values()),
            "dimension_squares_sum_to_60": sum(
                dimension * dimension for dimension in dimensions.values()
            )
            == 60,
        },
        "port_module": {
            "character": [_q5_text(value) for value in permutation_character],
            "decomposition": decomposition,
            "dimension": _q5_integer(permutation_character[0]),
            "multiplicity_free": all(multiplicity in (0, 1) for multiplicity in decomposition.values()),
        },
        "sm_adjoint_match": {
            "assignments": assignments,
            "both_opposite_triplet_assignments_match": all(
                assignment["matches_port_character"] for assignment in assignments
            ),
            "semantic_boundary": (
                "End_0(3) and End_0(3prime) give the eight-dimensional color adjoint "
                "characters 3+5 and 3prime+5 respectively. Pairing color with the opposite "
                "weak triplet and the trivial U(1) yields the port character. This is a "
                "linear-module character identity, not a pointwise 8+3+1 port partition."
            ),
        },
    }


def _adjacency_certificate(
    *,
    group: tuple[Permutation, ...],
    h5: frozenset[Permutation],
    actions: Mapping[Permutation, ActionRow],
    base_vertex: int,
    vertex_count: int,
) -> dict[str, Any]:
    unseen = set(range(vertex_count))
    stabilizer_orbits: list[tuple[int, ...]] = []
    while unseen:
        seed = min(unseen)
        orbit = tuple(sorted({actions[element][seed] for element in h5}))
        stabilizer_orbits.append(orbit)
        unseen.difference_update(orbit)

    candidates: list[tuple[tuple[tuple[int, int], ...], int]] = []
    for orbit in stabilizer_orbits:
        if len(orbit) != 5:
            continue
        seed = min(orbit)
        edges = {
            tuple(sorted((action[base_vertex], action[seed])))
            for action in actions.values()
        }
        edges.discard((base_vertex, base_vertex))
        candidates.append((tuple(sorted(edges)), seed))
    edges, orbital_seed = min(candidates)
    edge_set = frozenset(edges)
    adjacency = [[0 for _ in range(vertex_count)] for _ in range(vertex_count)]
    for left, right in edges:
        adjacency[left][right] = 1
        adjacency[right][left] = 1
    degrees = [sum(row) for row in adjacency]
    distance_profiles = [_distance_profile(adjacency, source) for source in range(vertex_count)]
    connected = all(sum(profile) == vertex_count and len(profile) == 4 for profile in distance_profiles)
    triangle_count = sum(
        1
        for i in range(vertex_count)
        for j in range(i + 1, vertex_count)
        for k in range(j + 1, vertex_count)
        if adjacency[i][j] and adjacency[i][k] and adjacency[j][k]
    )
    invariant = all(
        tuple(sorted((action[left], action[right]))) in edge_set
        for action in actions.values()
        for left, right in edges
    )

    characteristic_descending = _characteristic_polynomial(adjacency)
    expected_ascending = [1]
    for factor, multiplicity in (
        ([-5, 1], 1),
        ([1, 1], 5),
        ([-5, 0, 1], 3),
    ):
        for _ in range(multiplicity):
            expected_ascending = _poly_multiply(expected_ascending, factor)
    factorization_verified = list(reversed(characteristic_descending)) == expected_ascending

    identity = _identity_matrix(vertex_count)
    adjacency_squared = _matrix_multiply(adjacency, adjacency)
    annihilator = _matrix_multiply(
        _matrix_add(adjacency, _matrix_scale(identity, -5)),
        _matrix_multiply(
            _matrix_add(adjacency, identity),
            _matrix_add(adjacency_squared, _matrix_scale(identity, -5)),
        ),
    )
    annihilator_max_abs = max(abs(value) for row in annihilator for value in row)

    return {
        "construction": "lexicographically_first_valency_five_H5_suborbit_graph",
        "selected_orbital_seed_vertex": orbital_seed,
        "stabilizer_orbit_sizes": sorted(len(orbit) for orbit in stabilizer_orbits),
        "vertex_count": vertex_count,
        "edge_count": len(edges),
        "degree_sequence": sorted(degrees),
        "connected": connected,
        "a5_invariant": invariant,
        "triangle_count": triangle_count,
        "distance_layer_profiles": [
            list(profile) for profile in sorted({tuple(profile) for profile in distance_profiles})
        ],
        "canonical_edge_sha256": _sha256_json([list(edge) for edge in edges]),
        "characteristic_polynomial": {
            "variable": "x",
            "coefficients_descending": characteristic_descending,
            "factorization": "(x-5)(x+1)^5(x^2-5)^3",
            "factorization_verified": factorization_verified,
        },
        "minimal_polynomial": {
            "polynomial": "(x-5)(x+1)(x^2-5)",
            "annihilator_max_abs_entry": annihilator_max_abs,
            "annihilator_is_zero": annihilator_max_abs == 0,
        },
        "spectrum": {"5": 1, "sqrt(5)": 3, "-sqrt(5)": 3, "-1": 5},
        "canonical_rank_sequence": [1, 3, 3, 5],
        "triplet_assignment_boundary": (
            "The two square-root eigenvalues carry the Galois-conjugate triplets; "
            "exchanging the two valency-five orbitals exchanges their labeling."
        ),
    }


def _exterior_one_generation_witness() -> dict[str, Any]:
    """Exact conditional witness for Lambda^2(C+W) + Lambda^4(C+W)."""

    carrier_internal = {
        "C": {"dimension": 3, "su3": "3", "su2": "1", "hypercharge": Fraction(-1, 3)},
        "W": {"dimension": 2, "su3": "1", "su2": "2", "hypercharge": Fraction(1, 2)},
    }
    trace_balance = sum(
        row["dimension"] * row["hypercharge"] for row in carrier_internal.values()
    )
    field_specs = {
        "Q": {"c_power": 1, "w_power": 1, "su3": "3", "su2": "2"},
        "u_c": {"c_power": 2, "w_power": 0, "su3": "3bar", "su2": "1"},
        "e_c": {"c_power": 0, "w_power": 2, "su3": "1", "su2": "1"},
        "d_c": {"c_power": 2, "w_power": 2, "su3": "3bar", "su2": "1"},
        "L": {"c_power": 3, "w_power": 1, "su3": "1", "su2": "2"},
    }
    fields_internal: dict[str, dict[str, Any]] = {}
    for name, spec in field_specs.items():
        c_power = spec["c_power"]
        w_power = spec["w_power"]
        fields_internal[name] = {
            "origin": f"Lambda^{c_power}(C) tensor Lambda^{w_power}(W)",
            "exterior_degree": c_power + w_power,
            "c_power": c_power,
            "w_power": w_power,
            "su3": spec["su3"],
            "su2": spec["su2"],
            "hypercharge": c_power * carrier_internal["C"]["hypercharge"]
            + w_power * carrier_internal["W"]["hypercharge"],
            "dimension": math.comb(3, c_power) * math.comb(2, w_power),
        }

    expected_signatures = {
        "Q": ("3", "2", Fraction(1, 6), 6),
        "u_c": ("3bar", "1", Fraction(-2, 3), 3),
        "e_c": ("1", "1", Fraction(1), 1),
        "d_c": ("3bar", "1", Fraction(1, 3), 3),
        "L": ("1", "2", Fraction(-1, 2), 2),
    }
    field_signature_matches = all(
        (
            row["su3"],
            row["su2"],
            row["hypercharge"],
            row["dimension"],
        )
        == expected_signatures[name]
        for name, row in fields_internal.items()
    )
    conjugate_su3 = {"1": "1", "3": "3bar", "3bar": "3"}
    field_signatures = {
        (row["su3"], row["su2"], row["hypercharge"])
        for row in fields_internal.values()
    }
    dual_signatures = {
        (conjugate_su3[row["su3"]], row["su2"], -row["hypercharge"])
        for row in fields_internal.values()
    }
    dual_overlap = sorted(
        (
            {"su3": su3, "su2": su2, "hypercharge": _fraction_text(hypercharge)}
            for su3, su2, hypercharge in field_signatures.intersection(dual_signatures)
        ),
        key=lambda row: (row["su3"], row["su2"], row["hypercharge"]),
    )

    scalars = {
        "H": {"su3": "1", "su2": "2", "hypercharge": carrier_internal["W"]["hypercharge"]},
        "Hdag": {"su3": "1", "su2": "2", "hypercharge": -carrier_internal["W"]["hypercharge"]},
    }
    fermion_field_order = ("Q", "L", "u_c", "d_c", "e_c")
    scalar_order = ("H", "Hdag")
    canonical_channel_ids = (
        "Q_H_u_c",
        "Q_Hdag_d_c",
        "L_Hdag_e_c",
    )
    candidate_ledger: list[dict[str, Any]] = []
    for left, right in itertools.combinations_with_replacement(
        fermion_field_order,
        2,
    ):
        for scalar in scalar_order:
            factors = (
                fields_internal[left],
                scalars[scalar],
                fields_internal[right],
            )
            charge_sum = sum(
                (factor["hypercharge"] for factor in factors),
                Fraction(0),
            )
            su3_multiplicity = _singlet_multiplicity("su3", factors)
            su2_multiplicity = _singlet_multiplicity("su2", factors)
            multiplicity = int(
                charge_sum == 0
                and su3_multiplicity == 1
                and su2_multiplicity == 1
            )
            candidate_ledger.append(
                {
                    "candidate_id": f"{left}_{scalar}_{right}",
                    "fermion_pair": [left, right],
                    "scalar": scalar,
                    "factors": [left, scalar, right],
                    "hypercharge_sum": _fraction_text(charge_sum),
                    "su3_singlet_multiplicity": su3_multiplicity,
                    "su2_singlet_multiplicity": su2_multiplicity,
                    "invariant_line_multiplicity": multiplicity,
                }
            )
    candidate_by_id = {row["candidate_id"]: row for row in candidate_ledger}
    channels: dict[str, dict[str, Any]] = {}
    for name in canonical_channel_ids:
        candidate = candidate_by_id[name]
        channels[name] = {
            key: value
            for key, value in candidate.items()
            if key not in {"candidate_id", "fermion_pair", "scalar"}
        }
    nonzero_candidate_ids = [
        row["candidate_id"]
        for row in candidate_ledger
        if row["invariant_line_multiplicity"] > 0
    ]
    expected_candidate_count = (
        math.comb(len(fermion_field_order) + 1, 2) * len(scalar_order)
    )
    actual_candidate_specs = tuple(
        (
            row["candidate_id"],
            row["fermion_pair"][0],
            row["scalar"],
            row["fermion_pair"][1],
        )
        for row in candidate_ledger
    )
    exhaustive_candidate_ledger_complete = bool(
        len(candidate_ledger) == expected_candidate_count == 30
        and len(candidate_by_id) == len(candidate_ledger)
        and actual_candidate_specs == _EXPECTED_ONE_HIGGS_CANDIDATE_SPECS
    )
    nonzero_set_equals_canonical = set(nonzero_candidate_ids) == set(
        canonical_channel_ids
    )
    invariant_line_count = sum(
        row["invariant_line_multiplicity"] for row in candidate_ledger
    )
    nonabelian_singlet_candidate_count = sum(
        row["su3_singlet_multiplicity"] == 1
        and row["su2_singlet_multiplicity"] == 1
        for row in candidate_ledger
    )
    hypercharge_neutral_candidate_count = sum(
        row["hypercharge_sum"] == "0" for row in candidate_ledger
    )

    color_dimension = {"1": 1, "3": 3, "3bar": 3}
    weak_dimension = {"1": 1, "2": 2}
    color_cubic_index = {"1": Fraction(0), "3": Fraction(1), "3bar": Fraction(-1)}
    color_quadratic_index = {"1": Fraction(0), "3": Fraction(1, 2), "3bar": Fraction(1, 2)}
    weak_quadratic_index = {"1": Fraction(0), "2": Fraction(1, 2)}
    anomalies = {
        "SU3_cubed": sum(
            weak_dimension[row["su2"]] * color_cubic_index[row["su3"]]
            for row in fields_internal.values()
        ),
        "SU3_squared_U1": sum(
            weak_dimension[row["su2"]]
            * color_quadratic_index[row["su3"]]
            * row["hypercharge"]
            for row in fields_internal.values()
        ),
        "SU2_squared_U1": sum(
            color_dimension[row["su3"]]
            * weak_quadratic_index[row["su2"]]
            * row["hypercharge"]
            for row in fields_internal.values()
        ),
        "gravity_squared_U1": sum(
            row["dimension"] * row["hypercharge"] for row in fields_internal.values()
        ),
        "U1_cubed": sum(
            row["dimension"] * row["hypercharge"] ** 3 for row in fields_internal.values()
        ),
    }
    weak_doublet_count = sum(
        color_dimension[row["su3"]]
        for row in fields_internal.values()
        if row["su2"] == "2"
    )

    fields = {
        name: {
            **{key: value for key, value in row.items() if key != "hypercharge"},
            "hypercharge": _fraction_text(row["hypercharge"]),
        }
        for name, row in fields_internal.items()
    }
    carrier = {
        name: {
            **{key: value for key, value in row.items() if key != "hypercharge"},
            "hypercharge": _fraction_text(row["hypercharge"]),
        }
        for name, row in carrier_internal.items()
    }
    return {
        "status": "conditional_exact_representation_witness",
        "carrier": {
            "definition": "V=C+W",
            "summands": carrier,
            "trace_balance": _fraction_text(trace_balance),
            "trace_balanced": trace_balance == 0,
        },
        "matter_package": {
            "definition": "M1=Lambda^2(C+W)+Lambda^4(C+W)",
            "lambda2_dimension": math.comb(5, 2),
            "lambda4_dimension": math.comb(5, 4),
            "complex_dimension": math.comb(5, 2) + math.comb(5, 4),
            "fields": fields,
            "field_dimension_sum": sum(row["dimension"] for row in fields_internal.values()),
            "field_signature_matches_canonical_generation": field_signature_matches,
            "dual_overlap": dual_overlap,
            "disjoint_from_dual": not dual_overlap,
        },
        "one_higgs_invariant_lines": {
            "conditional_scalar_identification": "H=W",
            "fermion_pair_enumeration": "unordered_with_repetition",
            "fermion_field_order": list(fermion_field_order),
            "scalar_order": list(scalar_order),
            "fermion_pair_count": math.comb(len(fermion_field_order) + 1, 2),
            "candidate_count": len(candidate_ledger),
            "candidate_ids": [row["candidate_id"] for row in candidate_ledger],
            "candidate_ledger": candidate_ledger,
            "nonabelian_singlet_candidate_count": (
                nonabelian_singlet_candidate_count
            ),
            "hypercharge_neutral_candidate_count": (
                hypercharge_neutral_candidate_count
            ),
            "channels": channels,
            "canonical_invariant_candidate_ids": list(canonical_channel_ids),
            "nonzero_invariant_candidate_ids": nonzero_candidate_ids,
            "exhaustive_candidate_ledger_complete": (
                exhaustive_candidate_ledger_complete
            ),
            "nonzero_set_equals_three_canonical_channels": (
                nonzero_set_equals_canonical
            ),
            "invariant_line_count": invariant_line_count,
            "all_channels_exact_singlets": bool(
                invariant_line_count == len(channels) == 3
                and all(
                    row["invariant_line_multiplicity"] == 1
                    for row in channels.values()
                )
            ),
        },
        "anomalies": {
            "convention": "left_handed_Weyl",
            "coefficients": {name: _fraction_text(value) for name, value in anomalies.items()},
            "all_coefficients_zero": all(value == 0 for value in anomalies.values()),
            "su2_witten_doublet_count": weak_doublet_count,
            "su2_witten_parity_mod_2": weak_doublet_count % 2,
            "su2_witten_parity_even": weak_doublet_count % 2 == 0,
        },
        "weak_doublet_count": {
            "colored_Q_copies": color_dimension[fields_internal["Q"]["su3"]],
            "lepton_L_copies": color_dimension[fields_internal["L"]["su3"]],
            "multiplicity_per_generation": weak_doublet_count,
            "conditional_beta_EW": weak_doublet_count,
            "physical_load_identification": False,
        },
        "selection_boundary": {
            "full_even_exterior_dimension_including_lambda0": 1
            + math.comb(5, 2)
            + math.comb(5, 4),
            "vacuum_singlet_excluded_by_computation": False,
            "other_anomaly_free_light_sectors_excluded": False,
            "package_selected_by_native_dynamics": False,
            "H_equals_W_selected_by_native_dynamics": False,
            "physical_family_attachment": False,
            "continuum_spin_qft_realized": False,
            "claim_boundary": (
                "The exact exterior witness is conditional on V=C+W, selection of the "
                "non-vacuum even package, and H=W. Its exhaustive 30-row representation ledger "
                "proves exactly three one-Higgs gauge-invariant lines at this representation "
                "level. It neither excludes Lambda^0 or other "
                "anomaly-free sectors nor supplies family descent, continuum QFT, or a "
                "physical port-load map."
            ),
        },
    }


def _singlet_multiplicity(group: str, factors: tuple[Mapping[str, Any], ...]) -> int:
    representations = [str(row[group]) for row in factors if row[group] != "1"]
    if not representations:
        return 1
    if group == "su3" and sorted(representations) == ["3", "3bar"]:
        return 1
    if group == "su2" and representations == ["2", "2"]:
        return 1
    return 0


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _physical_promotion(
    *,
    structural_receipt: bool,
    sm_match: bool,
    partition_no_go: bool,
    exterior_witness: bool,
) -> dict[str, Any]:
    gates = {
        "FINITE_A5_STRUCTURAL_RECEIPT": bool(structural_receipt),
        "SM_ADJOINT_MODULE_CHARACTER_MATCH": bool(sm_match),
        "INVARIANT_SET_8_3_1_ANTIBRIDGE": bool(partition_no_go),
        "CONDITIONAL_EXTERIOR_ONE_GENERATION_WITNESS": bool(exterior_witness),
        **{name: False for name in _PHYSICAL_GATES},
    }
    blockers = [name for name, accepted in gates.items() if not accepted]
    return {
        "strongest_allowed_claim": "exact_finite_structural_certificate_only",
        "gate_status": gates,
        "open_gate_descriptions": dict(_PHYSICAL_GATES),
        "blockers": blockers,
        "promotion_allowed": False,
        "physical_port_current_identification": False,
        "physical_standard_model_derivation": False,
        "reason": (
            "Exact dimension and character agreement does not supply the source-selected port "
            "orbit, a physical current algebra, a refinement intertwiner, global-form descent, "
            "or selection of the conditional exterior matter package, realized family descent, "
            "continuum QFT, weak-load map, and exclusion of extra light sectors."
        ),
    }


def _normalize_imported_evidence(
    evidence: Sequence[Mapping[str, Any]] | None,
) -> dict[str, Any]:
    if evidence is None:
        raw_rows: list[Any] = []
    elif isinstance(evidence, Sequence) and not isinstance(evidence, (str, bytes)):
        raw_rows = list(evidence)
    else:
        raw_rows = [evidence]
    rows: list[dict[str, Any]] = []
    for raw in raw_rows:
        item = dict(raw) if isinstance(raw, Mapping) else {}
        source = str(item.get("source") or "")
        digest = str(item.get("content_sha256") or "").lower()
        theorem_values = item.get("theorems")
        theorems = (
            [str(value) for value in theorem_values]
            if isinstance(theorem_values, list)
            else []
        )
        hash_valid = bool(re.fullmatch(r"(?:sha256:)?[0-9a-f]{64}", digest))
        rows.append(
            {
                "source": source,
                "source_commit": str(item.get("source_commit") or "") or None,
                "content_sha256": digest or None,
                "theorems": theorems,
                "declared_verification_status": str(
                    item.get("verification_status") or "unverified_import"
                ),
                "hash_format_valid": hash_valid,
                "well_formed": bool(source and theorems and hash_valid),
                "evidence_class": "imported_theorem_provenance_only",
                "simulation_receipt_eligible": False,
                "used_for_native_receipts": False,
            }
        )
    return {
        "provided_count": len(rows),
        "records": rows,
        "all_well_formed": bool(rows) and all(row["well_formed"] for row in rows),
        "simulation_receipts_promoted_by_import": False,
        "used_for_native_computation": False,
    }


def _imports_are_nonpromoting(report: Mapping[str, Any]) -> bool:
    evidence_layers = report.get("evidence_layers")
    if not isinstance(evidence_layers, Mapping):
        return False
    imported = evidence_layers.get("imported_theorem_evidence")
    if not isinstance(imported, Mapping):
        return False
    rows = imported.get("records")
    if not isinstance(rows, list):
        return False
    return (
        imported.get("simulation_receipts_promoted_by_import") is False
        and imported.get("used_for_native_computation") is False
        and all(
            isinstance(row, Mapping)
            and row.get("simulation_receipt_eligible") is False
            and row.get("used_for_native_receipts") is False
            for row in rows
        )
    )


def _import_summary_is_consistent(report: Mapping[str, Any]) -> bool:
    evidence_layers = report.get("evidence_layers")
    if not isinstance(evidence_layers, Mapping):
        return False
    imported = evidence_layers.get("imported_theorem_evidence")
    if not isinstance(imported, Mapping):
        return False
    rows = imported.get("records")
    if not isinstance(rows, list) or not all(isinstance(row, Mapping) for row in rows):
        return False
    expected_all_well_formed = bool(rows) and all(row.get("well_formed") is True for row in rows)
    return (
        imported.get("provided_count") == len(rows)
        and imported.get("all_well_formed") is expected_all_well_formed
    )


def _label_conjugacy_classes(
    group: tuple[Permutation, ...],
    inverses: Mapping[Permutation, Permutation],
) -> dict[str, frozenset[Permutation]]:
    unseen = set(group)
    classes: list[frozenset[Permutation]] = []
    while unseen:
        element = min(unseen)
        conjugates = frozenset(
            _compose(_compose(conjugator, element), inverses[conjugator])
            for conjugator in group
        )
        classes.append(conjugates)
        unseen.difference_update(conjugates)
    by_order: dict[int, list[frozenset[Permutation]]] = {}
    for conjugacy_class in classes:
        order = _permutation_order(min(conjugacy_class))
        by_order.setdefault(order, []).append(conjugacy_class)
    five_classes = sorted(by_order[5], key=lambda row: tuple(sorted(row)))
    return {
        "1A": by_order[1][0],
        "2A": by_order[2][0],
        "3A": by_order[3][0],
        "5A": five_classes[0],
        "5B": five_classes[1],
    }


def _class_fixed_point_character(
    elements: frozenset[Permutation],
    actions: Mapping[Permutation, ActionRow],
) -> dict[str, Any]:
    values = [
        sum(index == image for index, image in enumerate(actions[element]))
        for element in elements
    ]
    return {
        "fixed_cosets": values[0],
        "constant_on_class": len(set(values)) == 1,
    }


def _left_cosets(
    group: tuple[Permutation, ...],
    subgroup: frozenset[Permutation],
) -> tuple[frozenset[Permutation], ...]:
    remaining = set(group)
    cosets: list[frozenset[Permutation]] = []
    while remaining:
        representative = min(remaining)
        coset = frozenset(_compose(representative, element) for element in subgroup)
        cosets.append(coset)
        remaining.difference_update(coset)
    return tuple(sorted(cosets, key=lambda coset: tuple(sorted(coset))))


def _invariant_subset_size_counts(
    actions: tuple[ActionRow, ...],
    vertex_count: int,
) -> Counter[int]:
    counts: Counter[int] = Counter()
    for mask in range(1 << vertex_count):
        invariant = True
        for action in actions:
            image_mask = 0
            for vertex in range(vertex_count):
                if mask & (1 << vertex):
                    image_mask |= 1 << action[vertex]
            if image_mask != mask:
                invariant = False
                break
        if invariant:
            counts[mask.bit_count()] += 1
    return counts


def _inner_product(left: tuple[_Q5, ...], right: tuple[_Q5, ...]) -> _Q5:
    total = sum(
        (size * lhs * rhs for size, lhs, rhs in zip(_CLASS_SIZES, left, right, strict=True)),
        _ZERO,
    )
    return total / 60


def _decompose(character: tuple[_Q5, ...]) -> dict[str, int]:
    return {
        name: _q5_integer(_inner_product(character, row))
        for name, row in _IRREPS.items()
    }


def _as_q5(value: _Q5 | int | Fraction) -> _Q5:
    return value if isinstance(value, _Q5) else _Q5(Fraction(value))


def _q5_integer(value: _Q5) -> int:
    if value.sqrt5 != 0 or value.rational.denominator != 1:
        raise ArithmeticError(f"expected an exact integer in Q(sqrt(5)), got {value}")
    return int(value.rational)


def _q5_text(value: _Q5) -> str:
    if value == _PHI:
        return "phi"
    if value == _PHI_BAR:
        return "phi_bar"
    if value.sqrt5 == 0:
        return str(value.rational.numerator) if value.rational.denominator == 1 else str(value.rational)
    return f"({value.rational})+({value.sqrt5})*sqrt(5)"


def _is_even(permutation: tuple[int, ...]) -> bool:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return inversions % 2 == 0


def _compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[index]] for index in range(5))  # type: ignore[return-value]


def _compose_action(left: ActionRow, right: ActionRow) -> ActionRow:
    return tuple(left[right[index]] for index in range(len(right)))


def _inverse(permutation: Permutation) -> Permutation:
    result = [0] * 5
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)  # type: ignore[return-value]


def _power(permutation: Permutation, exponent: int) -> Permutation:
    result: Permutation = (0, 1, 2, 3, 4)
    for _ in range(exponent):
        result = _compose(result, permutation)
    return result


def _permutation_order(permutation: Permutation) -> int:
    identity: Permutation = (0, 1, 2, 3, 4)
    value = identity
    for order in range(1, 61):
        value = _compose(value, permutation)
        if value == identity:
            return order
    raise ArithmeticError("permutation order exceeded A5 group order")


def _distance_profile(adjacency: list[list[int]], source: int) -> list[int]:
    distances = {source: 0}
    queue: deque[int] = deque([source])
    while queue:
        vertex = queue.popleft()
        for neighbor, adjacent in enumerate(adjacency[vertex]):
            if adjacent and neighbor not in distances:
                distances[neighbor] = distances[vertex] + 1
                queue.append(neighbor)
    if not distances:
        return []
    return [sum(distance == layer for distance in distances.values()) for layer in range(max(distances.values()) + 1)]


def _characteristic_polynomial(matrix: list[list[int]]) -> list[int]:
    dimension = len(matrix)
    power = _identity_matrix(dimension)
    traces: list[int] = []
    for _ in range(dimension):
        power = _matrix_multiply(power, matrix)
        traces.append(sum(power[index][index] for index in range(dimension)))
    coefficients = [Fraction(1)]
    for degree in range(1, dimension + 1):
        coefficient = -sum(
            coefficients[degree - power_degree] * traces[power_degree - 1]
            for power_degree in range(1, degree + 1)
        ) / degree
        coefficients.append(coefficient)
    if any(value.denominator != 1 for value in coefficients):
        raise ArithmeticError("integer matrix produced nonintegral characteristic coefficients")
    return [int(value) for value in coefficients]


def _identity_matrix(dimension: int) -> list[list[int]]:
    return [[int(row == column) for column in range(dimension)] for row in range(dimension)]


def _matrix_multiply(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    rows = len(left)
    columns = len(right[0])
    inner = len(right)
    return [
        [sum(left[row][index] * right[index][column] for index in range(inner)) for column in range(columns)]
        for row in range(rows)
    ]


def _matrix_add(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    return [
        [lhs + rhs for lhs, rhs in zip(left_row, right_row, strict=True)]
        for left_row, right_row in zip(left, right, strict=True)
    ]


def _matrix_scale(matrix: list[list[int]], scalar: int) -> list[list[int]]:
    return [[scalar * value for value in row] for row in matrix]


def _poly_multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            result[left_degree + right_degree] += left_value * right_value
    return result


def _payload_sha256(report: Mapping[str, Any]) -> str:
    payload = dict(report)
    payload.pop("certificate_payload_sha256", None)
    return _sha256_json(payload)


def _sha256_json(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()
