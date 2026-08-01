#!/usr/bin/env python3
"""Independent exact verifier for the issue #566 stage-two Jacobi reduction.

No producer module is imported.  This verifier rebuilds the general bracket,
forms Jacobi by multiplying linear coefficient forms, reconstructs the signed
target orbits, and performs separate exact elimination over Q and Q(sqrt(5)).
"""

from __future__ import annotations

import ast
import copy
import hashlib
import itertools
import json
import math
import runpy
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
STAGE1 = REPO / "code/a5_closure/issue_566_bracket_space_stage1"
STAGE1_BASIS = STAGE1 / "a5_alternating_bracket_reynolds_basis.json"
STAGE1_RECEIPT = STAGE1 / "a5_alternating_bracket_space_stage1.receipt.json"
SYSTEM_PATH = HERE / "a5_jacobi_system_reduction.json"
RECEIPT_PATH = HERE / "a5_jacobi_stage2.receipt.json"
PRODUCER_PATH = HERE / "certify.py"
TEST_PATH = HERE / "test_stage2.py"

N = 12
P = 14
MONOMIALS = [(a, b) for a in range(P) for b in range(a, P)]
MONOMIAL_INDEX = {item: i for i, item in enumerate(MONOMIALS)}
EXACT_IMPORT_ROOTS_BY_FILE = {
    "certify.py": {
        "__future__",
        "argparse",
        "ast",
        "copy",
        "dataclasses",
        "fractions",
        "hashlib",
        "itertools",
        "json",
        "math",
        "pathlib",
        "runpy",
        "sys",
        "typing",
    },
    "verify.py": {
        "__future__",
        "ast",
        "copy",
        "dataclasses",
        "fractions",
        "hashlib",
        "itertools",
        "json",
        "math",
        "pathlib",
        "runpy",
        "sys",
        "typing",
    },
    "test_stage2.py": {
        "__future__",
        "json",
        "pathlib",
        "subprocess",
        "sys",
        "unittest",
    },
}
RUNTIME_FORBIDDEN_IMPORT_ROOTS = {"certify.py": {"subprocess"}, "verify.py": {"subprocess"}}

STAGE1_BASIS_KEYS = {
    "ambient_coordinate_count", "basis", "coordinate_convention", "domain_dimension", "field",
    "group_action_convention", "issue", "normalization", "port_dimension", "proper_port_group_sha256", "schema",
}
STAGE1_BASIS_ROW_KEYS = {"basis_id", "entries", "integral_orbit_scale", "orbit_size", "seed_coordinate", "stabilizer_size"}
STAGE1_RECEIPT_KEYS = {
    "claim_boundary", "dimension_certificate", "implementation_pins", "issue", "later_gates", "mutation_tests",
    "proper_port_action", "receipt_sha256", "representation", "reynolds_basis", "schema", "semantic_input", "status", "target_firewall",
}
STAGE1_SEMANTIC_KEYS = {"accessed_value_paths", "carrier_projection_sha256", "file_sha256", "ignored_value_sections", "path"}
STAGE1_FIREWALL_KEYS = {
    "arbitrary_input_path_allowed", "audited_import_roots", "audited_import_roots_by_file", "basis_seed_rule",
    "command_line_target_parameters_allowed", "desired_coefficient_input_allowed", "enabled", "environment_target_read",
    "exact_manifest_keysets_enforced", "measurement_input_allowed", "network_target_read", "nonstdlib_imports",
    "semantic_input_count", "target_free_scope",
}
STAGE1_ACTION_KEYS = {
    "a5_isomorphism", "all_rows_bijective", "closure_checked_pairs", "identity_present", "incidence_automorphism_count",
    "incidence_preserved", "inverses_present", "orientation_preserving_count", "orientation_reversing_count",
    "oriented_faces_preserved", "permutation_rows", "permutation_rows_sha256", "port_count", "transitive_on_ports",
}
STAGE1_A5_KEYS = {
    "bijection_size", "generator_orders", "homomorphism_checked_pairs", "isomorphic", "mapping_sha256", "model",
    "model_generators", "model_order", "port_generators",
}
STAGE1_REPRESENTATION_KEYS = {"V_dimension", "action", "exterior_square_V_dimension", "field", "hom_ambient_dimension"}
STAGE1_DIMENSION_KEYS = {"class_rows", "dimension", "formula", "projector_trace_denominator", "projector_trace_numerator"}
STAGE1_CLASS_ROW_KEYS = {
    "chi_V", "chi_V_of_square", "chi_exterior_square_V", "count", "element_order",
    "projector_trace_contribution_each", "projector_trace_contribution_total",
}
STAGE1_REYNOLDS_KEYS = {
    "all_signed_orbits_orientable", "alternation_checks", "basis_count", "canonical_json_sha256", "covariance_checks",
    "exact_alternation_passed", "exact_covariance_passed", "normalization", "orbit_sizes", "path", "rank_certificate",
    "support_union_size", "supports_pairwise_disjoint",
}
STAGE1_RANK_KEYS = {
    "ambient_coordinate_count", "gram_determinant", "gram_diagonal", "matrix_shape", "pivot_coordinates",
    "pivot_minor_determinant", "rank_over_Q",
}
STAGE1_FRACTION_KEYS = {"denominator", "numerator"}
STAGE1_PIN_KEYS = {"independent_verifier_sha256", "producer_sha256", "test_sha256"}
MUTATION_ROW_KEYS = {"actual_error", "expected_error", "name", "passed"}

SYSTEM_KEYS = {"field", "fixed_line_reduction", "issue", "jacobi_tensor", "monomial_order", "parameter_variables", "schema", "upstream"}
UPSTREAM_KEYS = {
    "basis_canonical_json_sha256", "basis_file_sha256", "basis_path", "proper_group_sha256",
    "receipt_file_sha256", "receipt_path", "receipt_self_sha256",
}
JACOBI_KEYS = {
    "all_coordinate_equations_nonzero", "all_coordinate_equations_sha256", "convention", "coordinate_equation_count",
    "reduced_integer_equation_count", "reduced_integer_equations", "representative_coefficient_rank_over_Q",
    "representative_equations", "representative_rref_pivot_monomial_indices", "signed_target_orbit_count",
    "signed_target_orbit_sizes", "target_invariant_character_certificate", "target_space_dimension", "zero_locus_equivalence",
}
TARGET_CHARACTER_KEYS = {"class_rows", "dimension", "projector_trace_denominator", "projector_trace_numerator"}
TARGET_CLASS_ROW_KEYS = {
    "chi_V", "chi_V_of_cube", "chi_V_of_square", "chi_exterior_cube_V", "contribution_each",
    "contribution_total", "count", "element_order",
}
REPRESENTATIVE_KEYS = {"equation", "signed_orbit_size", "target_coordinate"}
FIXED_REDUCTION_KEYS = {
    "central_specialization", "channel_decomposition", "contracted_coefficient_rank_over_Q", "contracted_coordinate_count",
    "contracted_integer_equations_in_x", "contracted_rref_pivot_monomial_indices", "derivation_only_subspace",
    "derivation_product_equations_in_channel_variables", "derivation_weight_arrangement", "fixed_vector", "full_rank_split",
    "product_equation_rank_over_Qsqrt5", "product_rowspace_equals_contracted_jacobi_rowspace",
    "residual_equations_in_channel_variables", "residual_integer_equations_in_x", "residual_rank_over_Qsqrt5",
}
CHANNEL_DECOMPOSITION_KEYS = {"channel_variables", "inverse_transform_rows", "spectral_decomposition", "splitting_field", "transform_determinant", "transform_rows"}
D_CHANNEL_KEYS = {"id", "linear_form_in_x", "meaning"}
T_CHANNEL_KEYS = {"derivation_weight_in_d", "domain_sectors", "extraction_coordinate", "id", "linear_form_in_x", "output_sector"}
SPECTRAL_KEYS = {"canonical_valency_five_orbital", "minimal_polynomial", "operator_sha256", "sectors"}
SECTOR_KEYS = {"dimension", "id", "operator_eigenvalue"}
PRODUCT_DECLARATION_KEYS = {"channel_id", "equation", "weight_coefficients_on_d_plus_d_minus_d_five"}
DERIVATION_ONLY_KEYS = {"all_jacobi_equations_vanish_identically", "definition", "dimension", "scope"}
CENTRAL_SPECIALIZATION_KEYS = {"definition", "remaining_channel_variable_count", "residual_coefficient_rank_over_Qsqrt5", "solution_set_classified"}
ARRANGEMENT_KEYS = {"distinct_weight_hyperplane_count", "flat_count", "flat_dimension_counts", "relative_open_flats", "scope"}
FLAT_KEYS = {
    "channels_forced_zero_on_relative_open_flat", "channels_not_forced_zero_by_fixed_line_sector", "defining_rowspace_rref",
    "dimension", "flat_id", "scope", "weights_vanishing_identically",
}
RECEIPT_KEYS = {
    "attained", "claim_boundary", "implementation_pins", "issue", "mutation_tests", "not_attained", "receipt_sha256",
    "schema", "semantic_inputs", "status", "system_artifact", "target_firewall",
}
SEMANTIC_INPUT_KEYS = {"count", "paths", "stage1_later_gates_all_false"}
FIREWALL_KEYS = {
    "accepted_input_schemas", "arbitrary_input_path_allowed", "audited_import_roots_by_file", "coefficient_target_input_allowed",
    "enabled", "environment_target_read", "exact_import_allowlists_by_file", "exact_input_keysets_enforced",
    "measurement_input_allowed", "network_target_read", "nonstdlib_imports", "runtime_forbidden_import_roots_by_file",
}
ATTAINED_KEYS = {
    "derivation_only_three_parameter_family_certified", "derivation_weight_relative_open_flat_count",
    "equivariant_target_orbit_count", "exact_channel_coordinate_change", "fixed_line_product_quadrics",
    "general_parameter_count", "generic_weight_stratum_forces_all_non_derivation_channels_zero",
    "independent_jacobi_quadrics", "raw_jacobi_coordinate_count", "residual_independent_quadrics",
}
NOT_ATTAINED_KEYS = {
    "compactness_established", "full_jacobi_solution_variety_classified", "issue_566_closed", "preferred_bracket_selected",
    "rational_descent_components_classified", "source_selection", "special_weight_flat_residual_solutions_classified",
}
SYSTEM_ARTIFACT_KEYS = {"canonical_json_sha256", "path"}
IMPLEMENTATION_PIN_KEYS = {"independent_verifier_sha256", "producer_sha256", "test_sha256"}

STATUS = "EXACT_JACOBI_SYSTEM_AND_FIXED_LINE_REDUCTION__FULL_CLASSIFICATION_OPEN"
CLAIM_BOUNDARY = (
    "This packet constructs the complete Jacobi polynomial system for the stage-one fourteen-parameter equivariant alternating bracket and certifies an exact coefficient-row rank 38 system with an 11+27 rowspace decomposition. "
    "The rank-27 residual equation system remains applicable on every derivation-weight flat. The three-dimensional derivation-only family is a Jacobi statement only and is neutral about compactness and source selection. "
    "It does not select coefficients, classify the residual solution variety, establish compactness, identify a preferred Lie algebra or source, or close issue #566."
)
ARRANGEMENT_SCOPE = "product-sector forced-zero stratification only; the rank-27 residual equations still apply on every flat"
DERIVATION_ONLY_SCOPE = "Jacobi vanishing only; compactness and source selection are not established"

Permutation = tuple[int, ...]
Coordinate = tuple[int, int, int, int]


class VerificationError(RuntimeError):
    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


def check(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise VerificationError(code, message)


def check_keys(value: Any, expected: set[str], code: str) -> None:
    check(isinstance(value, Mapping), code, "value is not an object")
    check(set(value) == expected, code, f"keyset differs: {sorted(set(value) ^ expected)}")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def digest(value: Any) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_digest(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[i]] for i in range(len(left)))


def sign_of_sort(values: Sequence[int]) -> int:
    inversions = sum(values[i] > values[j] for i in range(len(values)) for j in range(i + 1, len(values)))
    return -1 if inversions % 2 else 1


def invoke_hardened_stage1_verifier() -> None:
    try:
        namespace = runpy.run_path(str(STAGE1 / "verify.py"), run_name="_oph_stage1_hardened_verify")
        result = namespace["verify"]()
    except Exception as exc:
        raise VerificationError("UPSTREAM_HARDENED_VERIFY", f"stage-one independent verifier failed: {exc}") from exc
    check(
        isinstance(result, Mapping)
        and result.get("verified") is True
        and result.get("basis_rank") == P
        and result.get("later_gates_all_false") is True,
        "UPSTREAM_HARDENED_VERIFY",
        "stage-one independent verifier returned an unexpected result",
    )


def validate_stage1_keysets(basis: Any, receipt: Any) -> None:
    check_keys(basis, STAGE1_BASIS_KEYS, "UPSTREAM_BASIS_KEYS")
    rows = basis.get("basis")
    check(isinstance(rows, list), "UPSTREAM_BASIS_ROWS", "basis rows are absent")
    for row in rows:
        check_keys(row, STAGE1_BASIS_ROW_KEYS, "UPSTREAM_BASIS_ROW_KEYS")
    check_keys(receipt, STAGE1_RECEIPT_KEYS, "UPSTREAM_RECEIPT_KEYS")
    check_keys(receipt.get("semantic_input"), STAGE1_SEMANTIC_KEYS, "UPSTREAM_SEMANTIC_KEYS")
    firewall = receipt.get("target_firewall")
    check_keys(firewall, STAGE1_FIREWALL_KEYS, "UPSTREAM_FIREWALL_KEYS")
    check_keys(firewall.get("audited_import_roots_by_file"), {"certify.py", "test_stage1.py", "verify.py"}, "UPSTREAM_IMPORT_FILE_KEYS")
    action = receipt.get("proper_port_action")
    check_keys(action, STAGE1_ACTION_KEYS, "UPSTREAM_ACTION_KEYS")
    check_keys(action.get("a5_isomorphism"), STAGE1_A5_KEYS, "UPSTREAM_A5_KEYS")
    check_keys(receipt.get("representation"), STAGE1_REPRESENTATION_KEYS, "UPSTREAM_REPRESENTATION_KEYS")
    dimension = receipt.get("dimension_certificate")
    check_keys(dimension, STAGE1_DIMENSION_KEYS, "UPSTREAM_DIMENSION_KEYS")
    check(isinstance(dimension.get("class_rows"), list), "UPSTREAM_CLASS_ROWS", "stage-one class rows are absent")
    for row in dimension["class_rows"]:
        check_keys(row, STAGE1_CLASS_ROW_KEYS, "UPSTREAM_CLASS_ROW_KEYS")
    reynolds = receipt.get("reynolds_basis")
    check_keys(reynolds, STAGE1_REYNOLDS_KEYS, "UPSTREAM_REYNOLDS_KEYS")
    rank = reynolds.get("rank_certificate")
    check_keys(rank, STAGE1_RANK_KEYS, "UPSTREAM_RANK_KEYS")
    check_keys(rank.get("pivot_minor_determinant"), STAGE1_FRACTION_KEYS, "UPSTREAM_FRACTION_KEYS")
    check_keys(rank.get("gram_determinant"), STAGE1_FRACTION_KEYS, "UPSTREAM_FRACTION_KEYS")
    check(isinstance(rank.get("gram_diagonal"), list), "UPSTREAM_GRAM_ROWS", "stage-one Gram rows are absent")
    for row in rank["gram_diagonal"]:
        check_keys(row, STAGE1_FRACTION_KEYS, "UPSTREAM_FRACTION_KEYS")
    later_keys = {
        "bracket_selected", "carrier_choice_derived_by_this_packet", "compactness_established", "issue_566_closed",
        "jacobi_solution_classification_complete", "jacobi_solution_found", "physical_current_identified", "source_selection",
    }
    check_keys(receipt.get("later_gates"), later_keys, "UPSTREAM_GATE_KEYS")
    check_keys(receipt.get("implementation_pins"), STAGE1_PIN_KEYS, "UPSTREAM_PIN_KEYS")
    mutations = receipt.get("mutation_tests")
    check(isinstance(mutations, list), "UPSTREAM_MUTATION_ROWS", "stage-one mutation rows are absent")
    for row in mutations:
        check_keys(row, MUTATION_ROW_KEYS, "UPSTREAM_MUTATION_KEYS")


def stage1_character_claim(group: Sequence[Permutation]) -> dict[str, Any]:
    bins: dict[tuple[int, int, int], int] = {}
    numerator = 0
    identity = tuple(range(N))
    for row in group:
        power = identity
        order = 0
        while True:
            order += 1
            power = compose(row, power)
            if power == identity:
                break
        square = compose(row, row)
        chi = sum(row[i] == i for i in range(N))
        chi2 = sum(square[i] == i for i in range(N))
        exterior2 = (chi * chi - chi2) // 2
        bins[(order, chi, chi2)] = bins.get((order, chi, chi2), 0) + 1
        numerator += chi * exterior2
    rows = []
    for (order, chi, chi2), count in sorted(bins.items()):
        exterior2 = (chi * chi - chi2) // 2
        rows.append(
            {
                "element_order": order,
                "count": count,
                "chi_V": chi,
                "chi_V_of_square": chi2,
                "chi_exterior_square_V": exterior2,
                "projector_trace_contribution_each": chi * exterior2,
                "projector_trace_contribution_total": count * chi * exterior2,
            }
        )
    return {
        "formula": "sum_g chi_V(g)*(chi_V(g)^2-chi_V(g^2))/2 divided by group_order",
        "class_rows": rows,
        "projector_trace_numerator": numerator,
        "projector_trace_denominator": len(group),
        "dimension": numerator // len(group),
    }


def validate_stage1_claims(
    basis: Mapping[str, Any], receipt: Mapping[str, Any], vectors: Sequence[Mapping[tuple[int, int, int], Fraction]], group: Sequence[Permutation]
) -> None:
    expected_metadata = {
        "schema": "oph.a5_alternating_bracket_reynolds_basis.v1",
        "issue": 566,
        "field": "Q",
        "port_dimension": 12,
        "domain_dimension": 66,
        "ambient_coordinate_count": 792,
        "coordinate_convention": "[output,left,right] with left<right; swapping inputs negates the coefficient",
        "group_action_convention": "(k,i,j) maps to (g(k),sort(g(i),g(j))) with the sorting sign",
        "normalization": "each vector is the exact Reynolds average (1/60)*sum_g g(seed)",
    }
    for key, expected in expected_metadata.items():
        check(basis.get(key) == expected, "UPSTREAM_BASIS_CLAIM", f"basis {key} differs")
    action = receipt["proper_port_action"]
    group_hash = digest([list(row) for row in group])
    check(basis.get("proper_port_group_sha256") == group_hash, "UPSTREAM_GROUP_HASH", "basis group pin differs")
    check(action.get("permutation_rows") == [list(row) for row in group], "UPSTREAM_GROUP_ROWS", "group serialization differs")
    check(action.get("permutation_rows_sha256") == group_hash, "UPSTREAM_GROUP_HASH", "group hash differs")
    expected_action = {
        "port_count": 12, "incidence_automorphism_count": 120, "orientation_preserving_count": 60,
        "orientation_reversing_count": 60, "identity_present": True, "all_rows_bijective": True,
        "incidence_preserved": True, "oriented_faces_preserved": True, "closure_checked_pairs": 3600,
        "inverses_present": True, "transitive_on_ports": True,
    }
    for key, expected in expected_action.items():
        check(action.get(key) == expected, "UPSTREAM_GROUP_CLAIM", f"action {key} differs")
    a5 = action["a5_isomorphism"]
    check(
        a5.get("model") == "even permutations on five unlabeled symbols" and a5.get("model_order") == 60
        and a5.get("bijection_size") == 60 and a5.get("generator_orders") == [2, 3, 5]
        and a5.get("homomorphism_checked_pairs") == 3600 and a5.get("isomorphic") is True,
        "UPSTREAM_A5_CLAIM", "A5 declaration differs",
    )
    check(all(tuple(row) in set(group) for row in a5.get("port_generators", [])), "UPSTREAM_A5_CLAIM", "port generator is absent")
    check(receipt.get("issue") == 566, "UPSTREAM_ISSUE", "wrong stage-one issue")
    check(receipt.get("status") == "EXACT_TARGET_FREE_SEARCH_SPACE_STAGE1_ONLY", "UPSTREAM_STATUS", "stage-one status differs")
    check(
        receipt.get("claim_boundary")
        == "Conditional on the pinned canonical oriented twelve-port carrier, this receipt certifies only the complete equivariant alternating-bracket search space on its permutation module. Target-free means that no desired gauge algebra, measurement, or coefficient target enters this stage; it does not mean that this packet derives the carrier choice. The receipt does not select a bracket, solve Jacobi, establish compactness, identify a physical current, or close issue #566.",
        "UPSTREAM_SCOPE", "stage-one scope differs",
    )
    check(
        receipt.get("representation")
        == {"field": "Q", "V_dimension": 12, "exterior_square_V_dimension": 66, "hom_ambient_dimension": 792,
            "action": "simultaneous output and signed unordered-input permutation"},
        "UPSTREAM_REPRESENTATION", "stage-one representation differs",
    )
    check(receipt.get("dimension_certificate") == stage1_character_claim(group), "UPSTREAM_DIMENSION", "stage-one character certificate differs")
    check(not any(receipt["later_gates"].values()), "UPSTREAM_GATE", "stage-one gate promoted")
    semantic = receipt["semantic_input"]
    check(
        semantic.get("path") == "code/a5_closure/manifests/echosahedral_federation_reference.json"
        and semantic.get("accessed_value_paths") == ["/schema", "/carrier/ports", "/carrier/edges", "/carrier/oriented_faces"]
        and semantic.get("ignored_value_sections")
        == ["/architecture", "/carrier/atoms_pairwise_orthogonal", "/carrier/atoms_sum_to_one",
            "/carrier/central_port_atoms", "/refinement_tower", "/source_readback"],
        "UPSTREAM_SEMANTIC_SCOPE", "stage-one semantic scope differs",
    )
    firewall = receipt["target_firewall"]
    check(
        firewall.get("enabled") is True and firewall.get("semantic_input_count") == 1
        and firewall.get("arbitrary_input_path_allowed") is False and firewall.get("command_line_target_parameters_allowed") is False
        and firewall.get("desired_coefficient_input_allowed") is False and firewall.get("measurement_input_allowed") is False
        and firewall.get("environment_target_read") is False and firewall.get("network_target_read") is False
        and firewall.get("nonstdlib_imports") == [] and firewall.get("exact_manifest_keysets_enforced") is True
        and firewall.get("basis_seed_rule") == "lexicographically first unassigned tensor coordinate"
        and firewall.get("target_free_scope") == "conditional_on_pinned_canonical_oriented_twelve_port_carrier",
        "UPSTREAM_FIREWALL", "stage-one firewall differs",
    )
    coordinates = sorted({coordinate for vector in vectors for coordinate in vector})
    supports = [set(vector) for vector in vectors]
    check(len(coordinates) == 792 and sum(map(len, supports)) == len(set().union(*supports)), "UPSTREAM_BASIS_SUPPORT", "support partition differs")
    basis_rank = q_rank([[vector.get(coordinate, Fraction(0)) for coordinate in coordinates] for vector in vectors])
    check(basis_rank == 14, "UPSTREAM_BASIS_RANK", "stage-one rank differs")
    reynolds = receipt["reynolds_basis"]
    check(
        reynolds.get("path") == "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
        and reynolds.get("canonical_json_sha256") == digest(basis) and reynolds.get("basis_count") == 14
        and reynolds.get("orbit_sizes") == sorted(map(len, vectors)) and reynolds.get("support_union_size") == 792
        and reynolds.get("supports_pairwise_disjoint") is True and reynolds.get("all_signed_orbits_orientable") is True
        and reynolds.get("normalization") == "R(seed)=(1/60)*sum_g g(seed); nonzero entries are signed reciprocals of orbit size"
        and reynolds.get("alternation_checks") == 24192 and reynolds.get("exact_alternation_passed") is True
        and reynolds.get("covariance_checks") == 665280 and reynolds.get("exact_covariance_passed") is True,
        "UPSTREAM_REYNOLDS_CLAIM", "stage-one Reynolds declaration differs",
    )
    rank = reynolds["rank_certificate"]
    gram = [sum((value * value for value in vector.values()), Fraction(0)) for vector in vectors]
    gram_product = math.prod(gram, start=Fraction(1))
    check(
        rank.get("ambient_coordinate_count") == 792 and rank.get("matrix_shape") == [792, 14]
        and rank.get("rank_over_Q") == basis_rank and rank.get("pivot_coordinates") == [list(min(vector)) for vector in vectors]
        and [Fraction(row["numerator"], row["denominator"]) for row in rank["gram_diagonal"]] == gram
        and Fraction(rank["gram_determinant"]["numerator"], rank["gram_determinant"]["denominator"]) == gram_product
        and Fraction(rank["pivot_minor_determinant"]["numerator"], rank["pivot_minor_determinant"]["denominator"]) == gram_product,
        "UPSTREAM_RANK_CLAIM", "stage-one rank certificate differs",
    )
    pins = receipt["implementation_pins"]
    check(
        pins.get("producer_sha256") == file_digest(STAGE1 / "certify.py")
        and pins.get("independent_verifier_sha256") == file_digest(STAGE1 / "verify.py")
        and pins.get("test_sha256") == file_digest(STAGE1 / "test_stage1.py"),
        "UPSTREAM_IMPLEMENTATION_PIN", "stage-one pin differs",
    )
    check(all(row.get("passed") is True for row in receipt["mutation_tests"]), "UPSTREAM_MUTATION", "stage-one mutation failed")


def decode_stage1() -> tuple[dict[str, Any], dict[str, Any], list[dict[tuple[int, int, int], Fraction]], list[Permutation]]:
    invoke_hardened_stage1_verifier()
    basis = json.loads(STAGE1_BASIS.read_text(encoding="utf-8"))
    receipt = json.loads(STAGE1_RECEIPT.read_text(encoding="utf-8"))
    validate_stage1_keysets(basis, receipt)
    check(basis.get("schema") == "oph.a5_alternating_bracket_reynolds_basis.v1", "UPSTREAM_SCHEMA", "bad basis schema")
    check(receipt.get("schema") == "oph.a5_alternating_bracket_space_stage1.receipt.v1", "UPSTREAM_SCHEMA", "bad receipt schema")
    unhashed = dict(receipt)
    stored = unhashed.pop("receipt_sha256", None)
    check(stored == digest(unhashed), "UPSTREAM_HASH", "stage-one receipt self-hash differs")
    check(digest(basis) == receipt["reynolds_basis"]["canonical_json_sha256"], "UPSTREAM_HASH", "stage-one basis hash differs")
    check(not any(receipt["later_gates"].values()), "UPSTREAM_GATE", "stage-one gate promoted")
    group = sorted(tuple(row) for row in receipt["proper_port_action"]["permutation_rows"])
    check(len(group) == len(set(group)) == 60, "UPSTREAM_GROUP", "wrong group row count")
    group_set = set(group)
    identity = tuple(range(N))
    check(identity in group_set, "UPSTREAM_GROUP", "identity is absent")
    for row in group:
        check(len(row) == N and set(row) == set(range(N)), "UPSTREAM_GROUP", "group row is not bijective")
        inverse = [0] * N
        for source, target in enumerate(row):
            inverse[target] = source
        check(tuple(inverse) in group_set, "UPSTREAM_GROUP", "group inverse is absent")
        for other in group:
            check(compose(row, other) in group_set, "UPSTREAM_GROUP", "group is not closed")
    check({row[0] for row in group} == set(range(N)), "UPSTREAM_GROUP", "group action is not transitive")
    vectors = []
    for number, row in enumerate(basis["basis"]):
        check(row["basis_id"] == f"R{number:02d}", "UPSTREAM_BASIS", "basis order differs")
        vector = {}
        for output, left, right, numerator, denominator in row["entries"]:
            coordinate = (output, left, right)
            check(0 <= output < N and 0 <= left < right < N and coordinate not in vector and denominator > 0, "UPSTREAM_BASIS", "bad basis entry")
            vector[coordinate] = Fraction(numerator, denominator)
        check(row["seed_coordinate"] == list(min(vector)), "UPSTREAM_BASIS", "basis seed differs")
        check(row["orbit_size"] == len(vector), "UPSTREAM_BASIS", "basis orbit size differs")
        check(row["stabilizer_size"] == 60 // len(vector), "UPSTREAM_BASIS", "basis stabilizer differs")
        check(row["integral_orbit_scale"] == len(vector), "UPSTREAM_BASIS", "basis scale differs")
        vectors.append(vector)
    check(len(vectors) == P and sum(map(len, vectors)) == 792, "UPSTREAM_BASIS", "basis shape differs")
    coordinates = [(output, left, right) for output in range(N) for left in range(N) for right in range(left + 1, N)]
    for vector in vectors:
        for row in group:
            for output, left, right in coordinates:
                a, b = row[left], row[right]
                image = (row[output], min(a, b), max(a, b))
                sign = 1 if a < b else -1
                check(vector.get(image, Fraction(0)) == sign * vector.get((output, left, right), Fraction(0)), "UPSTREAM_COVARIANCE", "stage-one covariance fails")
    validate_stage1_claims(basis, receipt, vectors, group)
    return basis, receipt, vectors, group


def coefficient_forms(vectors: Sequence[Mapping[tuple[int, int, int], Fraction]]) -> list[list[list[tuple[Fraction, ...]]]]:
    forms = [
        [[tuple(Fraction(0) for _ in range(P)) for _ in range(N)] for _ in range(N)]
        for _ in range(N)
    ]
    mutable = [[[[Fraction(0) for _ in range(P)] for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for parameter, vector in enumerate(vectors):
        for (output, left, right), value in vector.items():
            mutable[output][left][right][parameter] = value
            mutable[output][right][left][parameter] = -value
    for output in range(N):
        for left in range(N):
            for right in range(N):
                forms[output][left][right] = tuple(mutable[output][left][right])
    return forms


def multiply_forms(left: Sequence[Fraction], right: Sequence[Fraction]) -> tuple[Fraction, ...]:
    row = [Fraction(0)] * len(MONOMIALS)
    for a, value_a in enumerate(left):
        if not value_a:
            continue
        for b, value_b in enumerate(right):
            if value_b:
                row[MONOMIAL_INDEX[tuple(sorted((a, b)))]] += value_a * value_b
    return tuple(row)


def add_qrows(*rows: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(sum(values, Fraction(0)) for values in zip(*rows))


def jacobi_from_forms(forms: Sequence[Sequence[Sequence[Sequence[Fraction]]]]) -> dict[Coordinate, tuple[Fraction, ...]]:
    result = {}
    for i, j, k in itertools.combinations(range(N), 3):
        for output in range(N):
            total = tuple(Fraction(0) for _ in MONOMIALS)
            for middle in range(N):
                total = add_qrows(
                    total,
                    multiply_forms(forms[middle][i][j], forms[output][middle][k]),
                    multiply_forms(forms[middle][j][k], forms[output][middle][i]),
                    multiply_forms(forms[middle][k][i], forms[output][middle][j]),
                )
            check(any(total), "RAW_JACOBI", "unexpected zero raw component")
            result[(output, i, j, k)] = total
    check(len(result) == 2640, "RAW_JACOBI", "wrong raw component count")
    return result


def signed_orbits(group: Sequence[Permutation]) -> list[dict[Coordinate, int]]:
    remaining = {(o, i, j, k) for o in range(N) for i, j, k in itertools.combinations(range(N), 3)}
    result = []
    while remaining:
        seed = min(remaining)
        signs = {}
        for row in group:
            images = [row[seed[1]], row[seed[2]], row[seed[3]]]
            image = (row[seed[0]], *sorted(images))
            sign = sign_of_sort(images)
            check(image not in signs or signs[image] == sign, "TARGET_ORBIT", "signed stabilizer conflict")
            signs[image] = sign
        remaining -= set(signs)
        result.append(signs)
    check(len(result) == 44 and all(len(orbit) == 60 for orbit in result), "TARGET_ORBIT", "target orbit count differs")
    return result


def target_character_claim(group: Sequence[Permutation]) -> dict[str, Any]:
    bins: dict[tuple[int, int, int, int], int] = {}
    numerator = 0
    identity = tuple(range(N))
    for row in group:
        power = identity
        order = 0
        while True:
            order += 1
            power = compose(row, power)
            if power == identity:
                break
        square = compose(row, row)
        cube = compose(square, row)
        chi = sum(row[i] == i for i in range(N))
        chi2 = sum(square[i] == i for i in range(N))
        chi3 = sum(cube[i] == i for i in range(N))
        exterior3 = (chi**3 - 3 * chi * chi2 + 2 * chi3) // 6
        bins[(order, chi, chi2, chi3)] = bins.get((order, chi, chi2, chi3), 0) + 1
        numerator += chi * exterior3
    rows = []
    for (order, chi, chi2, chi3), count in sorted(bins.items()):
        exterior3 = (chi**3 - 3 * chi * chi2 + 2 * chi3) // 6
        rows.append(
            {
                "element_order": order,
                "count": count,
                "chi_V": chi,
                "chi_V_of_square": chi2,
                "chi_V_of_cube": chi3,
                "chi_exterior_cube_V": exterior3,
                "contribution_each": chi * exterior3,
                "contribution_total": count * chi * exterior3,
            }
        )
    check(numerator == 2640 and numerator // len(group) == 44, "TARGET_CHARACTER", "target character differs")
    return {"class_rows": rows, "projector_trace_numerator": numerator, "projector_trace_denominator": len(group), "dimension": numerator // len(group)}


def q_rref(rows: Sequence[Sequence[Fraction]]) -> tuple[list[tuple[Fraction, ...]], list[int]]:
    work = [list(row) for row in rows if any(row)]
    current = 0
    pivots = []
    if not work:
        return [], []
    for column in range(len(work[0])):
        pivot = next((r for r in range(current, len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        work[current], work[pivot] = work[pivot], work[current]
        scale = work[current][column]
        work[current] = [value / scale for value in work[current]]
        for r in range(len(work)):
            if r != current and work[r][column]:
                factor = work[r][column]
                work[r] = [x - factor * y for x, y in zip(work[r], work[current])]
        pivots.append(column)
        current += 1
        if current == len(work):
            break
    return [tuple(row) for row in work[:current]], pivots


def q_rank(rows: Sequence[Sequence[Fraction]]) -> int:
    return len(q_rref(rows)[0])


def lcm(left: int, right: int) -> int:
    return abs(left * right) // math.gcd(left, right) if left and right else 0


def primitive_integer_row(row: Sequence[Fraction]) -> tuple[int, ...]:
    denominator = 1
    for value in row:
        denominator = lcm(denominator, value.denominator)
    integers = [value.numerator * (denominator // value.denominator) for value in row]
    divisor = 0
    for value in integers:
        divisor = math.gcd(divisor, abs(value))
    check(divisor > 0, "PRIMITIVE_ROW", "cannot normalize zero row")
    integers = [value // divisor for value in integers]
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def encode_rational_equation(row: Sequence[Fraction]) -> list[list[int]]:
    return [[a, b, value.numerator, value.denominator] for (a, b), value in zip(MONOMIALS, row) if value]


def encode_integer_equation(row: Sequence[int]) -> list[list[int]]:
    return [[a, b, value] for (a, b), value in zip(MONOMIALS, row) if value]


def encode_k_equation(row: Sequence[K5]) -> list[list[int]]:
    return [[a, b, *encode_k5(value)] for (a, b), value in zip(MONOMIALS, row) if value]


def choose_complement(fixed_rows: Sequence[Sequence[Fraction]], full_rows: Sequence[Sequence[Fraction]]) -> list[tuple[Fraction, ...]]:
    current, _ = q_rref(fixed_rows)
    rank = len(current)
    result = []
    for row in full_rows:
        if q_rank([*current, row]) > rank:
            current.append(tuple(row))
            result.append(tuple(row))
            rank += 1
    check(rank == 38 and len(result) == 27, "RESIDUAL_COMPLEMENT", "canonical residual complement differs")
    return result


def decode_rational_equation(entries: Sequence[Sequence[int]]) -> tuple[Fraction, ...]:
    row = [Fraction(0)] * len(MONOMIALS)
    for a, b, numerator, denominator in entries:
        row[MONOMIAL_INDEX[(a, b)]] = Fraction(numerator, denominator)
    return tuple(row)


def decode_integer_equation(entries: Sequence[Sequence[int]]) -> tuple[Fraction, ...]:
    row = [Fraction(0)] * len(MONOMIALS)
    for a, b, value in entries:
        row[MONOMIAL_INDEX[(a, b)]] = Fraction(value)
    return tuple(row)


def raw_hash(raw: Mapping[Coordinate, Sequence[Fraction]]) -> str:
    encoded = []
    for coordinate in sorted(raw):
        equation = [
            [a, b, value.numerator, value.denominator]
            for (a, b), value in zip(MONOMIALS, raw[coordinate])
            if value
        ]
        encoded.append([list(coordinate), equation])
    return digest(encoded)


def arbitrary_component(raw: Mapping[Coordinate, tuple[Fraction, ...]], output: int, triple: Sequence[int]) -> tuple[Fraction, ...]:
    if len(set(triple)) < 3:
        return tuple(Fraction(0) for _ in MONOMIALS)
    row = raw[(output, *sorted(triple))]
    sign = sign_of_sort(triple)
    return tuple(sign * value for value in row)


def fixed_line_components(raw: Mapping[Coordinate, tuple[Fraction, ...]]) -> list[tuple[Fraction, ...]]:
    rows = []
    for output in range(N):
        for i, j in itertools.combinations(range(N), 2):
            total = tuple(Fraction(0) for _ in MONOMIALS)
            for source in range(N):
                total = add_qrows(total, arbitrary_component(raw, output, (source, i, j)))
            rows.append(total)
    check(len(rows) == 792 and q_rank(rows) == 11, "FIXED_LINE", "fixed-line rank differs")
    return rows


@dataclass(frozen=True)
class K5:
    rational: Fraction = Fraction(0)
    radical: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "rational", Fraction(self.rational))
        object.__setattr__(self, "radical", Fraction(self.radical))

    @staticmethod
    def lift(value: Any) -> "K5":
        return value if isinstance(value, K5) else K5(Fraction(value))

    def __add__(self, other: Any) -> "K5":
        value = K5.lift(other)
        return K5(self.rational + value.rational, self.radical + value.radical)

    __radd__ = __add__

    def __neg__(self) -> "K5":
        return K5(-self.rational, -self.radical)

    def __sub__(self, other: Any) -> "K5":
        return self + (-K5.lift(other))

    def __rsub__(self, other: Any) -> "K5":
        return K5.lift(other) - self

    def __mul__(self, other: Any) -> "K5":
        value = K5.lift(other)
        return K5(
            self.rational * value.rational + 5 * self.radical * value.radical,
            self.rational * value.radical + self.radical * value.rational,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "K5":
        norm = self.rational**2 - 5 * self.radical**2
        check(norm != 0, "K5_DIVISION", "division by zero")
        return K5(self.rational / norm, -self.radical / norm)

    def __truediv__(self, other: Any) -> "K5":
        return self * K5.lift(other).reciprocal()

    def __bool__(self) -> bool:
        return bool(self.rational or self.radical)


ZERO = K5()
ONE = K5(1)
ROOT = K5(0, 1)
KMatrix = list[list[K5]]


def decode_k5(values: Sequence[int]) -> K5:
    check(isinstance(values, list) and len(values) == 4 and all(type(value) is int for value in values), "K5_FORMAT", "bad Q(sqrt(5)) value")
    check(values[1] > 0 and values[3] > 0, "K5_FORMAT", "non-positive Q(sqrt(5)) denominator")
    return K5(Fraction(values[0], values[1]), Fraction(values[2], values[3]))


def encode_k5(value: K5) -> list[int]:
    return [value.rational.numerator, value.rational.denominator, value.radical.numerator, value.radical.denominator]


def decode_k_linear(entries: Sequence[Sequence[int]]) -> tuple[K5, ...]:
    row = [ZERO] * P
    for entry in entries:
        row[entry[0]] = decode_k5(entry[1:])
    return tuple(row)


def decode_k_equation(entries: Sequence[Sequence[int]]) -> tuple[K5, ...]:
    row = [ZERO] * len(MONOMIALS)
    for entry in entries:
        row[MONOMIAL_INDEX[(entry[0], entry[1])]] = decode_k5(entry[2:])
    return tuple(row)


def k_rank(rows: Sequence[Sequence[K5]]) -> int:
    work = [list(row) for row in rows if any(row)]
    if not work:
        return 0
    current = 0
    for column in range(len(work[0])):
        pivot = next((r for r in range(current, len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        work[current], work[pivot] = work[pivot], work[current]
        scale = work[current][column]
        work[current] = [value / scale for value in work[current]]
        for r in range(current + 1, len(work)):
            if work[r][column]:
                factor = work[r][column]
                work[r] = [x - factor * y for x, y in zip(work[r], work[current])]
        current += 1
    return current


def zero_matrix(rows: int, columns: int) -> KMatrix:
    return [[ZERO for _ in range(columns)] for _ in range(rows)]


def identity(size: int) -> KMatrix:
    result = zero_matrix(size, size)
    for i in range(size):
        result[i][i] = ONE
    return result


def matrix_add(left: KMatrix, right: KMatrix) -> KMatrix:
    return [[left[i][j] + right[i][j] for j in range(len(left[0]))] for i in range(len(left))]


def matrix_sub(left: KMatrix, right: KMatrix) -> KMatrix:
    return [[left[i][j] - right[i][j] for j in range(len(left[0]))] for i in range(len(left))]


def matrix_scale(value: K5, matrix: KMatrix) -> KMatrix:
    return [[value * entry for entry in row] for row in matrix]


def matrix_mul(left: KMatrix, right: KMatrix) -> KMatrix:
    result = zero_matrix(len(left), len(right[0]))
    for i in range(len(left)):
        for middle in range(len(right)):
            if not left[i][middle]:
                continue
            for j in range(len(right[0])):
                if right[middle][j]:
                    result[i][j] = result[i][j] + left[i][middle] * right[middle][j]
    return result


def matrix_determinant(matrix: Sequence[Sequence[K5]]) -> K5:
    work = [list(row) for row in matrix]
    result = ONE
    for column in range(len(work)):
        pivot = next((row for row in range(column, len(work)) if work[row][column]), None)
        if pivot is None:
            return ZERO
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        value = work[column][column]
        result = result * value
        for row in range(column + 1, len(work)):
            if work[row][column]:
                factor = work[row][column] / value
                work[row] = [left - factor * right for left, right in zip(work[row], work[column])]
    return result


def matrix_inverse(matrix: Sequence[Sequence[K5]]) -> KMatrix:
    size = len(matrix)
    work = [list(row) + identity(size)[number] for number, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        check(pivot is not None, "CHANNEL_INVERSE", "channel transform is singular")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(size):
            if row != column and work[row][column]:
                factor = work[row][column]
                work[row] = [left - factor * right for left, right in zip(work[row], work[column])]
    return [row[size:] for row in work]


def pair_orbits(group: Sequence[Permutation]) -> list[set[tuple[int, int]]]:
    remaining = {(i, j) for i in range(N) for j in range(N)}
    result = []
    while remaining:
        seed = min(remaining)
        orbit = {(row[seed[0]], row[seed[1]]) for row in group}
        remaining -= orbit
        result.append(orbit)
    return result


def rebuild_projectors(group: Sequence[Permutation]) -> tuple[dict[str, KMatrix], dict[str, Any]]:
    orbitals = sorted(
        (orbit for orbit in pair_orbits(group) if len(orbit) == 60),
        key=lambda orbit: tuple(sorted(orbit)),
    )
    check(len(orbitals) == 2, "PROJECTOR_ORBIT", "wrong valency-five orbital count")
    operator = zero_matrix(N, N)
    for i, j in orbitals[0]:
        operator[i][j] = ONE
    eigenvalues = {"fixed": K5(5), "three_plus": ROOT, "three_minus": -ROOT, "five": K5(-1)}
    projectors = {}
    unit = identity(N)
    square = matrix_mul(operator, operator)
    polynomial = matrix_mul(
        matrix_mul(matrix_sub(operator, matrix_scale(K5(5), unit)), matrix_add(operator, unit)),
        matrix_sub(square, matrix_scale(K5(5), unit)),
    )
    check(polynomial == zero_matrix(N, N), "PROJECTOR_POLYNOMIAL", "orbital minimal polynomial relation fails")
    for label, eigenvalue in eigenvalues.items():
        projector = unit
        for other_label, other in eigenvalues.items():
            if label != other_label:
                projector = matrix_scale((eigenvalue - other).reciprocal(), matrix_mul(projector, matrix_sub(operator, matrix_scale(other, unit))))
        projectors[label] = projector
    dimensions = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}
    total = zero_matrix(N, N)
    for label, projector in projectors.items():
        check(matrix_mul(projector, projector) == projector, "PROJECTOR_IDEMPOTENT", "projector is not idempotent")
        trace = sum((projector[i][i] for i in range(N)), ZERO)
        check(trace == K5(dimensions[label]), "PROJECTOR_TRACE", "projector trace differs")
        check(matrix_mul(operator, projector) == matrix_scale(eigenvalues[label], projector), "PROJECTOR_EIGENVALUE", "projector eigenvalue differs")
        total = matrix_add(total, projector)
    check(total == unit, "PROJECTOR_COMPLETE", "projectors are incomplete")
    labels = list(projectors)
    for number, left in enumerate(labels):
        for right in labels[number + 1 :]:
            check(matrix_mul(projectors[left], projectors[right]) == zero_matrix(N, N), "PROJECTOR_ORTHOGONAL", "projectors overlap")
    metadata = {
        "canonical_valency_five_orbital": [list(pair) for pair in sorted(orbitals[0])],
        "operator_sha256": digest([[encode_k5(value) for value in row] for row in operator]),
        "minimal_polynomial": "(t-5)(t+1)(t^2-5)",
        "sectors": [
            {"id": label, "dimension": dimensions[label], "operator_eigenvalue": encode_k5(eigenvalues[label])}
            for label in labels
        ],
    }
    return projectors, metadata


def constants_from_vectors(vectors: Sequence[Mapping[tuple[int, int, int], Fraction]]) -> list[list[list[list[Fraction]]]]:
    constants = [[[[Fraction(0) for _ in range(N)] for _ in range(N)] for _ in range(N)] for _ in range(P)]
    for parameter, vector in enumerate(vectors):
        for (output, left, right), value in vector.items():
            constants[parameter][output][left][right] = value
            constants[parameter][output][right][left] = -value
    return constants


CHANNEL_SPECS = [
    ("t_pp_to_p", "three_plus", "three_plus", "three_plus", (1, 0, 0)),
    ("t_mm_to_m", "three_minus", "three_minus", "three_minus", (0, 1, 0)),
    ("t_ff_to_p", "five", "five", "three_plus", (-1, 0, 2)),
    ("t_ff_to_m", "five", "five", "three_minus", (0, -1, 2)),
    ("t_pm_to_f", "three_plus", "three_minus", "five", (1, 1, -1)),
    ("t_pf_to_p", "three_plus", "five", "three_plus", (0, 0, 1)),
    ("t_pf_to_m", "three_plus", "five", "three_minus", (1, -1, 1)),
    ("t_pf_to_f", "three_plus", "five", "five", (1, 0, 0)),
    ("t_mf_to_p", "three_minus", "five", "three_plus", (-1, 1, 1)),
    ("t_mf_to_m", "three_minus", "five", "three_minus", (0, 0, 1)),
    ("t_mf_to_f", "three_minus", "five", "five", (0, 1, 0)),
]


def encode_k_linear(row: Sequence[K5]) -> list[list[int]]:
    return [[index, *encode_k5(value)] for index, value in enumerate(row) if value]


def channel_form(
    vectors: Sequence[Mapping[tuple[int, int, int], Fraction]],
    projectors: Mapping[str, KMatrix],
    left: str,
    right: str,
    output_sector: str,
    coordinate: tuple[int, int, int],
) -> tuple[K5, ...]:
    output, left_input, right_input = coordinate
    result = []
    for vector in vectors:
        value = ZERO
        for (source_output, p, q), coefficient in vector.items():
            output_factor = projectors[output_sector][output][source_output]
            input_factor = (
                projectors[left][p][left_input] * projectors[right][q][right_input]
                - projectors[left][q][left_input] * projectors[right][p][right_input]
            )
            value = value + output_factor * input_factor * coefficient
        result.append(value)
    return tuple(result)


def first_channel_form(
    vectors: Sequence[Mapping[tuple[int, int, int], Fraction]],
    projectors: Mapping[str, KMatrix],
    left: str,
    right: str,
    output_sector: str,
) -> tuple[tuple[int, int, int], tuple[K5, ...]]:
    for coordinate in itertools.product(range(N), repeat=3):
        form = channel_form(vectors, projectors, left, right, output_sector, coordinate)
        if any(form):
            return coordinate, form
    raise VerificationError("CHANNEL_ZERO", "canonical channel is zero")


def verify_channel_transform(
    vectors: Sequence[Mapping[tuple[int, int, int], Fraction]],
    group: Sequence[Permutation],
    channel_data: Mapping[str, Any],
) -> tuple[list[tuple[K5, ...]], KMatrix]:
    projectors, spectral_metadata = rebuild_projectors(group)
    constants = constants_from_vectors(vectors)
    derivations = []
    for parameter in range(P):
        matrix = zero_matrix(N, N)
        for output in range(N):
            for column in range(N):
                matrix[output][column] = K5(sum((constants[parameter][output][source][column] for source in range(N)), Fraction(0)))
        derivations.append(matrix)
    transform: list[tuple[K5, ...]] = []
    derivation_rows = []
    meanings = {
        "three_plus": "fixed-line adjoint eigenvalue on three_plus",
        "three_minus": "fixed-line adjoint eigenvalue on three_minus",
        "five": "fixed-line adjoint eigenvalue on five",
    }
    derivation_ids = {"three_plus": "d_plus", "three_minus": "d_minus", "five": "d_five"}
    for sector, dimension in (("three_plus", 3), ("three_minus", 3), ("five", 5)):
        expected: list[K5] = []
        projector = projectors[sector]
        for matrix in derivations:
            trace = ZERO
            for i in range(N):
                for j in range(N):
                    trace = trace + projector[i][j] * matrix[j][i]
            expected.append(trace / dimension)
        form = tuple(expected)
        transform.append(form)
        derivation_rows.append({"id": derivation_ids[sector], "meaning": meanings[sector], "linear_form_in_x": encode_k_linear(form)})
    channel_rows = []
    for channel_id, left, right, output_sector, weight in CHANNEL_SPECS:
        coordinate, form = first_channel_form(vectors, projectors, left, right, output_sector)
        transform.append(form)
        channel_rows.append(
            {
                "id": channel_id,
                "domain_sectors": [left, right],
                "output_sector": output_sector,
                "extraction_coordinate": list(coordinate),
                "linear_form_in_x": encode_k_linear(form),
                "derivation_weight_in_d": list(weight),
            }
        )
    determinant = matrix_determinant(transform)
    check(bool(determinant), "CHANNEL_DETERMINANT", "channel transform is singular")
    inverse_transform = matrix_inverse(transform)
    check(
        matrix_mul([list(row) for row in transform], inverse_transform) == identity(P)
        and matrix_mul(inverse_transform, [list(row) for row in transform]) == identity(P),
        "CHANNEL_INVERSE",
        "independently reconstructed channel inverse differs",
    )
    expected_channel_data = {
        "splitting_field": "Q(sqrt(5)) with sqrt(5)^2=5",
        "spectral_decomposition": spectral_metadata,
        "channel_variables": [*derivation_rows, *channel_rows],
        "transform_determinant": encode_k5(determinant),
        "transform_rows": [encode_k_linear(row) for row in transform],
        "inverse_transform_rows": [encode_k_linear(row) for row in inverse_transform],
    }
    check(channel_data == expected_channel_data, "CHANNEL_DECLARATIONS", "channel, field, determinant, dimension, or transform declaration differs")
    return transform, inverse_transform


def transform_quadratic(row: Sequence[Any], inverse_transform: KMatrix) -> tuple[K5, ...]:
    result = [ZERO] * len(MONOMIALS)
    for coefficient_raw, (i, j) in zip(row, MONOMIALS):
        coefficient = K5.lift(coefficient_raw)
        if not coefficient:
            continue
        for a in range(P):
            if not inverse_transform[i][a]:
                continue
            for b in range(P):
                if inverse_transform[j][b]:
                    index = MONOMIAL_INDEX[tuple(sorted((a, b)))]
                    result[index] = result[index] + coefficient * inverse_transform[i][a] * inverse_transform[j][b]
    return tuple(result)


WEIGHTS = [
    (1, 0, 0),
    (0, 1, 0),
    (-1, 0, 2),
    (0, -1, 2),
    (1, 1, -1),
    (0, 0, 1),
    (1, -1, 1),
    (1, 0, 0),
    (-1, 1, 1),
    (0, 0, 1),
    (0, 1, 0),
]

DISTINCT_WEIGHTS = {
    "d_plus": (1, 0, 0),
    "d_minus": (0, 1, 0),
    "d_five": (0, 0, 1),
    "2d_five_minus_d_plus": (-1, 0, 2),
    "2d_five_minus_d_minus": (0, -1, 2),
    "d_plus_plus_d_minus_minus_d_five": (1, 1, -1),
    "d_plus_plus_d_five_minus_d_minus": (1, -1, 1),
    "d_minus_plus_d_five_minus_d_plus": (-1, 1, 1),
}


def expected_products() -> list[tuple[K5, ...]]:
    rows = []
    for channel, weight in enumerate(WEIGHTS, start=3):
        row = [ZERO] * len(MONOMIALS)
        for d_variable, coefficient in enumerate(weight):
            if coefficient:
                index = MONOMIAL_INDEX[tuple(sorted((d_variable, channel)))]
                row[index] = row[index] + K5(coefficient)
        rows.append(tuple(row))
    return rows


def rref3(rows: Sequence[Sequence[Any]]) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    reduced, _ = q_rref([[Fraction(value) for value in row] for row in rows])
    return tuple(tuple(row) for row in reduced)  # type: ignore[return-value]


def contains3(space: Sequence[Sequence[Fraction]], row: Sequence[int]) -> bool:
    return len(rref3([*space, row])) == len(space)


def recompute_flats(channel_ids: Sequence[str]) -> list[dict[str, Any]]:
    names = list(DISTINCT_WEIGHTS)
    spaces = {
        rref3([DISTINCT_WEIGHTS[name] for name in subset])
        for size in range(len(names) + 1)
        for subset in itertools.combinations(names, size)
    }
    rows = []
    for number, space in enumerate(sorted(spaces, key=lambda value: (len(value), value))):
        vanish = [name for name, weight in DISTINCT_WEIGHTS.items() if contains3(space, weight)]
        allowed = [channel_id for channel_id, weight in zip(channel_ids, WEIGHTS) if contains3(space, weight)]
        rows.append(
            {
                "flat_id": f"F{number:02d}",
                "dimension": 3 - len(space),
                "defining_rowspace_rref": [
                    [[value.numerator, value.denominator] for value in row]
                    for row in space
                ],
                "weights_vanishing_identically": vanish,
                "channels_not_forced_zero_by_fixed_line_sector": allowed,
                "channels_forced_zero_on_relative_open_flat": [channel_id for channel_id in channel_ids if channel_id not in allowed],
                "scope": "relative-open part of this flat; larger subflats are separate rows",
            }
        )
    check(len(rows) == 28, "FLAT_COUNT", "weight arrangement count differs")
    return rows


def validate_integer_entries(rows: Any, width: int, code: str) -> None:
    check(isinstance(rows, list), code, "equation collection is not a list")
    for row in rows:
        check(isinstance(row, list), code, "equation row is not a list")
        for entry in row:
            check(isinstance(entry, list) and len(entry) == width and all(type(value) is int for value in entry), code, "bad encoded equation entry")


def validate_stage2_keysets(system: Any, receipt: Any) -> None:
    check_keys(system, SYSTEM_KEYS, "SYSTEM_KEYS")
    check_keys(system.get("upstream"), UPSTREAM_KEYS, "SYSTEM_UPSTREAM_KEYS")
    jacobi = system.get("jacobi_tensor")
    check_keys(jacobi, JACOBI_KEYS, "SYSTEM_JACOBI_KEYS")
    character = jacobi.get("target_invariant_character_certificate")
    check_keys(character, TARGET_CHARACTER_KEYS, "SYSTEM_TARGET_CHARACTER_KEYS")
    check(isinstance(character.get("class_rows"), list), "SYSTEM_TARGET_CLASS_ROWS", "target character rows are absent")
    for row in character["class_rows"]:
        check_keys(row, TARGET_CLASS_ROW_KEYS, "SYSTEM_TARGET_CLASS_ROW_KEYS")
    representatives = jacobi.get("representative_equations")
    check(isinstance(representatives, list), "SYSTEM_REPRESENTATIVE_ROWS", "representatives are absent")
    for row in representatives:
        check_keys(row, REPRESENTATIVE_KEYS, "SYSTEM_REPRESENTATIVE_KEYS")
        validate_integer_entries([row.get("equation")], 4, "SYSTEM_RATIONAL_EQUATION_FORMAT")
    validate_integer_entries(jacobi.get("reduced_integer_equations"), 3, "SYSTEM_INTEGER_EQUATION_FORMAT")
    fixed = system.get("fixed_line_reduction")
    check_keys(fixed, FIXED_REDUCTION_KEYS, "SYSTEM_FIXED_KEYS")
    validate_integer_entries(fixed.get("contracted_integer_equations_in_x"), 3, "SYSTEM_INTEGER_EQUATION_FORMAT")
    validate_integer_entries(fixed.get("residual_integer_equations_in_x"), 3, "SYSTEM_INTEGER_EQUATION_FORMAT")
    validate_integer_entries(fixed.get("residual_equations_in_channel_variables"), 6, "SYSTEM_Q5_EQUATION_FORMAT")
    channel = fixed.get("channel_decomposition")
    check_keys(channel, CHANNEL_DECOMPOSITION_KEYS, "SYSTEM_CHANNEL_KEYS")
    variables = channel.get("channel_variables")
    check(isinstance(variables, list) and len(variables) == P, "SYSTEM_CHANNEL_VARIABLES", "channel-variable list differs")
    for number, row in enumerate(variables):
        check_keys(row, D_CHANNEL_KEYS if number < 3 else T_CHANNEL_KEYS, "SYSTEM_CHANNEL_VARIABLE_KEYS")
        validate_integer_entries([row.get("linear_form_in_x")], 5, "SYSTEM_Q5_LINEAR_FORMAT")
    validate_integer_entries(channel.get("transform_rows"), 5, "SYSTEM_Q5_LINEAR_FORMAT")
    validate_integer_entries(channel.get("inverse_transform_rows"), 5, "SYSTEM_Q5_LINEAR_FORMAT")
    spectral = channel.get("spectral_decomposition")
    check_keys(spectral, SPECTRAL_KEYS, "SYSTEM_SPECTRAL_KEYS")
    sectors = spectral.get("sectors")
    check(isinstance(sectors, list), "SYSTEM_SECTOR_ROWS", "spectral sectors are absent")
    for row in sectors:
        check_keys(row, SECTOR_KEYS, "SYSTEM_SECTOR_KEYS")
    products = fixed.get("derivation_product_equations_in_channel_variables")
    check(isinstance(products, list), "SYSTEM_PRODUCT_ROWS", "derivation products are absent")
    for row in products:
        check_keys(row, PRODUCT_DECLARATION_KEYS, "SYSTEM_PRODUCT_KEYS")
        validate_integer_entries([row.get("equation")], 6, "SYSTEM_Q5_EQUATION_FORMAT")
    check_keys(fixed.get("derivation_only_subspace"), DERIVATION_ONLY_KEYS, "SYSTEM_DERIVATION_KEYS")
    check_keys(fixed.get("central_specialization"), CENTRAL_SPECIALIZATION_KEYS, "SYSTEM_CENTRAL_KEYS")
    arrangement = fixed.get("derivation_weight_arrangement")
    check_keys(arrangement, ARRANGEMENT_KEYS, "SYSTEM_ARRANGEMENT_KEYS")
    check_keys(arrangement.get("flat_dimension_counts"), {"0", "1", "2", "3"}, "SYSTEM_FLAT_COUNT_KEYS")
    flats = arrangement.get("relative_open_flats")
    check(isinstance(flats, list), "SYSTEM_FLAT_ROWS", "arrangement flats are absent")
    for row in flats:
        check_keys(row, FLAT_KEYS, "SYSTEM_FLAT_KEYS")
    check_keys(receipt, RECEIPT_KEYS, "RECEIPT_KEYS")
    check_keys(receipt.get("semantic_inputs"), SEMANTIC_INPUT_KEYS, "RECEIPT_SEMANTIC_KEYS")
    check_keys(receipt["semantic_inputs"].get("paths"), UPSTREAM_KEYS, "RECEIPT_SEMANTIC_PATH_KEYS")
    firewall = receipt.get("target_firewall")
    check_keys(firewall, FIREWALL_KEYS, "RECEIPT_FIREWALL_KEYS")
    for field in ("audited_import_roots_by_file", "exact_import_allowlists_by_file"):
        check_keys(firewall.get(field), set(EXACT_IMPORT_ROOTS_BY_FILE), "RECEIPT_IMPORT_FILE_KEYS")
    check_keys(firewall.get("runtime_forbidden_import_roots_by_file"), set(RUNTIME_FORBIDDEN_IMPORT_ROOTS), "RECEIPT_FORBIDDEN_FILE_KEYS")
    check_keys(receipt.get("attained"), ATTAINED_KEYS, "RECEIPT_ATTAINED_KEYS")
    check_keys(receipt.get("not_attained"), NOT_ATTAINED_KEYS, "RECEIPT_NOT_ATTAINED_KEYS")
    check_keys(receipt.get("system_artifact"), SYSTEM_ARTIFACT_KEYS, "RECEIPT_SYSTEM_ARTIFACT_KEYS")
    check_keys(receipt.get("implementation_pins"), IMPLEMENTATION_PIN_KEYS, "RECEIPT_PIN_KEYS")
    mutations = receipt.get("mutation_tests")
    check(isinstance(mutations, list), "RECEIPT_MUTATION_ROWS", "mutation rows are absent")
    for row in mutations:
        check_keys(row, MUTATION_ROW_KEYS, "RECEIPT_MUTATION_KEYS")


def expected_producer_mutations() -> list[dict[str, Any]]:
    rows = [
        ("inject_unregistered_upstream_target", "FIREWALL_BASIS_KEYS"),
        ("tamper_raw_jacobi_coordinate", "JACOBI_COVARIANCE"),
        ("duplicate_channel_coordinate", "CHANNEL_TRANSFORM_SINGULAR"),
        ("alter_derivation_product", "PRODUCT_ROWSPACE"),
        ("drop_residual_equation", "RESIDUAL_COUNT"),
        ("duplicate_residual_equation", "RESIDUAL_SPAN"),
        ("promote_verifier_runtime_import_boundary", "FIREWALL_RUNTIME_IMPORT"),
        ("rehashed_false_full_classification_claim", "RECEIPT_BOUNDARY"),
        ("rehashed_false_preferred_algebra_claim", "RECEIPT_ATTAINED_KEYS"),
        ("rehashed_false_source_selection_claim", "RECEIPT_BOUNDARY"),
        ("rehashed_false_compactness_claim", "RECEIPT_BOUNDARY"),
        ("rehashed_receipt_issue_change", "RECEIPT_ISSUE"),
        ("rehashed_receipt_count_change", "RECEIPT_COUNTS"),
        ("rehashed_semantic_path_change", "RECEIPT_SEMANTIC_INPUTS"),
        ("rehashed_firewall_change", "RECEIPT_FIREWALL"),
        ("rehashed_receipt_extra_key", "RECEIPT_KEYS"),
        ("rehashed_system_issue_change", "SYSTEM_ISSUE"),
        ("rehashed_system_count_change", "SYSTEM_JACOBI_CLAIMS"),
        ("rehashed_system_rank_change", "SYSTEM_RANKS"),
        ("rehashed_system_split_change", "SYSTEM_SPLIT"),
        ("rehashed_system_field_change", "SYSTEM_FIELD"),
        ("rehashed_minimal_polynomial_change", "SYSTEM_MINIMAL_POLYNOMIAL"),
        ("rehashed_transform_determinant_change", "SYSTEM_DETERMINANT"),
        ("rehashed_sector_dimension_change", "SYSTEM_SECTOR_DIMENSIONS"),
        ("rehashed_system_extra_key", "SYSTEM_KEYS"),
    ]
    return [{"name": name, "expected_error": code, "actual_error": code, "passed": True} for name, code in rows]


def validate_system_claims(
    system: Mapping[str, Any],
    basis: Mapping[str, Any],
    stage1_receipt: Mapping[str, Any],
    group: Sequence[Permutation],
    raw: Mapping[Coordinate, Sequence[Fraction]],
    orbits: Sequence[Mapping[Coordinate, int]],
    representatives: Sequence[Sequence[Fraction]],
    full_rref: Sequence[Sequence[Fraction]],
    full_pivots: Sequence[int],
    fixed_components: Sequence[Sequence[Fraction]],
    fixed_rref: Sequence[Sequence[Fraction]],
    fixed_pivots: Sequence[int],
    residual_rref: Sequence[Sequence[Fraction]],
    products: Sequence[Sequence[K5]],
    residual_y: Sequence[Sequence[K5]],
    central_rank: int,
    expected_flats: Sequence[Mapping[str, Any]],
) -> None:
    check(system.get("schema") == "oph.a5_jacobi_system_reduction.v1", "SYSTEM_SCHEMA", "wrong system schema")
    check(system.get("issue") == 566, "SYSTEM_ISSUE", "wrong system issue")
    check(system.get("field") == "Q", "SYSTEM_FIELD", "wrong system field")
    check(system.get("parameter_variables") == [f"x{i:02d}" for i in range(P)], "SYSTEM_VARIABLES", "parameter variables differ")
    check(system.get("monomial_order") == [list(row) for row in MONOMIALS], "SYSTEM_MONOMIALS", "monomial order differs")
    expected_upstream = {
        "basis_path": "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json",
        "basis_file_sha256": file_digest(STAGE1_BASIS),
        "basis_canonical_json_sha256": digest(basis),
        "receipt_path": "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json",
        "receipt_file_sha256": file_digest(STAGE1_RECEIPT),
        "receipt_self_sha256": stage1_receipt["receipt_sha256"],
        "proper_group_sha256": stage1_receipt["proper_port_action"]["permutation_rows_sha256"],
    }
    check(system.get("upstream") == expected_upstream, "SYSTEM_UPSTREAM", "upstream path or pin differs")
    jacobi = system["jacobi_tensor"]
    expected_representatives = []
    for orbit in orbits:
        seed = min(orbit)
        expected_representatives.append(
            {"target_coordinate": list(seed), "signed_orbit_size": len(orbit), "equation": encode_rational_equation(raw[seed])}
        )
    expected_reduced = [encode_integer_equation(primitive_integer_row(row)) for row in full_rref]
    check(
        jacobi.get("convention") == "J(u,v,w)=B(B(u,v),w)+B(B(v,w),u)+B(B(w,u),v)"
        and jacobi.get("coordinate_equation_count") == len(raw) == 2640
        and jacobi.get("all_coordinate_equations_nonzero") is all(any(row) for row in raw.values()) is True
        and jacobi.get("all_coordinate_equations_sha256") == raw_hash(raw)
        and jacobi.get("target_space_dimension") == N * math.comb(N, 3) == 2640
        and jacobi.get("target_invariant_character_certificate") == target_character_claim(group)
        and jacobi.get("signed_target_orbit_count") == len(orbits) == 44
        and jacobi.get("signed_target_orbit_sizes") == [len(orbit) for orbit in orbits] == [60] * 44
        and jacobi.get("representative_equations") == expected_representatives
        and jacobi.get("representative_coefficient_rank_over_Q") == len(full_rref) == q_rank(representatives) == 38
        and jacobi.get("representative_rref_pivot_monomial_indices") == list(full_pivots)
        and jacobi.get("reduced_integer_equations") == expected_reduced
        and jacobi.get("reduced_integer_equation_count") == len(expected_reduced) == 38
        and jacobi.get("zero_locus_equivalence")
        == "all 2640 components vanish iff the 44 orbit representatives vanish iff the coefficient-row rank 38 reduced integer quadrics vanish",
        "SYSTEM_JACOBI_DECLARATIONS",
        "Jacobi count, coefficient-row rank, character, or scope declaration differs",
    )
    fixed = system["fixed_line_reduction"]
    expected_fixed = [encode_integer_equation(primitive_integer_row(row)) for row in fixed_rref]
    expected_residual = [encode_integer_equation(primitive_integer_row(row)) for row in residual_rref]
    check(
        fixed.get("fixed_vector") == [1] * N
        and fixed.get("contracted_coordinate_count") == len(fixed_components) == 792
        and fixed.get("contracted_coefficient_rank_over_Q") == len(fixed_rref) == 11
        and fixed.get("contracted_rref_pivot_monomial_indices") == list(fixed_pivots)
        and fixed.get("contracted_integer_equations_in_x") == expected_fixed,
        "SYSTEM_FIXED_DECLARATIONS",
        "fixed-line count, coefficient-row rank, pivots, or equations differ",
    )
    expected_product_rows = [
        {"channel_id": spec[0], "weight_coefficients_on_d_plus_d_minus_d_five": list(spec[4]), "equation": encode_k_equation(row)}
        for spec, row in zip(CHANNEL_SPECS, products)
    ]
    check(
        fixed.get("derivation_product_equations_in_channel_variables") == expected_product_rows
        and fixed.get("product_equation_rank_over_Qsqrt5") == k_rank(products) == 11
        and fixed.get("product_rowspace_equals_contracted_jacobi_rowspace") is True
        and fixed.get("residual_integer_equations_in_x") == expected_residual
        and fixed.get("residual_equations_in_channel_variables") == [encode_k_equation(row) for row in residual_y]
        and fixed.get("residual_rank_over_Qsqrt5") == k_rank(residual_y) == 27,
        "SYSTEM_REDUCTION_DECLARATIONS",
        "product or residual declaration differs",
    )
    check(fixed.get("full_rank_split") == [11, 27, 38], "SYSTEM_SPLIT", "11+27 rowspace split differs")
    check(
        fixed.get("derivation_only_subspace")
        == {"dimension": 3, "definition": "all eleven non-derivation channel variables are zero",
            "all_jacobi_equations_vanish_identically": True, "scope": DERIVATION_ONLY_SCOPE},
        "SYSTEM_DERIVATION_SCOPE",
        "derivation-only scope differs",
    )
    check(
        fixed.get("central_specialization")
        == {"definition": "d_plus=d_minus=d_five=0", "remaining_channel_variable_count": 11,
            "residual_coefficient_rank_over_Qsqrt5": central_rank, "solution_set_classified": False},
        "SYSTEM_CENTRAL_SPECIALIZATION",
        "central specialization declaration differs",
    )
    arrangement = fixed["derivation_weight_arrangement"]
    expected_counts = {str(dimension): sum(row["dimension"] == dimension for row in expected_flats) for dimension in range(4)}
    check(
        arrangement.get("distinct_weight_hyperplane_count") == len(DISTINCT_WEIGHTS) == 8
        and arrangement.get("flat_count") == len(expected_flats) == 28
        and arrangement.get("flat_dimension_counts") == expected_counts == {"0": 1, "1": 18, "2": 8, "3": 1}
        and arrangement.get("relative_open_flats") == list(expected_flats)
        and arrangement.get("scope") == ARRANGEMENT_SCOPE,
        "SYSTEM_ARRANGEMENT_DECLARATIONS",
        "arrangement count, flats, or residual scope differs",
    )


def validate_receipt_claims(
    system: Mapping[str, Any], receipt: Mapping[str, Any], imports: Mapping[str, Sequence[str]]
) -> None:
    unhashed = dict(receipt)
    stored = unhashed.pop("receipt_sha256", None)
    check(stored == digest(unhashed), "RECEIPT_HASH", "receipt self-hash differs")
    check(receipt.get("schema") == "oph.a5_jacobi_stage2.receipt.v1", "RECEIPT_SCHEMA", "wrong receipt schema")
    check(receipt.get("issue") == 566, "RECEIPT_ISSUE", "wrong receipt issue")
    check(receipt.get("status") == STATUS, "RECEIPT_STATUS", "status differs")
    check(receipt.get("claim_boundary") == CLAIM_BOUNDARY, "RECEIPT_SCOPE", "claim boundary differs")
    check(
        receipt.get("semantic_inputs") == {"count": 2, "paths": system["upstream"], "stage1_later_gates_all_false": True},
        "RECEIPT_SEMANTIC_INPUTS", "semantic-input declaration differs",
    )
    expected_firewall = {
        "enabled": True,
        "accepted_input_schemas": ["oph.a5_alternating_bracket_reynolds_basis.v1", "oph.a5_alternating_bracket_space_stage1.receipt.v1"],
        "exact_input_keysets_enforced": True,
        "arbitrary_input_path_allowed": False,
        "coefficient_target_input_allowed": False,
        "measurement_input_allowed": False,
        "environment_target_read": False,
        "network_target_read": False,
        "nonstdlib_imports": [],
        "audited_import_roots_by_file": dict(imports),
        "exact_import_allowlists_by_file": {filename: sorted(roots) for filename, roots in sorted(EXACT_IMPORT_ROOTS_BY_FILE.items())},
        "runtime_forbidden_import_roots_by_file": {filename: sorted(roots) for filename, roots in sorted(RUNTIME_FORBIDDEN_IMPORT_ROOTS.items())},
    }
    check(receipt.get("target_firewall") == expected_firewall, "RECEIPT_FIREWALL", "firewall declaration differs")
    expected_attained = {
        "general_parameter_count": 14, "raw_jacobi_coordinate_count": 2640, "equivariant_target_orbit_count": 44,
        "independent_jacobi_quadrics": 38, "fixed_line_product_quadrics": 11, "residual_independent_quadrics": 27,
        "exact_channel_coordinate_change": True, "derivation_weight_relative_open_flat_count": 28,
        "generic_weight_stratum_forces_all_non_derivation_channels_zero": True,
        "derivation_only_three_parameter_family_certified": True,
    }
    check(receipt.get("attained") == expected_attained, "RECEIPT_COUNTS", "attained counts or ranks differ")
    check(receipt.get("not_attained") == {key: False for key in NOT_ATTAINED_KEYS}, "RECEIPT_BOUNDARY", "non-attained boundary differs")
    check(
        receipt.get("system_artifact")
        == {"path": "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_system_reduction.json", "canonical_json_sha256": digest(system)},
        "RECEIPT_SYSTEM_ARTIFACT", "system artifact path or hash differs",
    )
    check(receipt.get("mutation_tests") == expected_producer_mutations(), "RECEIPT_MUTATIONS", "producer mutation declaration differs")
    pins = receipt["implementation_pins"]
    check(
        pins.get("producer_sha256") == file_digest(PRODUCER_PATH)
        and pins.get("independent_verifier_sha256") == file_digest(Path(__file__).resolve())
        and pins.get("test_sha256") == file_digest(TEST_PATH),
        "RECEIPT_IMPLEMENTATION_PIN", "implementation pin differs",
    )


def validate_reduced_rowspace(representatives: Sequence[Sequence[Fraction]], reduced: Sequence[Sequence[Fraction]]) -> None:
    check(len(reduced) == 38 and q_rank(reduced) == 38, "REDUCED_RANK", "reduced Jacobi system rank differs")
    check(q_rank([*representatives, *reduced]) == 38, "REDUCED_ROWSPACE", "reduced Jacobi rowspace differs")


def validate_k_decomposition(
    full: Sequence[Sequence[K5]], fixed: Sequence[Sequence[K5]], products: Sequence[Sequence[K5]], residual: Sequence[Sequence[K5]]
) -> None:
    check(k_rank(fixed) == k_rank(products) == 11, "PRODUCT_RANK", "fixed or product rank differs")
    check(k_rank([*fixed, *products]) == 11, "PRODUCT_ROWSPACE", "product rowspace differs")
    check(len(residual) == 27 and k_rank(residual) == 27, "RESIDUAL_RANK", "residual rank differs")
    check(k_rank(full) == k_rank([*products, *residual]) == k_rank([*full, *products, *residual]) == 38, "RESIDUAL_SPAN", "11+27 split changes full rowspace")


def validate_boundary(receipt: Mapping[str, Any]) -> None:
    check(
        receipt.get("status") == "EXACT_JACOBI_SYSTEM_AND_FIXED_LINE_REDUCTION__FULL_CLASSIFICATION_OPEN",
        "BOUNDARY_STATUS",
        "status overstates the result",
    )
    not_attained = receipt.get("not_attained")
    check(isinstance(not_attained, Mapping) and not_attained and not any(not_attained.values()), "BOUNDARY_PROMOTION", "a non-attained gate was promoted")
    firewall = receipt.get("target_firewall")
    check(
        isinstance(firewall, Mapping)
        and firewall.get("enabled") is True
        and firewall.get("arbitrary_input_path_allowed") is False
        and firewall.get("coefficient_target_input_allowed") is False
        and firewall.get("measurement_input_allowed") is False,
        "BOUNDARY_FIREWALL",
        "target firewall is open",
    )


def validate_import_contract(filename: str, roots: set[str]) -> None:
    forbidden = RUNTIME_FORBIDDEN_IMPORT_ROOTS.get(filename, set())
    check(not (roots & forbidden), "IMPORT_RUNTIME_FORBIDDEN", f"{filename} imports a forbidden runtime root")
    check(roots == EXACT_IMPORT_ROOTS_BY_FILE[filename], "IMPORT_EXACT_CONTRACT", f"{filename} roots differ from exact contract")


def audit_imports() -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for path in (PRODUCER_PATH, Path(__file__).resolve(), TEST_PATH):
        roots: set[str] = set()
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                roots.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                roots.add(node.module.split(".")[0])
        validate_import_contract(path.name, roots)
        result[path.name] = sorted(roots)
    check(set(result) == set(EXACT_IMPORT_ROOTS_BY_FILE), "IMPORT_FILES", "per-file import audit is incomplete")
    return dict(sorted(result.items()))


def validate_import_receipt(receipt: Mapping[str, Any], actual: Mapping[str, Sequence[str]]) -> None:
    firewall = receipt.get("target_firewall")
    check(isinstance(firewall, Mapping), "IMPORT_RECEIPT", "target firewall is absent")
    expected = {filename: sorted(roots) for filename, roots in sorted(EXACT_IMPORT_ROOTS_BY_FILE.items())}
    forbidden = {filename: sorted(roots) for filename, roots in sorted(RUNTIME_FORBIDDEN_IMPORT_ROOTS.items())}
    reported = firewall.get("audited_import_roots_by_file")
    check(isinstance(reported, Mapping), "IMPORT_RECEIPT", "per-file audited roots are absent")
    check(
        "subprocess" not in reported.get("certify.py", []) and "subprocess" not in reported.get("verify.py", []),
        "IMPORT_BOUNDARY_PROMOTION",
        "reported runtime import boundary includes subprocess",
    )
    check(reported == dict(actual), "IMPORT_RECEIPT", "per-file audited roots differ")
    check(firewall.get("exact_import_allowlists_by_file") == expected, "IMPORT_RECEIPT", "per-file exact allowlists differ")
    check(firewall.get("runtime_forbidden_import_roots_by_file") == forbidden, "IMPORT_RECEIPT", "runtime forbidden roots differ")
    check(
        "subprocess" not in actual["certify.py"] and "subprocess" not in actual["verify.py"],
        "IMPORT_BOUNDARY_PROMOTION",
        "runtime file gained subprocess",
    )


def catches(code: str, action: Any) -> bool:
    try:
        action()
    except VerificationError as exc:
        return exc.code == code
    return False


def independent_mutations(
    system: Mapping[str, Any],
    receipt: Mapping[str, Any],
    representatives: Sequence[tuple[Fraction, ...]],
    reduced: Sequence[tuple[Fraction, ...]],
    full_y: Sequence[tuple[K5, ...]],
    fixed_y: Sequence[tuple[K5, ...]],
    products: Sequence[tuple[K5, ...]],
    residual_y: Sequence[tuple[K5, ...]],
    transform: Sequence[tuple[K5, ...]],
    inverse_transform: KMatrix,
    artifact_validator: Any,
) -> dict[str, bool]:
    result: dict[str, bool] = {}

    def rehash(changed_receipt: dict[str, Any]) -> None:
        unhashed = dict(changed_receipt)
        unhashed.pop("receipt_sha256", None)
        changed_receipt["receipt_sha256"] = digest(unhashed)

    def receipt_case(name: str, expected: str, mutate: Any) -> None:
        changed_receipt = copy.deepcopy(receipt)
        mutate(changed_receipt)
        rehash(changed_receipt)
        result[name] = catches(expected, lambda: artifact_validator(system, changed_receipt))

    def system_case(name: str, expected: str, mutate: Any) -> None:
        changed_system = copy.deepcopy(system)
        changed_receipt = copy.deepcopy(receipt)
        mutate(changed_system)
        changed_receipt["system_artifact"]["canonical_json_sha256"] = digest(changed_system)
        changed_receipt["semantic_inputs"]["paths"] = copy.deepcopy(changed_system["upstream"])
        rehash(changed_receipt)
        result[name] = catches(expected, lambda: artifact_validator(changed_system, changed_receipt))

    receipt_case(
        "rehashed_false_full_classification_caught", "RECEIPT_BOUNDARY",
        lambda value: value["not_attained"].__setitem__("full_jacobi_solution_variety_classified", True),
    )
    receipt_case(
        "rehashed_false_preferred_algebra_caught", "RECEIPT_ATTAINED_KEYS",
        lambda value: value["attained"].__setitem__("standard_model_identified", True),
    )
    receipt_case(
        "rehashed_false_source_selection_caught", "RECEIPT_BOUNDARY",
        lambda value: value["not_attained"].__setitem__("source_selection", True),
    )
    receipt_case(
        "rehashed_false_compactness_caught", "RECEIPT_BOUNDARY",
        lambda value: value["not_attained"].__setitem__("compactness_established", True),
    )
    receipt_case("rehashed_receipt_issue_change_caught", "RECEIPT_ISSUE", lambda value: value.__setitem__("issue", 567))
    receipt_case(
        "rehashed_receipt_count_change_caught", "RECEIPT_COUNTS",
        lambda value: value["attained"].__setitem__("raw_jacobi_coordinate_count", 2639),
    )
    receipt_case(
        "rehashed_semantic_path_change_caught", "RECEIPT_SEMANTIC_INPUTS",
        lambda value: value["semantic_inputs"]["paths"].__setitem__("basis_path", "unregistered.json"),
    )
    receipt_case(
        "rehashed_firewall_change_caught", "RECEIPT_FIREWALL",
        lambda value: value["target_firewall"].__setitem__("measurement_input_allowed", True),
    )
    receipt_case("rehashed_receipt_extra_key_caught", "RECEIPT_KEYS", lambda value: value.__setitem__("extra_claim", True))
    receipt_case(
        "rehashed_import_boundary_promotion_caught", "RECEIPT_FIREWALL",
        lambda value: value["target_firewall"]["audited_import_roots_by_file"]["verify.py"].append("subprocess"),
    )

    system_case("rehashed_system_issue_change_caught", "SYSTEM_ISSUE", lambda value: value.__setitem__("issue", 567))
    system_case("rehashed_system_field_change_caught", "SYSTEM_FIELD", lambda value: value.__setitem__("field", "R"))
    system_case(
        "rehashed_system_path_change_caught", "SYSTEM_UPSTREAM",
        lambda value: value["upstream"].__setitem__("basis_path", "unregistered.json"),
    )
    system_case(
        "rehashed_system_count_change_caught", "SYSTEM_JACOBI_DECLARATIONS",
        lambda value: value["jacobi_tensor"].__setitem__("coordinate_equation_count", 2639),
    )
    system_case(
        "rehashed_system_rank_change_caught", "SYSTEM_REDUCTION_DECLARATIONS",
        lambda value: value["fixed_line_reduction"].__setitem__("residual_rank_over_Qsqrt5", 26),
    )
    system_case(
        "rehashed_system_split_change_caught", "SYSTEM_SPLIT",
        lambda value: value["fixed_line_reduction"].__setitem__("full_rank_split", [11, 26, 37]),
    )
    system_case(
        "rehashed_minimal_polynomial_change_caught", "CHANNEL_DECLARATIONS",
        lambda value: value["fixed_line_reduction"]["channel_decomposition"]["spectral_decomposition"].__setitem__("minimal_polynomial", "(t-5)(t+1)"),
    )
    system_case(
        "rehashed_transform_determinant_change_caught", "CHANNEL_DECLARATIONS",
        lambda value: value["fixed_line_reduction"]["channel_decomposition"].__setitem__("transform_determinant", [1, 1, 0, 1]),
    )
    system_case(
        "rehashed_sector_dimension_change_caught", "CHANNEL_DECLARATIONS",
        lambda value: value["fixed_line_reduction"]["channel_decomposition"]["spectral_decomposition"]["sectors"][1].__setitem__("dimension", 4),
    )
    system_case(
        "rehashed_reduced_equation_tamper_caught", "SYSTEM_JACOBI_DECLARATIONS",
        lambda value: value["jacobi_tensor"]["reduced_integer_equations"][0].append([13, 13, 1]),
    )
    system_case(
        "rehashed_product_equation_tamper_caught", "SYSTEM_REDUCTION_DECLARATIONS",
        lambda value: value["fixed_line_reduction"]["derivation_product_equations_in_channel_variables"][0]["equation"].append([3, 3, 1, 1, 0, 1]),
    )
    system_case(
        "rehashed_residual_drop_caught", "SYSTEM_REDUCTION_DECLARATIONS",
        lambda value: value["fixed_line_reduction"]["residual_equations_in_channel_variables"].pop(),
    )
    system_case(
        "rehashed_channel_inverse_tamper_caught", "CHANNEL_DECLARATIONS",
        lambda value: value["fixed_line_reduction"]["channel_decomposition"]["inverse_transform_rows"][0].append([13, 1, 1, 0, 1]),
    )
    system_case("rehashed_system_extra_key_caught", "SYSTEM_KEYS", lambda value: value.__setitem__("extra_claim", True))

    # Separate algebra-only mutations make sure the rowspace checks themselves,
    # not only exact serialization comparisons, remain active.
    bad_products = [tuple(row) for row in products]
    row = list(bad_products[0])
    row[MONOMIAL_INDEX[(3, 3)]] = row[MONOMIAL_INDEX[(3, 3)]] + ONE
    bad_products[0] = tuple(row)
    result["algebraic_product_rowspace_tamper_caught"] = catches(
        "PRODUCT_ROWSPACE", lambda: validate_k_decomposition(full_y, fixed_y, bad_products, residual_y)
    )

    bad_inverse = [list(row) for row in inverse_transform]
    bad_inverse[0][0] = bad_inverse[0][0] + ONE
    result["algebraic_channel_inverse_tamper_caught"] = matrix_mul([list(row) for row in transform], bad_inverse) != identity(P)
    result["algebraic_residual_drop_caught"] = catches(
        "RESIDUAL_RANK", lambda: validate_k_decomposition(full_y, fixed_y, products, residual_y[:-1])
    )
    mutated_reduced = [tuple(row) for row in reduced]
    changed = list(mutated_reduced[0])
    changed[-1] += 1
    mutated_reduced[0] = tuple(changed)
    result["algebraic_reduced_rowspace_tamper_caught"] = catches(
        "REDUCED_ROWSPACE", lambda: validate_reduced_rowspace(representatives, mutated_reduced)
    )
    check(all(result.values()), "INDEPENDENT_MUTATION", "independent semantic mutation survived")
    return result


def verify() -> dict[str, Any]:
    basis_raw, stage1_receipt, vectors, group = decode_stage1()
    system = json.loads(SYSTEM_PATH.read_text(encoding="utf-8"))
    receipt = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))
    check(isinstance(system, dict) and isinstance(receipt, dict), "ARTIFACT_JSON", "stage-two roots must be objects")
    validate_stage2_keysets(system, receipt)
    unhashed = dict(receipt)
    stored_receipt_hash = unhashed.pop("receipt_sha256", None)
    check(stored_receipt_hash == digest(unhashed), "RECEIPT_HASH", "receipt self-hash differs")
    check(receipt["system_artifact"]["canonical_json_sha256"] == digest(system), "SYSTEM_HASH", "system hash differs")
    validate_boundary(receipt)

    upstream = system["upstream"]
    check(upstream["basis_file_sha256"] == file_digest(STAGE1_BASIS), "UPSTREAM_FILE_HASH", "basis file hash differs")
    check(upstream["receipt_file_sha256"] == file_digest(STAGE1_RECEIPT), "UPSTREAM_FILE_HASH", "receipt file hash differs")
    check(upstream["basis_canonical_json_sha256"] == digest(basis_raw), "UPSTREAM_OBJECT_HASH", "basis object hash differs")
    check(upstream["receipt_self_sha256"] == stage1_receipt["receipt_sha256"], "UPSTREAM_OBJECT_HASH", "receipt self pin differs")

    forms = coefficient_forms(vectors)
    raw = jacobi_from_forms(forms)
    jacobi_artifact = system["jacobi_tensor"]
    check(raw_hash(raw) == jacobi_artifact["all_coordinate_equations_sha256"], "RAW_HASH", "all-component hash differs")
    orbits = signed_orbits(group)
    check(len(jacobi_artifact["representative_equations"]) == len(orbits) == 44, "REPRESENTATIVE_COUNT", "representative count differs")
    representatives = []
    for orbit, artifact_row in zip(orbits, jacobi_artifact["representative_equations"]):
        seed = min(orbit)
        check(artifact_row["target_coordinate"] == list(seed), "REPRESENTATIVE_COORDINATE", "orbit representative differs")
        decoded = decode_rational_equation(artifact_row["equation"])
        check(decoded == raw[seed], "REPRESENTATIVE_EQUATION", "representative equation differs")
        for coordinate, sign in orbit.items():
            check(raw[coordinate] == tuple(sign * value for value in raw[seed]), "RAW_COVARIANCE", "raw tensor covariance fails")
        representatives.append(decoded)
    check(len(representatives) == 44 and q_rank(representatives) == 38, "REPRESENTATIVE_RANK", "representative rank differs")
    full_rref, full_pivots = q_rref(representatives)
    reduced = [decode_integer_equation(row) for row in jacobi_artifact["reduced_integer_equations"]]
    validate_reduced_rowspace(representatives, reduced)
    check(
        jacobi_artifact["reduced_integer_equations"]
        == [encode_integer_equation(primitive_integer_row(row)) for row in full_rref],
        "REDUCED_CANONICAL",
        "reduced integer equations are not the canonical exact RREF normalization",
    )

    fixed_components = fixed_line_components(raw)
    fixed_rref, fixed_pivots = q_rref(fixed_components)
    fixed_artifact = system["fixed_line_reduction"]
    fixed_x = [decode_integer_equation(row) for row in fixed_artifact["contracted_integer_equations_in_x"]]
    check(q_rank(fixed_x) == 11 and q_rank([*fixed_components, *fixed_x]) == 11, "FIXED_ROWSPACE", "contracted integer rowspace differs")
    check(
        fixed_artifact["contracted_integer_equations_in_x"]
        == [encode_integer_equation(primitive_integer_row(row)) for row in fixed_rref],
        "FIXED_CANONICAL",
        "contracted equations are not the canonical exact RREF normalization",
    )
    residual_rref = choose_complement(fixed_rref, full_rref)
    residual_x = [decode_integer_equation(row) for row in fixed_artifact["residual_integer_equations_in_x"]]
    check(len(residual_x) == q_rank(residual_x) == 27, "RESIDUAL_X_RANK", "x-coordinate residual rank differs")
    check(q_rank([*fixed_x, *residual_x]) == q_rank([*reduced, *fixed_x, *residual_x]) == 38, "RESIDUAL_X_SPAN", "x-coordinate split differs")
    check(
        fixed_artifact["residual_integer_equations_in_x"]
        == [encode_integer_equation(primitive_integer_row(row)) for row in residual_rref],
        "RESIDUAL_CANONICAL",
        "residual complement is not the canonical exact selection",
    )

    channel_data = fixed_artifact["channel_decomposition"]
    transform, inverse_transform = verify_channel_transform(vectors, group, channel_data)
    full_y = [transform_quadratic(row, inverse_transform) for row in full_rref]
    fixed_y = [transform_quadratic(row, inverse_transform) for row in fixed_rref]
    residual_y_computed = [transform_quadratic(row, inverse_transform) for row in residual_x]
    residual_y_artifact = [decode_k_equation(row) for row in fixed_artifact["residual_equations_in_channel_variables"]]
    check(residual_y_artifact == residual_y_computed, "RESIDUAL_TRANSFORM", "residual channel equations differ")
    products = [decode_k_equation(row["equation"]) for row in fixed_artifact["derivation_product_equations_in_channel_variables"]]
    check(products == expected_products(), "PRODUCT_EQUATIONS", "simple derivation products differ")
    validate_k_decomposition(full_y, fixed_y, products, residual_y_artifact)

    pure_d = [MONOMIAL_INDEX[(a, b)] for a in range(3) for b in range(a, 3)]
    check(all(not row[index] for row in residual_y_artifact for index in pure_d), "DERIVATION_ONLY", "residual does not vanish on d-only family")
    central_columns = [MONOMIAL_INDEX[(a, b)] for a in range(3, P) for b in range(a, P)]
    central_rank = k_rank([tuple(row[index] for index in central_columns) for row in residual_y_artifact])
    check(central_rank == fixed_artifact["central_specialization"]["residual_coefficient_rank_over_Qsqrt5"], "CENTRAL_RANK", "central specialization rank differs")

    channel_ids = [spec[0] for spec in CHANNEL_SPECS]
    expected_flats = recompute_flats(channel_ids)
    arrangement = fixed_artifact["derivation_weight_arrangement"]
    check(expected_flats == arrangement["relative_open_flats"], "FLAT_ARTIFACT", "arrangement flats differ")
    check(arrangement["flat_dimension_counts"] == {"0": 1, "1": 18, "2": 8, "3": 1}, "FLAT_COUNTS", "flat dimension counts differ")

    imports = audit_imports()
    validate_import_receipt(receipt, imports)

    validate_system_claims(
        system,
        basis_raw,
        stage1_receipt,
        group,
        raw,
        orbits,
        representatives,
        full_rref,
        full_pivots,
        fixed_components,
        fixed_rref,
        fixed_pivots,
        residual_rref,
        products,
        residual_y_artifact,
        central_rank,
        expected_flats,
    )
    validate_receipt_claims(system, receipt, imports)

    def validate_artifact_pair(changed_system: Mapping[str, Any], changed_receipt: Mapping[str, Any]) -> None:
        validate_stage2_keysets(changed_system, changed_receipt)
        validate_receipt_claims(changed_system, changed_receipt, imports)
        validate_system_claims(
            changed_system,
            basis_raw,
            stage1_receipt,
            group,
            raw,
            orbits,
            representatives,
            full_rref,
            full_pivots,
            fixed_components,
            fixed_rref,
            fixed_pivots,
            residual_rref,
            products,
            residual_y_artifact,
            central_rank,
            expected_flats,
        )
        verify_channel_transform(vectors, group, changed_system["fixed_line_reduction"]["channel_decomposition"])

    mutations = independent_mutations(
        system,
        receipt,
        representatives,
        reduced,
        full_y,
        fixed_y,
        products,
        residual_y_artifact,
        transform,
        inverse_transform,
        validate_artifact_pair,
    )
    return {
        "verified": True,
        "raw_jacobi_coordinates": len(raw),
        "target_orbits": len(orbits),
        "jacobi_rank": q_rank(reduced),
        "fixed_line_rank": q_rank(fixed_x),
        "residual_rank": k_rank(residual_y_artifact),
        "arrangement_flats": len(expected_flats),
        "central_specialization_residual_rank": central_rank,
        "independent_mutations": mutations,
        "full_classification_open": not receipt["not_attained"]["full_jacobi_solution_variety_classified"],
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), sort_keys=True))
    except VerificationError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
