#!/usr/bin/env python3
"""Exact target-free Jacobi reduction for the issue #566 stage-one basis.

This producer accepts exactly two semantic inputs: the stage-one rational
Reynolds basis and its receipt.  It constructs the Jacobi tensor of the most
general bracket in that fourteen-dimensional space.  No preferred coefficient,
external fixture, measurement, or physical target is accepted.

The attained reduction is deliberately algebraic rather than interpretive:

* 2,640 coordinate Jacobi quadrics reduce by equivariance to 44 representatives;
* those representatives span an exact rank-38 quadratic system over Q;
* contraction with the invariant port sum has rank 11;
* after extending scalars to Q(sqrt(5)), those eleven equations become explicit
  derivation-weight times channel-coordinate products;
* a certified rank-27 complement is retained without claiming its solution set
  has been classified.
"""

from __future__ import annotations

import argparse
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
from typing import Any, Iterable, Mapping, Sequence


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
STAGE1 = REPO / "code/a5_closure/issue_566_bracket_space_stage1"
STAGE1_BASIS = STAGE1 / "a5_alternating_bracket_reynolds_basis.json"
STAGE1_RECEIPT = STAGE1 / "a5_alternating_bracket_space_stage1.receipt.json"
SYSTEM_PATH = HERE / "a5_jacobi_system_reduction.json"
RECEIPT_PATH = HERE / "a5_jacobi_stage2.receipt.json"
VERIFIER_PATH = HERE / "verify.py"
TEST_PATH = HERE / "test_stage2.py"

PORTS = 12
PARAMETERS = 14
MONOMIALS = [(a, b) for a in range(PARAMETERS) for b in range(a, PARAMETERS)]
MONOMIAL_INDEX = {monomial: i for i, monomial in enumerate(MONOMIALS)}
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
EXPECTED_STAGE1_BASIS_KEYS = {
    "ambient_coordinate_count",
    "basis",
    "coordinate_convention",
    "domain_dimension",
    "field",
    "group_action_convention",
    "issue",
    "normalization",
    "port_dimension",
    "proper_port_group_sha256",
    "schema",
}
EXPECTED_STAGE1_RECEIPT_KEYS = {
    "claim_boundary",
    "dimension_certificate",
    "implementation_pins",
    "issue",
    "later_gates",
    "mutation_tests",
    "proper_port_action",
    "receipt_sha256",
    "representation",
    "reynolds_basis",
    "schema",
    "semantic_input",
    "status",
    "target_firewall",
}
STAGE1_BASIS_ROW_KEYS = {
    "basis_id",
    "entries",
    "integral_orbit_scale",
    "orbit_size",
    "seed_coordinate",
    "stabilizer_size",
}
STAGE1_SEMANTIC_KEYS = {
    "accessed_value_paths",
    "carrier_projection_sha256",
    "file_sha256",
    "ignored_value_sections",
    "path",
}
STAGE1_FIREWALL_KEYS = {
    "arbitrary_input_path_allowed",
    "audited_import_roots",
    "audited_import_roots_by_file",
    "basis_seed_rule",
    "command_line_target_parameters_allowed",
    "desired_coefficient_input_allowed",
    "enabled",
    "environment_target_read",
    "exact_manifest_keysets_enforced",
    "measurement_input_allowed",
    "network_target_read",
    "nonstdlib_imports",
    "semantic_input_count",
    "target_free_scope",
}
STAGE1_ACTION_KEYS = {
    "a5_isomorphism",
    "all_rows_bijective",
    "closure_checked_pairs",
    "identity_present",
    "incidence_automorphism_count",
    "incidence_preserved",
    "inverses_present",
    "orientation_preserving_count",
    "orientation_reversing_count",
    "oriented_faces_preserved",
    "permutation_rows",
    "permutation_rows_sha256",
    "port_count",
    "transitive_on_ports",
}
STAGE1_A5_KEYS = {
    "bijection_size",
    "generator_orders",
    "homomorphism_checked_pairs",
    "isomorphic",
    "mapping_sha256",
    "model",
    "model_generators",
    "model_order",
    "port_generators",
}
STAGE1_REPRESENTATION_KEYS = {
    "V_dimension",
    "action",
    "exterior_square_V_dimension",
    "field",
    "hom_ambient_dimension",
}
STAGE1_DIMENSION_KEYS = {
    "class_rows",
    "dimension",
    "formula",
    "projector_trace_denominator",
    "projector_trace_numerator",
}
STAGE1_CLASS_ROW_KEYS = {
    "chi_V",
    "chi_V_of_square",
    "chi_exterior_square_V",
    "count",
    "element_order",
    "projector_trace_contribution_each",
    "projector_trace_contribution_total",
}
STAGE1_REYNOLDS_KEYS = {
    "all_signed_orbits_orientable",
    "alternation_checks",
    "basis_count",
    "canonical_json_sha256",
    "covariance_checks",
    "exact_alternation_passed",
    "exact_covariance_passed",
    "normalization",
    "orbit_sizes",
    "path",
    "rank_certificate",
    "support_union_size",
    "supports_pairwise_disjoint",
}
STAGE1_RANK_KEYS = {
    "ambient_coordinate_count",
    "gram_determinant",
    "gram_diagonal",
    "matrix_shape",
    "pivot_coordinates",
    "pivot_minor_determinant",
    "rank_over_Q",
}
STAGE1_FRACTION_KEYS = {"denominator", "numerator"}
STAGE1_PIN_KEYS = {"independent_verifier_sha256", "producer_sha256", "test_sha256"}
MUTATION_ROW_KEYS = {"actual_error", "expected_error", "name", "passed"}

SYSTEM_KEYS = {"field", "fixed_line_reduction", "issue", "jacobi_tensor", "monomial_order", "parameter_variables", "schema", "upstream"}
UPSTREAM_KEYS = {
    "basis_canonical_json_sha256",
    "basis_file_sha256",
    "basis_path",
    "proper_group_sha256",
    "receipt_file_sha256",
    "receipt_path",
    "receipt_self_sha256",
}
JACOBI_KEYS = {
    "all_coordinate_equations_nonzero",
    "all_coordinate_equations_sha256",
    "convention",
    "coordinate_equation_count",
    "reduced_integer_equation_count",
    "reduced_integer_equations",
    "representative_coefficient_rank_over_Q",
    "representative_equations",
    "representative_rref_pivot_monomial_indices",
    "signed_target_orbit_count",
    "signed_target_orbit_sizes",
    "target_invariant_character_certificate",
    "target_space_dimension",
    "zero_locus_equivalence",
}
TARGET_CHARACTER_KEYS = {"class_rows", "dimension", "projector_trace_denominator", "projector_trace_numerator"}
TARGET_CLASS_ROW_KEYS = {
    "chi_V",
    "chi_V_of_cube",
    "chi_V_of_square",
    "chi_exterior_cube_V",
    "contribution_each",
    "contribution_total",
    "count",
    "element_order",
}
REPRESENTATIVE_KEYS = {"equation", "signed_orbit_size", "target_coordinate"}
FIXED_REDUCTION_KEYS = {
    "central_specialization",
    "channel_decomposition",
    "contracted_coefficient_rank_over_Q",
    "contracted_coordinate_count",
    "contracted_integer_equations_in_x",
    "contracted_rref_pivot_monomial_indices",
    "derivation_only_subspace",
    "derivation_product_equations_in_channel_variables",
    "derivation_weight_arrangement",
    "fixed_vector",
    "full_rank_split",
    "product_equation_rank_over_Qsqrt5",
    "product_rowspace_equals_contracted_jacobi_rowspace",
    "residual_equations_in_channel_variables",
    "residual_integer_equations_in_x",
    "residual_rank_over_Qsqrt5",
}
CHANNEL_DECOMPOSITION_KEYS = {
    "channel_variables",
    "inverse_transform_rows",
    "spectral_decomposition",
    "splitting_field",
    "transform_determinant",
    "transform_rows",
}
D_CHANNEL_KEYS = {"id", "linear_form_in_x", "meaning"}
T_CHANNEL_KEYS = {"derivation_weight_in_d", "domain_sectors", "extraction_coordinate", "id", "linear_form_in_x", "output_sector"}
SPECTRAL_KEYS = {"canonical_valency_five_orbital", "minimal_polynomial", "operator_sha256", "sectors"}
SECTOR_KEYS = {"dimension", "id", "operator_eigenvalue"}
PRODUCT_DECLARATION_KEYS = {"channel_id", "equation", "weight_coefficients_on_d_plus_d_minus_d_five"}
DERIVATION_ONLY_KEYS = {"all_jacobi_equations_vanish_identically", "definition", "dimension", "scope"}
CENTRAL_SPECIALIZATION_KEYS = {"definition", "remaining_channel_variable_count", "residual_coefficient_rank_over_Qsqrt5", "solution_set_classified"}
ARRANGEMENT_KEYS = {"distinct_weight_hyperplane_count", "flat_count", "flat_dimension_counts", "relative_open_flats", "scope"}
FLAT_KEYS = {
    "channels_forced_zero_on_relative_open_flat",
    "channels_not_forced_zero_by_fixed_line_sector",
    "defining_rowspace_rref",
    "dimension",
    "flat_id",
    "scope",
    "weights_vanishing_identically",
}
RECEIPT_KEYS = {
    "attained",
    "claim_boundary",
    "implementation_pins",
    "issue",
    "mutation_tests",
    "not_attained",
    "receipt_sha256",
    "schema",
    "semantic_inputs",
    "status",
    "system_artifact",
    "target_firewall",
}
SEMANTIC_INPUT_KEYS = {"count", "paths", "stage1_later_gates_all_false"}
FIREWALL_KEYS = {
    "accepted_input_schemas",
    "arbitrary_input_path_allowed",
    "audited_import_roots_by_file",
    "coefficient_target_input_allowed",
    "enabled",
    "environment_target_read",
    "exact_import_allowlists_by_file",
    "exact_input_keysets_enforced",
    "measurement_input_allowed",
    "network_target_read",
    "nonstdlib_imports",
    "runtime_forbidden_import_roots_by_file",
}
ATTAINED_KEYS = {
    "derivation_only_three_parameter_family_certified",
    "derivation_weight_relative_open_flat_count",
    "equivariant_target_orbit_count",
    "exact_channel_coordinate_change",
    "fixed_line_product_quadrics",
    "general_parameter_count",
    "generic_weight_stratum_forces_all_non_derivation_channels_zero",
    "independent_jacobi_quadrics",
    "raw_jacobi_coordinate_count",
    "residual_independent_quadrics",
}
NOT_ATTAINED_KEYS = {
    "compactness_established",
    "full_jacobi_solution_variety_classified",
    "issue_566_closed",
    "preferred_bracket_selected",
    "rational_descent_components_classified",
    "source_selection",
    "special_weight_flat_residual_solutions_classified",
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
Coordinate2 = tuple[int, int, int]
Coordinate3 = tuple[int, int, int, int]
QRow = tuple[Fraction, ...]


class CertificateError(RuntimeError):
    def __init__(self, code: str, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


def require_keys(value: Any, expected: set[str], code: str) -> None:
    require(isinstance(value, Mapping), code, "value is not an object")
    require(set(value) == expected, code, f"keyset differs: {sorted(set(value) ^ expected)}")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def object_sha256(value: Any) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def compose(p: Permutation, q: Permutation) -> Permutation:
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p: Permutation) -> Permutation:
    result = [0] * len(p)
    for i, image in enumerate(p):
        result[image] = i
    return tuple(result)


def parity(values: Sequence[int]) -> int:
    return sum(values[i] > values[j] for i in range(len(values)) for j in range(i + 1, len(values))) % 2


def invoke_hardened_stage1_verifier() -> None:
    try:
        namespace = runpy.run_path(str(STAGE1 / "verify.py"), run_name="_oph_stage1_hardened_verify")
        result = namespace["verify"]()
    except Exception as exc:
        raise CertificateError("UPSTREAM_HARDENED_VERIFY", f"stage-one independent verifier failed: {exc}") from exc
    require(
        isinstance(result, Mapping)
        and result.get("verified") is True
        and result.get("basis_rank") == PARAMETERS
        and result.get("later_gates_all_false") is True,
        "UPSTREAM_HARDENED_VERIFY",
        "stage-one independent verifier returned an unexpected result",
    )


def load_stage1() -> tuple[dict[str, Any], dict[str, Any], list[dict[Coordinate2, Fraction]], list[Permutation]]:
    invoke_hardened_stage1_verifier()
    basis_raw = json.loads(STAGE1_BASIS.read_text(encoding="utf-8"))
    receipt_raw = json.loads(STAGE1_RECEIPT.read_text(encoding="utf-8"))
    require(isinstance(basis_raw, dict) and isinstance(receipt_raw, dict), "UPSTREAM_JSON", "stage-one artifacts must be objects")
    validate_upstream_keysets(basis_raw, receipt_raw)
    require(
        basis_raw.get("schema") == "oph.a5_alternating_bracket_reynolds_basis.v1",
        "UPSTREAM_BASIS_SCHEMA",
        "unexpected stage-one basis schema",
    )
    require(
        receipt_raw.get("schema") == "oph.a5_alternating_bracket_space_stage1.receipt.v1",
        "UPSTREAM_RECEIPT_SCHEMA",
        "unexpected stage-one receipt schema",
    )
    receipt_copy = dict(receipt_raw)
    stored_self_hash = receipt_copy.pop("receipt_sha256", None)
    require(stored_self_hash == object_sha256(receipt_copy), "UPSTREAM_RECEIPT_HASH", "stage-one receipt self-hash differs")
    require(receipt_raw.get("issue") == 566 and basis_raw.get("issue") == 566, "UPSTREAM_ISSUE", "wrong upstream issue")
    require(not any(receipt_raw.get("later_gates", {}).values()), "UPSTREAM_GATE", "stage-one later gate was promoted")
    require(
        object_sha256(basis_raw) == receipt_raw["reynolds_basis"]["canonical_json_sha256"],
        "UPSTREAM_BASIS_HASH",
        "stage-one basis hash differs",
    )

    group_rows = receipt_raw["proper_port_action"].get("permutation_rows")
    require(isinstance(group_rows, list) and len(group_rows) == 60, "UPSTREAM_GROUP", "stage-one group rows are absent")
    group = [tuple(row) for row in group_rows]
    require(len(set(group)) == 60, "UPSTREAM_GROUP", "stage-one group rows are duplicated")
    group_set = set(group)
    identity = tuple(range(PORTS))
    require(identity in group_set, "UPSTREAM_GROUP", "identity row is absent")
    for row in group:
        require(len(row) == PORTS and set(row) == set(range(PORTS)), "UPSTREAM_GROUP", "group row is not a bijection")
        require(inverse(row) in group_set, "UPSTREAM_GROUP", "group inverse is absent")
        for other in group:
            require(compose(row, other) in group_set, "UPSTREAM_GROUP", "group closure fails")
    require(
        object_sha256([list(row) for row in sorted(group)])
        == receipt_raw["proper_port_action"]["permutation_rows_sha256"],
        "UPSTREAM_GROUP_HASH",
        "stage-one group hash differs",
    )

    rows = basis_raw.get("basis")
    require(isinstance(rows, list) and len(rows) == PARAMETERS, "UPSTREAM_BASIS_COUNT", "expected fourteen basis rows")
    vectors: list[dict[Coordinate2, Fraction]] = []
    valid_coordinates = {
        (output, left, right)
        for output in range(PORTS)
        for left in range(PORTS)
        for right in range(left + 1, PORTS)
    }
    for number, row in enumerate(rows):
        require(isinstance(row, Mapping) and row.get("basis_id") == f"R{number:02d}", "UPSTREAM_BASIS_ROW", "basis row id differs")
        vector: dict[Coordinate2, Fraction] = {}
        for entry in row.get("entries", []):
            require(isinstance(entry, list) and len(entry) == 5 and all(type(x) is int for x in entry), "UPSTREAM_BASIS_ENTRY", "bad basis entry")
            coordinate = tuple(entry[:3])
            require(coordinate in valid_coordinates and coordinate not in vector, "UPSTREAM_BASIS_ENTRY", "bad or duplicate basis coordinate")
            vector[coordinate] = Fraction(entry[3], entry[4])
        require(len(vector) == row.get("orbit_size"), "UPSTREAM_BASIS_ROW", "basis orbit size differs")
        vectors.append(vector)
    require(sum(len(vector) for vector in vectors) == 792, "UPSTREAM_BASIS_SUPPORT", "basis support count differs")

    # Recheck the mathematical covariance, not only the upstream hashes.
    for vector in vectors:
        for row in group:
            for coordinate in valid_coordinates:
                output, left, right = coordinate
                a, b = row[left], row[right]
                image = (row[output], min(a, b), max(a, b))
                sign = 1 if a < b else -1
                require(vector.get(image, Fraction(0)) == sign * vector.get(coordinate, Fraction(0)), "UPSTREAM_COVARIANCE", "stage-one basis covariance fails")
    validate_stage1_claims(basis_raw, receipt_raw, vectors, group)
    return basis_raw, receipt_raw, vectors, sorted(group)


def validate_upstream_keysets(basis_raw: Mapping[str, Any], receipt_raw: Mapping[str, Any]) -> None:
    require_keys(basis_raw, EXPECTED_STAGE1_BASIS_KEYS, "FIREWALL_BASIS_KEYS")
    basis_rows = basis_raw.get("basis")
    require(isinstance(basis_rows, list), "FIREWALL_BASIS_ROWS", "stage-one basis rows are not a list")
    for row in basis_rows:
        require_keys(row, STAGE1_BASIS_ROW_KEYS, "FIREWALL_BASIS_ROW_KEYS")

    require_keys(receipt_raw, EXPECTED_STAGE1_RECEIPT_KEYS, "FIREWALL_RECEIPT_KEYS")
    require_keys(receipt_raw.get("semantic_input"), STAGE1_SEMANTIC_KEYS, "FIREWALL_STAGE1_SEMANTIC_KEYS")
    firewall = receipt_raw.get("target_firewall")
    require_keys(firewall, STAGE1_FIREWALL_KEYS, "FIREWALL_STAGE1_FIREWALL_KEYS")
    audited = firewall.get("audited_import_roots_by_file")
    require_keys(audited, {"certify.py", "test_stage1.py", "verify.py"}, "FIREWALL_STAGE1_IMPORT_FILES")
    require_keys(receipt_raw.get("proper_port_action"), STAGE1_ACTION_KEYS, "FIREWALL_STAGE1_ACTION_KEYS")
    require_keys(receipt_raw["proper_port_action"].get("a5_isomorphism"), STAGE1_A5_KEYS, "FIREWALL_STAGE1_A5_KEYS")
    require_keys(receipt_raw.get("representation"), STAGE1_REPRESENTATION_KEYS, "FIREWALL_STAGE1_REPRESENTATION_KEYS")
    dimension = receipt_raw.get("dimension_certificate")
    require_keys(dimension, STAGE1_DIMENSION_KEYS, "FIREWALL_STAGE1_DIMENSION_KEYS")
    require(isinstance(dimension.get("class_rows"), list), "FIREWALL_STAGE1_CLASS_ROWS", "stage-one class rows are absent")
    for row in dimension["class_rows"]:
        require_keys(row, STAGE1_CLASS_ROW_KEYS, "FIREWALL_STAGE1_CLASS_ROW_KEYS")
    reynolds = receipt_raw.get("reynolds_basis")
    require_keys(reynolds, STAGE1_REYNOLDS_KEYS, "FIREWALL_STAGE1_REYNOLDS_KEYS")
    rank = reynolds.get("rank_certificate")
    require_keys(rank, STAGE1_RANK_KEYS, "FIREWALL_STAGE1_RANK_KEYS")
    require_keys(rank.get("pivot_minor_determinant"), STAGE1_FRACTION_KEYS, "FIREWALL_STAGE1_FRACTION_KEYS")
    require_keys(rank.get("gram_determinant"), STAGE1_FRACTION_KEYS, "FIREWALL_STAGE1_FRACTION_KEYS")
    require(isinstance(rank.get("gram_diagonal"), list), "FIREWALL_STAGE1_GRAM_ROWS", "stage-one Gram diagonal is absent")
    for row in rank["gram_diagonal"]:
        require_keys(row, STAGE1_FRACTION_KEYS, "FIREWALL_STAGE1_FRACTION_KEYS")
    later_keys = {
        "bracket_selected",
        "carrier_choice_derived_by_this_packet",
        "compactness_established",
        "issue_566_closed",
        "jacobi_solution_classification_complete",
        "jacobi_solution_found",
        "physical_current_identified",
        "source_selection",
    }
    require_keys(receipt_raw.get("later_gates"), later_keys, "FIREWALL_STAGE1_GATE_KEYS")
    require_keys(receipt_raw.get("implementation_pins"), STAGE1_PIN_KEYS, "FIREWALL_STAGE1_PIN_KEYS")
    mutations = receipt_raw.get("mutation_tests")
    require(isinstance(mutations, list), "FIREWALL_STAGE1_MUTATION_ROWS", "stage-one mutations are not a list")
    for row in mutations:
        require_keys(row, MUTATION_ROW_KEYS, "FIREWALL_STAGE1_MUTATION_KEYS")


def stage1_character_claim(group: Sequence[Permutation]) -> dict[str, Any]:
    bins: dict[tuple[int, int, int], int] = {}
    numerator = 0
    identity = tuple(range(PORTS))
    for row in group:
        power = identity
        order = 0
        while True:
            order += 1
            power = compose(row, power)
            if power == identity:
                break
        square = compose(row, row)
        chi = sum(row[i] == i for i in range(PORTS))
        chi2 = sum(square[i] == i for i in range(PORTS))
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
    basis_raw: Mapping[str, Any],
    receipt_raw: Mapping[str, Any],
    vectors: Sequence[Mapping[Coordinate2, Fraction]],
    group: Sequence[Permutation],
) -> None:
    expected_basis_metadata = {
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
    for key, expected in expected_basis_metadata.items():
        require(basis_raw.get(key) == expected, "UPSTREAM_BASIS_CLAIM", f"stage-one basis {key} differs")
    action = receipt_raw["proper_port_action"]
    group_hash = object_sha256([list(row) for row in group])
    require(basis_raw.get("proper_port_group_sha256") == group_hash, "UPSTREAM_GROUP_HASH", "basis group pin differs")
    require(action.get("permutation_rows") == [list(row) for row in group], "UPSTREAM_GROUP_ROWS", "stage-one group rows differ")
    require(action.get("permutation_rows_sha256") == group_hash, "UPSTREAM_GROUP_HASH", "stage-one group hash differs")
    expected_action = {
        "port_count": 12,
        "incidence_automorphism_count": 120,
        "orientation_preserving_count": 60,
        "orientation_reversing_count": 60,
        "identity_present": True,
        "all_rows_bijective": True,
        "incidence_preserved": True,
        "oriented_faces_preserved": True,
        "closure_checked_pairs": 3600,
        "inverses_present": True,
        "transitive_on_ports": True,
    }
    for key, expected in expected_action.items():
        require(action.get(key) == expected, "UPSTREAM_GROUP_CLAIM", f"stage-one action {key} differs")
    a5 = action["a5_isomorphism"]
    require(
        a5.get("model") == "even permutations on five unlabeled symbols"
        and a5.get("model_order") == 60
        and a5.get("bijection_size") == 60
        and a5.get("generator_orders") == [2, 3, 5]
        and a5.get("homomorphism_checked_pairs") == 3600
        and a5.get("isomorphic") is True,
        "UPSTREAM_A5_CLAIM",
        "stage-one A5 witness declaration differs",
    )
    require(all(tuple(row) in set(group) for row in a5.get("port_generators", [])), "UPSTREAM_A5_CLAIM", "A5 port generator is absent")

    require(receipt_raw.get("issue") == 566, "UPSTREAM_ISSUE", "wrong stage-one issue")
    require(receipt_raw.get("status") == "EXACT_TARGET_FREE_SEARCH_SPACE_STAGE1_ONLY", "UPSTREAM_STATUS", "stage-one status differs")
    require(
        receipt_raw.get("claim_boundary")
        == "Conditional on the pinned canonical oriented twelve-port carrier, this receipt certifies only the complete equivariant alternating-bracket search space on its permutation module. Target-free means that no desired gauge algebra, measurement, or coefficient target enters this stage; it does not mean that this packet derives the carrier choice. The receipt does not select a bracket, solve Jacobi, establish compactness, identify a physical current, or close issue #566.",
        "UPSTREAM_SCOPE",
        "stage-one claim boundary differs",
    )
    expected_representation = {
        "field": "Q",
        "V_dimension": 12,
        "exterior_square_V_dimension": 66,
        "hom_ambient_dimension": 792,
        "action": "simultaneous output and signed unordered-input permutation",
    }
    require(receipt_raw.get("representation") == expected_representation, "UPSTREAM_REPRESENTATION", "stage-one representation claim differs")
    require(receipt_raw.get("dimension_certificate") == stage1_character_claim(group), "UPSTREAM_DIMENSION", "stage-one dimension certificate differs")

    later = receipt_raw["later_gates"]
    require(not any(later.values()), "UPSTREAM_GATE", "stage-one gate was promoted")
    semantic = receipt_raw["semantic_input"]
    require(
        semantic.get("path") == "code/a5_closure/manifests/echosahedral_federation_reference.json"
        and semantic.get("accessed_value_paths") == ["/schema", "/carrier/ports", "/carrier/edges", "/carrier/oriented_faces"]
        and semantic.get("ignored_value_sections")
        == [
            "/architecture",
            "/carrier/atoms_pairwise_orthogonal",
            "/carrier/atoms_sum_to_one",
            "/carrier/central_port_atoms",
            "/refinement_tower",
            "/source_readback",
        ],
        "UPSTREAM_SEMANTIC_SCOPE",
        "stage-one semantic-input scope differs",
    )
    firewall = receipt_raw["target_firewall"]
    require(
        firewall.get("enabled") is True
        and firewall.get("semantic_input_count") == 1
        and firewall.get("arbitrary_input_path_allowed") is False
        and firewall.get("command_line_target_parameters_allowed") is False
        and firewall.get("desired_coefficient_input_allowed") is False
        and firewall.get("measurement_input_allowed") is False
        and firewall.get("environment_target_read") is False
        and firewall.get("network_target_read") is False
        and firewall.get("nonstdlib_imports") == []
        and firewall.get("exact_manifest_keysets_enforced") is True
        and firewall.get("basis_seed_rule") == "lexicographically first unassigned tensor coordinate"
        and firewall.get("target_free_scope") == "conditional_on_pinned_canonical_oriented_twelve_port_carrier",
        "UPSTREAM_FIREWALL",
        "stage-one firewall declaration differs",
    )

    coordinates = sorted({coordinate for vector in vectors for coordinate in vector})
    require(len(coordinates) == 792, "UPSTREAM_BASIS_SUPPORT", "stage-one support union differs")
    supports = [set(vector) for vector in vectors]
    require(sum(map(len, supports)) == len(set().union(*supports)), "UPSTREAM_BASIS_SUPPORT", "stage-one supports overlap")
    matrix = [[vector.get(coordinate, Fraction(0)) for coordinate in coordinates] for vector in vectors]
    basis_rank = rational_rank(matrix)
    require(basis_rank == 14, "UPSTREAM_BASIS_RANK", "stage-one basis rank differs")
    reynolds = receipt_raw["reynolds_basis"]
    require(
        reynolds.get("path") == "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
        and reynolds.get("canonical_json_sha256") == object_sha256(basis_raw)
        and reynolds.get("basis_count") == 14
        and reynolds.get("orbit_sizes") == sorted(map(len, vectors))
        and reynolds.get("support_union_size") == 792
        and reynolds.get("supports_pairwise_disjoint") is True
        and reynolds.get("all_signed_orbits_orientable") is True
        and reynolds.get("normalization") == "R(seed)=(1/60)*sum_g g(seed); nonzero entries are signed reciprocals of orbit size"
        and reynolds.get("alternation_checks") == 24192
        and reynolds.get("exact_alternation_passed") is True
        and reynolds.get("covariance_checks") == 665280
        and reynolds.get("exact_covariance_passed") is True,
        "UPSTREAM_REYNOLDS_CLAIM",
        "stage-one Reynolds receipt differs",
    )
    rank = reynolds["rank_certificate"]
    gram = [sum((value * value for value in vector.values()), Fraction(0)) for vector in vectors]
    gram_product = math.prod(gram, start=Fraction(1))
    decoded_gram = [Fraction(row["numerator"], row["denominator"]) for row in rank["gram_diagonal"]]
    require(
        rank.get("ambient_coordinate_count") == 792
        and rank.get("matrix_shape") == [792, 14]
        and rank.get("rank_over_Q") == basis_rank
        and rank.get("pivot_coordinates") == [list(min(vector)) for vector in vectors]
        and decoded_gram == gram
        and Fraction(rank["gram_determinant"]["numerator"], rank["gram_determinant"]["denominator"]) == gram_product
        and Fraction(rank["pivot_minor_determinant"]["numerator"], rank["pivot_minor_determinant"]["denominator"]) == gram_product,
        "UPSTREAM_RANK_CLAIM",
        "stage-one rank certificate differs",
    )
    pins = receipt_raw["implementation_pins"]
    require(
        pins.get("producer_sha256") == file_sha256(STAGE1 / "certify.py")
        and pins.get("independent_verifier_sha256") == file_sha256(STAGE1 / "verify.py")
        and pins.get("test_sha256") == file_sha256(STAGE1 / "test_stage1.py"),
        "UPSTREAM_IMPLEMENTATION_PIN",
        "stage-one implementation pin differs",
    )
    require(all(row.get("passed") is True for row in receipt_raw["mutation_tests"]), "UPSTREAM_MUTATION", "stage-one mutation failed")


def full_structure_constants(vectors: Sequence[Mapping[Coordinate2, Fraction]]) -> list[list[list[list[Fraction]]]]:
    constants = [
        [[[Fraction(0) for _ in range(PORTS)] for _ in range(PORTS)] for _ in range(PORTS)]
        for _ in range(PARAMETERS)
    ]
    for parameter, vector in enumerate(vectors):
        for (output, left, right), value in vector.items():
            constants[parameter][output][left][right] = value
            constants[parameter][output][right][left] = -value
    return constants


def jacobi_component(
    constants: Sequence[Sequence[Sequence[Sequence[Fraction]]]],
    output: int,
    i: int,
    j: int,
    k: int,
) -> QRow:
    result = [Fraction(0)] * len(MONOMIALS)
    for middle in range(PORTS):
        for a in range(PARAMETERS):
            first = constants[a][middle][i][j]
            second = constants[a][middle][j][k]
            third = constants[a][middle][k][i]
            if first == second == third == 0:
                continue
            for b in range(PARAMETERS):
                value = (
                    first * constants[b][output][middle][k]
                    + second * constants[b][output][middle][i]
                    + third * constants[b][output][middle][j]
                )
                if value:
                    result[MONOMIAL_INDEX[tuple(sorted((a, b)))]] += value
    return tuple(result)


def build_raw_jacobi(
    constants: Sequence[Sequence[Sequence[Sequence[Fraction]]]],
) -> dict[Coordinate3, QRow]:
    result: dict[Coordinate3, QRow] = {}
    for i, j, k in itertools.combinations(range(PORTS), 3):
        for output in range(PORTS):
            row = jacobi_component(constants, output, i, j, k)
            require(any(row), "JACOBI_COMPONENT_ZERO", "unexpected identically zero coordinate equation")
            result[(output, i, j, k)] = row
    require(len(result) == 2640, "JACOBI_COMPONENT_COUNT", "wrong Jacobi coordinate count")
    return result


def target_action(row: Permutation, coordinate: Coordinate3) -> tuple[Coordinate3, int]:
    output, i, j, k = coordinate
    images = [row[i], row[j], row[k]]
    sign = -1 if parity(images) else 1
    return (row[output], *sorted(images)), sign


def target_orbits(group: Sequence[Permutation]) -> list[dict[Coordinate3, int]]:
    unassigned = {
        (output, i, j, k)
        for output in range(PORTS)
        for i, j, k in itertools.combinations(range(PORTS), 3)
    }
    result: list[dict[Coordinate3, int]] = []
    while unassigned:
        seed = min(unassigned)
        signs: dict[Coordinate3, int] = {}
        for row in group:
            image, sign = target_action(row, seed)
            require(image not in signs or signs[image] == sign, "TARGET_SIGNED_ORBIT", "signed target orbit has a negative stabilizer")
            signs[image] = sign
        unassigned -= set(signs)
        result.append(signs)
    require(len(result) == 44 and all(len(orbit) == 60 for orbit in result), "TARGET_ORBIT_COUNT", "expected forty-four size-sixty target orbits")
    return result


def scalar_row(row: QRow, scalar: Fraction) -> QRow:
    return tuple(scalar * value for value in row)


def add_rows(left: QRow, right: QRow) -> QRow:
    return tuple(a + b for a, b in zip(left, right))


def character_target(group: Sequence[Permutation]) -> dict[str, Any]:
    bins: dict[tuple[int, int, int, int], int] = {}
    numerator = 0
    for row in group:
        identity = tuple(range(PORTS))
        power = identity
        order = 0
        while True:
            order += 1
            power = compose(row, power)
            if power == identity:
                break
        square = compose(row, row)
        cube = compose(square, row)
        chi = sum(row[i] == i for i in range(PORTS))
        chi2 = sum(square[i] == i for i in range(PORTS))
        chi3 = sum(cube[i] == i for i in range(PORTS))
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
    require(numerator == 2640 and numerator // 60 == 44, "TARGET_CHARACTER", "target invariant dimension is not forty-four")
    return {"class_rows": rows, "projector_trace_numerator": numerator, "projector_trace_denominator": 60, "dimension": 44}


def rref(rows: Sequence[Sequence[Fraction]]) -> tuple[list[QRow], list[int]]:
    work = [list(row) for row in rows if any(row)]
    if not work:
        return [], []
    pivot_row = 0
    pivots: list[int] = []
    for column in range(len(work[0])):
        pivot = next((r for r in range(pivot_row, len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for r in range(len(work)):
            if r != pivot_row and work[r][column]:
                factor = work[r][column]
                work[r] = [x - factor * y for x, y in zip(work[r], work[pivot_row])]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(work):
            break
    return [tuple(row) for row in work[:pivot_row]], pivots


def rational_rank(rows: Sequence[Sequence[Fraction]]) -> int:
    return len(rref(rows)[0])


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
    require(divisor > 0, "PRIMITIVE_ZERO", "cannot normalize a zero row")
    integers = [value // divisor for value in integers]
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def encode_rational_equation(row: Sequence[Fraction]) -> list[list[int]]:
    return [
        [a, b, value.numerator, value.denominator]
        for (a, b), value in zip(MONOMIALS, row)
        if value
    ]


def encode_integer_equation(row: Sequence[int]) -> list[list[int]]:
    return [[a, b, value] for (a, b), value in zip(MONOMIALS, row) if value]


def raw_component_hash(raw: Mapping[Coordinate3, QRow]) -> str:
    encoded = [
        [list(coordinate), encode_rational_equation(raw[coordinate])]
        for coordinate in sorted(raw)
    ]
    return object_sha256(encoded)


def signed_raw_component(raw: Mapping[Coordinate3, QRow], output: int, values: Sequence[int]) -> QRow:
    if len(set(values)) < 3:
        return tuple(Fraction(0) for _ in MONOMIALS)
    ordered = sorted(values)
    row = raw[(output, *ordered)]
    return scalar_row(row, Fraction(-1 if parity(values) else 1))


def invariant_line_rows(raw: Mapping[Coordinate3, QRow]) -> list[QRow]:
    rows: list[QRow] = []
    for output in range(PORTS):
        for i, j in itertools.combinations(range(PORTS), 2):
            total = tuple(Fraction(0) for _ in MONOMIALS)
            for source in range(PORTS):
                total = add_rows(total, signed_raw_component(raw, output, (source, i, j)))
            require(any(total), "FIXED_LINE_COMPONENT_ZERO", "unexpected zero fixed-line component")
            rows.append(total)
    require(len(rows) == 792 and rational_rank(rows) == 11, "FIXED_LINE_RANK", "fixed-line Jacobi sector does not have rank eleven")
    return rows


@dataclass(frozen=True)
class Q5:
    """The exact quadratic field Q(sqrt(5)), represented as a+b*sqrt(5)."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", Fraction(self.a))
        object.__setattr__(self, "b", Fraction(self.b))

    @staticmethod
    def make(value: Any) -> "Q5":
        if isinstance(value, Q5):
            return value
        return Q5(Fraction(value), Fraction(0))

    def __add__(self, other: Any) -> "Q5":
        value = Q5.make(other)
        return Q5(self.a + value.a, self.b + value.b)

    __radd__ = __add__

    def __neg__(self) -> "Q5":
        return Q5(-self.a, -self.b)

    def __sub__(self, other: Any) -> "Q5":
        return self + (-Q5.make(other))

    def __rsub__(self, other: Any) -> "Q5":
        return Q5.make(other) - self

    def __mul__(self, other: Any) -> "Q5":
        value = Q5.make(other)
        return Q5(self.a * value.a + 5 * self.b * value.b, self.a * value.b + self.b * value.a)

    __rmul__ = __mul__

    def inverse(self) -> "Q5":
        norm = self.a * self.a - 5 * self.b * self.b
        require(norm != 0, "Q5_DIVISION", "division by zero in Q(sqrt(5))")
        return Q5(self.a / norm, -self.b / norm)

    def __truediv__(self, other: Any) -> "Q5":
        return self * Q5.make(other).inverse()

    def __bool__(self) -> bool:
        return bool(self.a or self.b)

    def conjugate(self) -> "Q5":
        return Q5(self.a, -self.b)


Q5_ZERO = Q5()
Q5_ONE = Q5(Fraction(1), Fraction(0))
Q5_SQRT = Q5(Fraction(0), Fraction(1))
Q5Matrix = list[list[Q5]]
Q5Row = tuple[Q5, ...]


def q5_matrix(rows: int, columns: int, fill: Q5 = Q5_ZERO) -> Q5Matrix:
    return [[fill for _ in range(columns)] for _ in range(rows)]


def q5_identity(size: int) -> Q5Matrix:
    result = q5_matrix(size, size)
    for i in range(size):
        result[i][i] = Q5_ONE
    return result


def q5_matrix_add(left: Q5Matrix, right: Q5Matrix) -> Q5Matrix:
    return [[left[i][j] + right[i][j] for j in range(len(left[0]))] for i in range(len(left))]


def q5_matrix_sub(left: Q5Matrix, right: Q5Matrix) -> Q5Matrix:
    return [[left[i][j] - right[i][j] for j in range(len(left[0]))] for i in range(len(left))]


def q5_matrix_scale(scalar: Q5, matrix: Q5Matrix) -> Q5Matrix:
    return [[scalar * value for value in row] for row in matrix]


def q5_matrix_mul(left: Q5Matrix, right: Q5Matrix) -> Q5Matrix:
    rows, inner, columns = len(left), len(right), len(right[0])
    require(len(left[0]) == inner, "Q5_MATRIX_SHAPE", "matrix product shape mismatch")
    result = q5_matrix(rows, columns)
    for i in range(rows):
        for k in range(inner):
            if not left[i][k]:
                continue
            for j in range(columns):
                if right[k][j]:
                    result[i][j] = result[i][j] + left[i][k] * right[k][j]
    return result


def q5_matrix_equal(left: Q5Matrix, right: Q5Matrix) -> bool:
    return left == right


def q5_trace(matrix: Q5Matrix) -> Q5:
    return sum((matrix[i][i] for i in range(len(matrix))), Q5_ZERO)


def q5_matrix_polynomial_factor(matrix: Q5Matrix, eigenvalue: Q5) -> Q5Matrix:
    return q5_matrix_sub(matrix, q5_matrix_scale(eigenvalue, q5_identity(len(matrix))))


def pair_orbits(group: Sequence[Permutation]) -> list[set[tuple[int, int]]]:
    unassigned = {(i, j) for i in range(PORTS) for j in range(PORTS)}
    result: list[set[tuple[int, int]]] = []
    while unassigned:
        seed = min(unassigned)
        orbit = {(row[seed[0]], row[seed[1]]) for row in group}
        unassigned -= orbit
        result.append(orbit)
    require(sorted(len(orbit) for orbit in result) == [12, 12, 60, 60], "PAIR_ORBITS", "unexpected ordered-pair orbits")
    return result


def commutant_generator(group: Sequence[Permutation]) -> tuple[Q5Matrix, set[tuple[int, int]]]:
    candidates = sorted(
        (orbit for orbit in pair_orbits(group) if len(orbit) == 60),
        key=lambda orbit: tuple(sorted(orbit)),
    )
    require(len(candidates) == 2, "COMMUTANT_ORBIT", "expected two valency-five orbitals")
    chosen = candidates[0]
    matrix = q5_matrix(PORTS, PORTS)
    for i, j in chosen:
        matrix[i][j] = Q5_ONE
    require(matrix == [list(row) for row in zip(*matrix)], "COMMUTANT_SYMMETRY", "orbital matrix is not symmetric")
    require(all(sum((value.a for value in row), Fraction(0)) == 5 for row in matrix), "COMMUTANT_VALENCY", "orbital matrix is not valency five")
    return matrix, chosen


def spectral_projectors(group: Sequence[Permutation]) -> tuple[Q5Matrix, dict[str, Q5Matrix], dict[str, Any]]:
    operator, orbital = commutant_generator(group)
    identity = q5_identity(PORTS)
    eigenvalues = {
        "fixed": Q5(5),
        "three_plus": Q5_SQRT,
        "three_minus": -Q5_SQRT,
        "five": Q5(-1),
    }
    # Exact minimal-polynomial check: (A-5I)(A+I)(A^2-5I)=0.
    square = q5_matrix_mul(operator, operator)
    polynomial = q5_matrix_mul(
        q5_matrix_mul(q5_matrix_sub(operator, q5_matrix_scale(Q5(5), identity)), q5_matrix_add(operator, identity)),
        q5_matrix_sub(square, q5_matrix_scale(Q5(5), identity)),
    )
    require(polynomial == q5_matrix(PORTS, PORTS), "COMMUTANT_POLYNOMIAL", "orbital operator has wrong polynomial")

    projectors: dict[str, Q5Matrix] = {}
    for label, eigenvalue in eigenvalues.items():
        projector = identity
        for other_label, other_eigenvalue in eigenvalues.items():
            if other_label == label:
                continue
            projector = q5_matrix_scale(
                (eigenvalue - other_eigenvalue).inverse(),
                q5_matrix_mul(projector, q5_matrix_polynomial_factor(operator, other_eigenvalue)),
            )
        projectors[label] = projector

    dimensions = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}
    total = q5_matrix(PORTS, PORTS)
    for label, projector in projectors.items():
        require(q5_matrix_mul(projector, projector) == projector, "PROJECTOR_IDEMPOTENT", f"{label} projector is not idempotent")
        require(q5_trace(projector) == Q5(dimensions[label]), "PROJECTOR_TRACE", f"{label} projector has wrong trace")
        require(
            q5_matrix_mul(operator, projector) == q5_matrix_scale(eigenvalues[label], projector),
            "PROJECTOR_EIGENVALUE",
            f"{label} projector has wrong eigenvalue",
        )
        total = q5_matrix_add(total, projector)
    require(total == identity, "PROJECTOR_COMPLETE", "spectral projectors do not sum to identity")
    labels = list(projectors)
    for i, left in enumerate(labels):
        for right in labels[i + 1 :]:
            require(
                q5_matrix_mul(projectors[left], projectors[right]) == q5_matrix(PORTS, PORTS),
                "PROJECTOR_ORTHOGONAL",
                "spectral projectors overlap",
            )
    require(
        all(projectors["fixed"][i][j] == Q5(Fraction(1, 12)) for i in range(PORTS) for j in range(PORTS)),
        "FIXED_PROJECTOR",
        "fixed projector is not the uniform line",
    )
    metadata = {
        "canonical_valency_five_orbital": [list(pair) for pair in sorted(orbital)],
        "operator_sha256": object_sha256(
            [[[value.a.numerator, value.a.denominator, value.b.numerator, value.b.denominator] for value in row] for row in operator]
        ),
        "minimal_polynomial": "(t-5)(t+1)(t^2-5)",
        "sectors": [
            {"id": label, "dimension": dimensions[label], "operator_eigenvalue": encode_q5(eigenvalues[label])}
            for label in labels
        ],
    }
    return operator, projectors, metadata


def derivation_matrices(constants: Sequence[Sequence[Sequence[Sequence[Fraction]]]]) -> list[Q5Matrix]:
    result: list[Q5Matrix] = []
    for parameter in range(PARAMETERS):
        matrix = q5_matrix(PORTS, PORTS)
        for output in range(PORTS):
            for column in range(PORTS):
                matrix[output][column] = Q5(
                    sum((constants[parameter][output][source][column] for source in range(PORTS)), Fraction(0))
                )
        result.append(matrix)
    return result


def eigenvalue_form(projector: Q5Matrix, dimension: int, matrices: Sequence[Q5Matrix]) -> Q5Row:
    values = []
    for matrix in matrices:
        trace = Q5_ZERO
        for i in range(PORTS):
            for j in range(PORTS):
                trace = trace + projector[i][j] * matrix[j][i]
        values.append(trace / dimension)
    return tuple(values)


def channel_coordinate_form(
    vectors: Sequence[Mapping[Coordinate2, Fraction]],
    left_projector: Q5Matrix,
    right_projector: Q5Matrix,
    output_projector: Q5Matrix,
    coordinate: tuple[int, int, int],
) -> Q5Row:
    output, left_input, right_input = coordinate
    form: list[Q5] = []
    for vector in vectors:
        value = Q5_ZERO
        for (source_output, p, q), coefficient in vector.items():
            output_factor = output_projector[output][source_output]
            if not output_factor:
                continue
            input_factor = (
                left_projector[p][left_input] * right_projector[q][right_input]
                - left_projector[q][left_input] * right_projector[p][right_input]
            )
            if input_factor:
                value = value + output_factor * input_factor * coefficient
        form.append(value)
    return tuple(form)


def first_channel_coordinate(
    vectors: Sequence[Mapping[Coordinate2, Fraction]],
    projectors: Mapping[str, Q5Matrix],
    left: str,
    right: str,
    output: str,
) -> tuple[tuple[int, int, int], Q5Row]:
    for coordinate in itertools.product(range(PORTS), repeat=3):
        form = channel_coordinate_form(
            vectors,
            projectors[left],
            projectors[right],
            projectors[output],
            coordinate,
        )
        if any(form):
            return coordinate, form
    raise CertificateError("CHANNEL_ZERO", f"channel {left},{right}->{output} has no coordinate")


def q5_determinant(matrix: Sequence[Sequence[Q5]]) -> Q5:
    work = [list(row) for row in matrix]
    result = Q5_ONE
    for column in range(len(work)):
        pivot = next((r for r in range(column, len(work)) if work[r][column]), None)
        if pivot is None:
            return Q5_ZERO
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        value = work[column][column]
        result = result * value
        for r in range(column + 1, len(work)):
            if work[r][column]:
                factor = work[r][column] / value
                work[r] = [x - factor * y for x, y in zip(work[r], work[column])]
    return result


def q5_inverse(matrix: Sequence[Sequence[Q5]]) -> Q5Matrix:
    size = len(matrix)
    work = [list(row) + q5_identity(size)[i] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((r for r in range(column, size) if work[r][column]), None)
        require(pivot is not None, "CHANNEL_TRANSFORM_SINGULAR", "channel transform is singular")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for r in range(size):
            if r != column and work[r][column]:
                factor = work[r][column]
                work[r] = [x - factor * y for x, y in zip(work[r], work[column])]
    return [row[size:] for row in work]


def q5_rank(rows: Sequence[Sequence[Q5]]) -> int:
    work = [list(row) for row in rows if any(row)]
    if not work:
        return 0
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next((r for r in range(pivot_row, len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for r in range(pivot_row + 1, len(work)):
            if work[r][column]:
                factor = work[r][column]
                work[r] = [x - factor * y for x, y in zip(work[r], work[pivot_row])]
        pivot_row += 1
    return pivot_row


def encode_q5(value: Q5) -> list[int]:
    return [value.a.numerator, value.a.denominator, value.b.numerator, value.b.denominator]


def encode_q5_linear_form(form: Sequence[Q5]) -> list[list[int]]:
    return [[index, *encode_q5(value)] for index, value in enumerate(form) if value]


def encode_q5_equation(row: Sequence[Q5]) -> list[list[int]]:
    return [[a, b, *encode_q5(value)] for (a, b), value in zip(MONOMIALS, row) if value]


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


def channel_transform(
    vectors: Sequence[Mapping[Coordinate2, Fraction]],
    constants: Sequence[Sequence[Sequence[Sequence[Fraction]]]],
    group: Sequence[Permutation],
) -> tuple[list[Q5Row], Q5Matrix, dict[str, Any]]:
    _, projectors, spectral_metadata = spectral_projectors(group)
    derivations = derivation_matrices(constants)
    dimensions = {"fixed": 1, "three_plus": 3, "three_minus": 3, "five": 5}
    eigenforms: dict[str, Q5Row] = {}
    for label in dimensions:
        form = eigenvalue_form(projectors[label], dimensions[label], derivations)
        eigenforms[label] = form
        for parameter, matrix in enumerate(derivations):
            require(
                q5_matrix_mul(matrix, projectors[label])
                == q5_matrix_scale(form[parameter], projectors[label]),
                "DERIVATION_EIGENFORM",
                f"fixed-line adjoint is not scalar on {label}",
            )
    require(not any(eigenforms["fixed"]), "DERIVATION_FIXED_LINE", "fixed-line adjoint does not kill the fixed line")

    transform: list[Q5Row] = [
        eigenforms["three_plus"],
        eigenforms["three_minus"],
        eigenforms["five"],
    ]
    channel_rows = []
    for channel_id, left, right, output, weight in CHANNEL_SPECS:
        coordinate, form = first_channel_coordinate(vectors, projectors, left, right, output)
        transform.append(form)
        channel_rows.append(
            {
                "id": channel_id,
                "domain_sectors": [left, right],
                "output_sector": output,
                "extraction_coordinate": list(coordinate),
                "linear_form_in_x": encode_q5_linear_form(form),
                "derivation_weight_in_d": list(weight),
            }
        )
    determinant = q5_determinant(transform)
    require(bool(determinant), "CHANNEL_TRANSFORM_SINGULAR", "fourteen channel forms are dependent")
    inverse_transform = q5_inverse(transform)
    require(
        q5_matrix_mul([list(row) for row in transform], inverse_transform) == q5_identity(PARAMETERS),
        "CHANNEL_TRANSFORM_INVERSE",
        "channel inverse is not exact",
    )
    metadata = {
        "splitting_field": "Q(sqrt(5)) with sqrt(5)^2=5",
        "spectral_decomposition": spectral_metadata,
        "channel_variables": [
            {"id": "d_plus", "meaning": "fixed-line adjoint eigenvalue on three_plus", "linear_form_in_x": encode_q5_linear_form(transform[0])},
            {"id": "d_minus", "meaning": "fixed-line adjoint eigenvalue on three_minus", "linear_form_in_x": encode_q5_linear_form(transform[1])},
            {"id": "d_five", "meaning": "fixed-line adjoint eigenvalue on five", "linear_form_in_x": encode_q5_linear_form(transform[2])},
            *channel_rows,
        ],
        "transform_determinant": encode_q5(determinant),
        "transform_rows": [encode_q5_linear_form(row) for row in transform],
        "inverse_transform_rows": [encode_q5_linear_form(row) for row in inverse_transform],
    }
    return transform, inverse_transform, metadata


def multiply_linear_forms(left: Sequence[Q5], right: Sequence[Q5]) -> Q5Row:
    result = [Q5_ZERO] * len(MONOMIALS)
    for i, a in enumerate(left):
        if not a:
            continue
        for j, b in enumerate(right):
            if b:
                index = MONOMIAL_INDEX[tuple(sorted((i, j)))]
                result[index] = result[index] + a * b
    return tuple(result)


def linear_combination(forms: Sequence[Sequence[Q5]], coefficients: Sequence[int]) -> Q5Row:
    return tuple(
        sum((Q5(coefficient) * forms[i][column] for i, coefficient in enumerate(coefficients) if coefficient), Q5_ZERO)
        for column in range(PARAMETERS)
    )


def product_equations_in_x(transform: Sequence[Q5Row]) -> list[Q5Row]:
    result = []
    for channel_number, (_, _, _, _, weight) in enumerate(CHANNEL_SPECS):
        weight_form = linear_combination(transform[:3], weight)
        result.append(multiply_linear_forms(weight_form, transform[3 + channel_number]))
    require(q5_rank(result) == 11, "PRODUCT_RANK", "derivation product equations are dependent")
    return result


def transform_quadratic(row: Sequence[Any], inverse_transform: Q5Matrix) -> Q5Row:
    """Substitute x = inverse_transform * y in a quadratic row."""

    result = [Q5_ZERO] * len(MONOMIALS)
    for coefficient_raw, (i, j) in zip(row, MONOMIALS):
        coefficient = Q5.make(coefficient_raw)
        if not coefficient:
            continue
        for a in range(PARAMETERS):
            if not inverse_transform[i][a]:
                continue
            for b in range(PARAMETERS):
                if inverse_transform[j][b]:
                    index = MONOMIAL_INDEX[tuple(sorted((a, b)))]
                    result[index] = result[index] + coefficient * inverse_transform[i][a] * inverse_transform[j][b]
    return tuple(result)


def simple_product_equations() -> list[Q5Row]:
    result: list[Q5Row] = []
    for channel_number, (_, _, _, _, weight) in enumerate(CHANNEL_SPECS):
        row = [Q5_ZERO] * len(MONOMIALS)
        channel_variable = 3 + channel_number
        for d_variable, coefficient in enumerate(weight):
            if coefficient:
                index = MONOMIAL_INDEX[tuple(sorted((d_variable, channel_variable)))]
                row[index] = row[index] + Q5(coefficient)
        result.append(tuple(row))
    return result


def choose_complement(fixed_rows: Sequence[QRow], full_rows: Sequence[QRow]) -> list[QRow]:
    basis, _ = rref(fixed_rows)
    current = list(basis)
    rank = len(current)
    complement: list[QRow] = []
    for row in full_rows:
        candidate_rank = rational_rank(current + [row])
        if candidate_rank > rank:
            current.append(row)
            complement.append(row)
            rank = candidate_rank
    require(rank == 38 and len(complement) == 27, "RESIDUAL_COMPLEMENT", "expected a rank-27 Jacobi complement")
    return complement


def rational_rows_as_q5(rows: Sequence[Sequence[Fraction]]) -> list[Q5Row]:
    return [tuple(Q5(value) for value in row) for row in rows]


def canonical_rref3(rows: Sequence[Sequence[int]]) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    reduced, _ = rref([[Fraction(value) for value in row] for row in rows])
    return tuple(tuple(row) for row in reduced)  # type: ignore[return-value]


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


def rowspace_contains3(space: Sequence[Sequence[Fraction]], row: Sequence[int]) -> bool:
    return len(canonical_rref3([*space, row])) == len(space)


def arrangement_flats() -> list[dict[str, Any]]:
    names = list(DISTINCT_WEIGHTS)
    spaces = {
        canonical_rref3([DISTINCT_WEIGHTS[name] for name in subset])
        for size in range(len(names) + 1)
        for subset in itertools.combinations(names, size)
    }
    result = []
    for number, space in enumerate(sorted(spaces, key=lambda rows: (len(rows), rows))):
        vanishing_weights = [name for name, row in DISTINCT_WEIGHTS.items() if rowspace_contains3(space, row)]
        allowed_channels = [
            spec[0]
            for spec in CHANNEL_SPECS
            if rowspace_contains3(space, spec[4])
        ]
        result.append(
            {
                "flat_id": f"F{number:02d}",
                "dimension": 3 - len(space),
                "defining_rowspace_rref": [
                    [[value.numerator, value.denominator] for value in row]
                    for row in space
                ],
                "weights_vanishing_identically": vanishing_weights,
                "channels_not_forced_zero_by_fixed_line_sector": allowed_channels,
                "channels_forced_zero_on_relative_open_flat": [
                    spec[0] for spec in CHANNEL_SPECS if spec[0] not in allowed_channels
                ],
                "scope": "relative-open part of this flat; larger subflats are separate rows",
            }
        )
    counts: dict[int, int] = {}
    for row in result:
        counts[row["dimension"]] = counts.get(row["dimension"], 0) + 1
    require(len(result) == 28 and counts == {3: 1, 2: 8, 1: 18, 0: 1}, "ARRANGEMENT_FLATS", "unexpected derivation-weight arrangement")
    return result


def validate_import_contract(filename: str, roots: set[str]) -> None:
    expected = EXACT_IMPORT_ROOTS_BY_FILE[filename]
    forbidden = RUNTIME_FORBIDDEN_IMPORT_ROOTS.get(filename, set())
    require(not (roots & forbidden), "FIREWALL_RUNTIME_IMPORT", f"{filename} imports forbidden runtime roots")
    require(
        roots == expected,
        "FIREWALL_IMPORT_CONTRACT",
        f"{filename} import roots differ from its exact per-file contract",
    )


def audit_imports(paths: Sequence[Path]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for path in paths:
        roots: set[str] = set()
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                roots.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                roots.add(node.module.split(".")[0])
        validate_import_contract(path.name, roots)
        result[path.name] = sorted(roots)
    require(set(result) == set(EXACT_IMPORT_ROOTS_BY_FILE), "FIREWALL_IMPORT_FILES", "per-file import audit is incomplete")
    return dict(sorted(result.items()))


def expect_error(name: str, expected: str, action: Any) -> dict[str, Any]:
    actual = "NO_ERROR"
    try:
        action()
    except CertificateError as exc:
        actual = exc.code
    return {"name": name, "expected_error": expected, "actual_error": actual, "passed": actual == expected}


def validate_integer_entries(rows: Any, width: int, code: str) -> None:
    require(isinstance(rows, list), code, "equation collection is not a list")
    for row in rows:
        require(isinstance(row, list), code, "equation row is not a list")
        for entry in row:
            require(isinstance(entry, list) and len(entry) == width and all(type(value) is int for value in entry), code, "bad encoded equation entry")


def validate_stage2_keysets(system: Any, receipt: Any) -> None:
    require_keys(system, SYSTEM_KEYS, "SYSTEM_KEYS")
    require_keys(system.get("upstream"), UPSTREAM_KEYS, "SYSTEM_UPSTREAM_KEYS")
    jacobi = system.get("jacobi_tensor")
    require_keys(jacobi, JACOBI_KEYS, "SYSTEM_JACOBI_KEYS")
    character = jacobi.get("target_invariant_character_certificate")
    require_keys(character, TARGET_CHARACTER_KEYS, "SYSTEM_TARGET_CHARACTER_KEYS")
    require(isinstance(character.get("class_rows"), list), "SYSTEM_TARGET_CLASS_ROWS", "target character rows are absent")
    for row in character["class_rows"]:
        require_keys(row, TARGET_CLASS_ROW_KEYS, "SYSTEM_TARGET_CLASS_ROW_KEYS")
    representatives = jacobi.get("representative_equations")
    require(isinstance(representatives, list), "SYSTEM_REPRESENTATIVE_ROWS", "representatives are absent")
    for row in representatives:
        require_keys(row, REPRESENTATIVE_KEYS, "SYSTEM_REPRESENTATIVE_KEYS")
        validate_integer_entries([row.get("equation")], 4, "SYSTEM_RATIONAL_EQUATION_FORMAT")
    validate_integer_entries(jacobi.get("reduced_integer_equations"), 3, "SYSTEM_INTEGER_EQUATION_FORMAT")

    fixed = system.get("fixed_line_reduction")
    require_keys(fixed, FIXED_REDUCTION_KEYS, "SYSTEM_FIXED_KEYS")
    validate_integer_entries(fixed.get("contracted_integer_equations_in_x"), 3, "SYSTEM_INTEGER_EQUATION_FORMAT")
    validate_integer_entries(fixed.get("residual_integer_equations_in_x"), 3, "SYSTEM_INTEGER_EQUATION_FORMAT")
    validate_integer_entries(fixed.get("residual_equations_in_channel_variables"), 6, "SYSTEM_Q5_EQUATION_FORMAT")
    channel = fixed.get("channel_decomposition")
    require_keys(channel, CHANNEL_DECOMPOSITION_KEYS, "SYSTEM_CHANNEL_KEYS")
    variables = channel.get("channel_variables")
    require(isinstance(variables, list) and len(variables) == PARAMETERS, "SYSTEM_CHANNEL_VARIABLES", "channel-variable list differs")
    for number, row in enumerate(variables):
        require_keys(row, D_CHANNEL_KEYS if number < 3 else T_CHANNEL_KEYS, "SYSTEM_CHANNEL_VARIABLE_KEYS")
        validate_integer_entries([row.get("linear_form_in_x")], 5, "SYSTEM_Q5_LINEAR_FORMAT")
    validate_integer_entries(channel.get("transform_rows"), 5, "SYSTEM_Q5_LINEAR_FORMAT")
    validate_integer_entries(channel.get("inverse_transform_rows"), 5, "SYSTEM_Q5_LINEAR_FORMAT")
    spectral = channel.get("spectral_decomposition")
    require_keys(spectral, SPECTRAL_KEYS, "SYSTEM_SPECTRAL_KEYS")
    sectors = spectral.get("sectors")
    require(isinstance(sectors, list), "SYSTEM_SECTOR_ROWS", "spectral sectors are absent")
    for row in sectors:
        require_keys(row, SECTOR_KEYS, "SYSTEM_SECTOR_KEYS")
    products = fixed.get("derivation_product_equations_in_channel_variables")
    require(isinstance(products, list), "SYSTEM_PRODUCT_ROWS", "derivation products are absent")
    for row in products:
        require_keys(row, PRODUCT_DECLARATION_KEYS, "SYSTEM_PRODUCT_KEYS")
        validate_integer_entries([row.get("equation")], 6, "SYSTEM_Q5_EQUATION_FORMAT")
    require_keys(fixed.get("derivation_only_subspace"), DERIVATION_ONLY_KEYS, "SYSTEM_DERIVATION_KEYS")
    require_keys(fixed.get("central_specialization"), CENTRAL_SPECIALIZATION_KEYS, "SYSTEM_CENTRAL_KEYS")
    arrangement = fixed.get("derivation_weight_arrangement")
    require_keys(arrangement, ARRANGEMENT_KEYS, "SYSTEM_ARRANGEMENT_KEYS")
    require_keys(arrangement.get("flat_dimension_counts"), {"0", "1", "2", "3"}, "SYSTEM_FLAT_COUNT_KEYS")
    flats = arrangement.get("relative_open_flats")
    require(isinstance(flats, list), "SYSTEM_FLAT_ROWS", "arrangement flats are absent")
    for row in flats:
        require_keys(row, FLAT_KEYS, "SYSTEM_FLAT_KEYS")

    require_keys(receipt, RECEIPT_KEYS, "RECEIPT_KEYS")
    require_keys(receipt.get("semantic_inputs"), SEMANTIC_INPUT_KEYS, "RECEIPT_SEMANTIC_KEYS")
    require_keys(receipt["semantic_inputs"].get("paths"), UPSTREAM_KEYS, "RECEIPT_SEMANTIC_PATH_KEYS")
    firewall = receipt.get("target_firewall")
    require_keys(firewall, FIREWALL_KEYS, "RECEIPT_FIREWALL_KEYS")
    for field in ("audited_import_roots_by_file", "exact_import_allowlists_by_file"):
        require_keys(firewall.get(field), set(EXACT_IMPORT_ROOTS_BY_FILE), "RECEIPT_IMPORT_FILE_KEYS")
    require_keys(firewall.get("runtime_forbidden_import_roots_by_file"), set(RUNTIME_FORBIDDEN_IMPORT_ROOTS), "RECEIPT_FORBIDDEN_FILE_KEYS")
    require_keys(receipt.get("attained"), ATTAINED_KEYS, "RECEIPT_ATTAINED_KEYS")
    require_keys(receipt.get("not_attained"), NOT_ATTAINED_KEYS, "RECEIPT_NOT_ATTAINED_KEYS")
    require_keys(receipt.get("system_artifact"), SYSTEM_ARTIFACT_KEYS, "RECEIPT_SYSTEM_ARTIFACT_KEYS")
    require_keys(receipt.get("implementation_pins"), IMPLEMENTATION_PIN_KEYS, "RECEIPT_PIN_KEYS")
    mutations = receipt.get("mutation_tests")
    require(isinstance(mutations, list), "RECEIPT_MUTATION_ROWS", "mutation rows are absent")
    for row in mutations:
        require_keys(row, MUTATION_ROW_KEYS, "RECEIPT_MUTATION_KEYS")


def decode_q5_value(raw: Sequence[int]) -> Q5:
    require(isinstance(raw, list) and len(raw) == 4 and all(type(value) is int for value in raw), "SYSTEM_Q5_VALUE", "bad Q(sqrt(5)) value")
    require(raw[1] > 0 and raw[3] > 0, "SYSTEM_Q5_VALUE", "non-positive Q(sqrt(5)) denominator")
    return Q5(Fraction(raw[0], raw[1]), Fraction(raw[2], raw[3]))


def decode_q5_linear(entries: Sequence[Sequence[int]]) -> Q5Row:
    result = [Q5_ZERO] * PARAMETERS
    seen: set[int] = set()
    for index, an, ad, bn, bd in entries:
        require(0 <= index < PARAMETERS and index not in seen and ad > 0 and bd > 0, "SYSTEM_Q5_LINEAR", "bad Q(sqrt(5)) linear entry")
        seen.add(index)
        result[index] = Q5(Fraction(an, ad), Fraction(bn, bd))
    return tuple(result)


def validate_system_declarations(system: Mapping[str, Any]) -> None:
    validate_stage2_keysets(system, {
        "schema": "oph.a5_jacobi_stage2.receipt.v1",
        "issue": 566,
        "status": STATUS,
        "claim_boundary": CLAIM_BOUNDARY,
        "semantic_inputs": {"count": 2, "paths": system["upstream"], "stage1_later_gates_all_false": True},
        "target_firewall": {
            "enabled": True,
            "accepted_input_schemas": ["oph.a5_alternating_bracket_reynolds_basis.v1", "oph.a5_alternating_bracket_space_stage1.receipt.v1"],
            "exact_input_keysets_enforced": True,
            "arbitrary_input_path_allowed": False,
            "coefficient_target_input_allowed": False,
            "measurement_input_allowed": False,
            "environment_target_read": False,
            "network_target_read": False,
            "nonstdlib_imports": [],
            "audited_import_roots_by_file": {name: sorted(roots) for name, roots in EXACT_IMPORT_ROOTS_BY_FILE.items()},
            "exact_import_allowlists_by_file": {name: sorted(roots) for name, roots in EXACT_IMPORT_ROOTS_BY_FILE.items()},
            "runtime_forbidden_import_roots_by_file": {name: sorted(roots) for name, roots in RUNTIME_FORBIDDEN_IMPORT_ROOTS.items()},
        },
        "attained": {key: False for key in ATTAINED_KEYS},
        "not_attained": {key: False for key in NOT_ATTAINED_KEYS},
        "system_artifact": {"path": "", "canonical_json_sha256": ""},
        "mutation_tests": [],
        "implementation_pins": {"producer_sha256": "", "independent_verifier_sha256": "", "test_sha256": ""},
        "receipt_sha256": "",
    })
    require(system.get("schema") == "oph.a5_jacobi_system_reduction.v1", "SYSTEM_SCHEMA", "system schema differs")
    require(system.get("issue") == 566, "SYSTEM_ISSUE", "system issue differs")
    require(system.get("field") == "Q", "SYSTEM_FIELD", "system field differs")
    require(system.get("parameter_variables") == [f"x{i:02d}" for i in range(PARAMETERS)], "SYSTEM_VARIABLES", "parameter variables differ")
    require(system.get("monomial_order") == [list(row) for row in MONOMIALS], "SYSTEM_MONOMIALS", "monomial order differs")
    upstream = system["upstream"]
    require(
        upstream.get("basis_path") == "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json"
        and upstream.get("receipt_path") == "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json",
        "SYSTEM_UPSTREAM_PATH",
        "upstream path differs",
    )
    current_basis = json.loads(STAGE1_BASIS.read_text(encoding="utf-8"))
    current_receipt = json.loads(STAGE1_RECEIPT.read_text(encoding="utf-8"))
    require(
        upstream.get("basis_file_sha256") == file_sha256(STAGE1_BASIS)
        and upstream.get("basis_canonical_json_sha256") == object_sha256(current_basis)
        and upstream.get("receipt_file_sha256") == file_sha256(STAGE1_RECEIPT)
        and upstream.get("receipt_self_sha256") == current_receipt.get("receipt_sha256")
        and upstream.get("proper_group_sha256") == current_receipt["proper_port_action"].get("permutation_rows_sha256"),
        "SYSTEM_UPSTREAM_PINS",
        "upstream pin differs",
    )
    jacobi = system["jacobi_tensor"]
    require(
        jacobi.get("convention") == "J(u,v,w)=B(B(u,v),w)+B(B(v,w),u)+B(B(w,u),v)"
        and jacobi.get("coordinate_equation_count") == 2640
        and jacobi.get("all_coordinate_equations_nonzero") is True
        and jacobi.get("target_space_dimension") == 2640
        and jacobi.get("signed_target_orbit_count") == 44
        and jacobi.get("signed_target_orbit_sizes") == [60] * 44
        and jacobi.get("representative_coefficient_rank_over_Q") == 38
        and jacobi.get("reduced_integer_equation_count") == 38
        and len(jacobi.get("representative_equations", [])) == 44
        and len(jacobi.get("reduced_integer_equations", [])) == 38
        and jacobi.get("zero_locus_equivalence") == "all 2640 components vanish iff the 44 orbit representatives vanish iff the coefficient-row rank 38 reduced integer quadrics vanish",
        "SYSTEM_JACOBI_CLAIMS",
        "Jacobi count, rank, or scope declaration differs",
    )
    fixed = system["fixed_line_reduction"]
    require(
        fixed.get("fixed_vector") == [1] * PORTS
        and fixed.get("contracted_coordinate_count") == 792
        and fixed.get("contracted_coefficient_rank_over_Q") == 11
        and len(fixed.get("contracted_integer_equations_in_x", [])) == 11
        and fixed.get("product_equation_rank_over_Qsqrt5") == 11
        and fixed.get("product_rowspace_equals_contracted_jacobi_rowspace") is True
        and len(fixed.get("residual_integer_equations_in_x", [])) == 27
        and len(fixed.get("residual_equations_in_channel_variables", [])) == 27
        and fixed.get("residual_rank_over_Qsqrt5") == 27,
        "SYSTEM_RANKS",
        "coefficient-row rank declaration differs",
    )
    require(fixed.get("full_rank_split") == [11, 27, 38], "SYSTEM_SPLIT", "rowspace split declaration differs")
    channel = fixed["channel_decomposition"]
    require(channel.get("splitting_field") == "Q(sqrt(5)) with sqrt(5)^2=5", "SYSTEM_SPLITTING_FIELD", "splitting field differs")
    spectral = channel["spectral_decomposition"]
    require(spectral.get("minimal_polynomial") == "(t-5)(t+1)(t^2-5)", "SYSTEM_MINIMAL_POLYNOMIAL", "minimal polynomial differs")
    require(
        [(row["id"], row["dimension"], row["operator_eigenvalue"]) for row in spectral["sectors"]]
        == [
            ("fixed", 1, [5, 1, 0, 1]),
            ("three_plus", 3, [0, 1, 1, 1]),
            ("three_minus", 3, [0, 1, -1, 1]),
            ("five", 5, [-1, 1, 0, 1]),
        ],
        "SYSTEM_SECTOR_DIMENSIONS",
        "spectral sector dimensions or eigenvalues differ",
    )
    transform = [decode_q5_linear(row) for row in channel["transform_rows"]]
    inverse_transform = [list(decode_q5_linear(row)) for row in channel["inverse_transform_rows"]]
    require(q5_matrix_mul([list(row) for row in transform], inverse_transform) == q5_identity(PARAMETERS), "SYSTEM_TRANSFORM", "channel transform inverse differs")
    require(decode_q5_value(channel["transform_determinant"]) == q5_determinant(transform), "SYSTEM_DETERMINANT", "transform determinant declaration differs")
    require(bool(q5_determinant(transform)), "SYSTEM_DETERMINANT", "transform determinant is zero")
    variables = channel["channel_variables"]
    require([row["id"] for row in variables[:3]] == ["d_plus", "d_minus", "d_five"], "SYSTEM_CHANNEL_IDS", "derivation channel ids differ")
    for row, spec in zip(variables[3:], CHANNEL_SPECS):
        require(
            row["id"] == spec[0]
            and row["domain_sectors"] == [spec[1], spec[2]]
            and row["output_sector"] == spec[3]
            and row["derivation_weight_in_d"] == list(spec[4]),
            "SYSTEM_CHANNEL_METADATA",
            "channel declaration differs",
        )
    products = fixed["derivation_product_equations_in_channel_variables"]
    require(
        [(row["channel_id"], row["weight_coefficients_on_d_plus_d_minus_d_five"]) for row in products]
        == [(spec[0], list(spec[4])) for spec in CHANNEL_SPECS],
        "SYSTEM_DERIVATION_PRODUCTS",
        "derivation product declaration differs",
    )
    require(
        fixed["derivation_only_subspace"]
        == {
            "dimension": 3,
            "definition": "all eleven non-derivation channel variables are zero",
            "all_jacobi_equations_vanish_identically": True,
            "scope": DERIVATION_ONLY_SCOPE,
        },
        "SYSTEM_DERIVATION_SCOPE",
        "derivation-only scope differs",
    )
    require(
        fixed["central_specialization"]
        == {
            "definition": "d_plus=d_minus=d_five=0",
            "remaining_channel_variable_count": 11,
            "residual_coefficient_rank_over_Qsqrt5": 27,
            "solution_set_classified": False,
        },
        "SYSTEM_CENTRAL_SPECIALIZATION",
        "central specialization differs",
    )
    arrangement = fixed["derivation_weight_arrangement"]
    require(
        arrangement.get("distinct_weight_hyperplane_count") == 8
        and arrangement.get("flat_count") == 28
        and arrangement.get("flat_dimension_counts") == {"0": 1, "1": 18, "2": 8, "3": 1}
        and len(arrangement.get("relative_open_flats", [])) == 28
        and arrangement.get("scope") == ARRANGEMENT_SCOPE,
        "SYSTEM_ARRANGEMENT_CLAIMS",
        "arrangement count or scope differs",
    )


def validate_receipt_declarations(system: Mapping[str, Any], receipt: Mapping[str, Any], imports: Mapping[str, Sequence[str]]) -> None:
    validate_stage2_keysets(system, receipt)
    unhashed = dict(receipt)
    stored = unhashed.pop("receipt_sha256", None)
    require(stored == object_sha256(unhashed), "RECEIPT_HASH", "receipt self-hash differs")
    require(receipt.get("schema") == "oph.a5_jacobi_stage2.receipt.v1", "RECEIPT_SCHEMA", "receipt schema differs")
    require(receipt.get("issue") == 566, "RECEIPT_ISSUE", "receipt issue differs")
    require(receipt.get("status") == STATUS, "RECEIPT_STATUS", "receipt status differs")
    require(receipt.get("claim_boundary") == CLAIM_BOUNDARY, "RECEIPT_SCOPE", "receipt claim boundary differs")
    require(
        receipt.get("semantic_inputs")
        == {"count": 2, "paths": system["upstream"], "stage1_later_gates_all_false": True},
        "RECEIPT_SEMANTIC_INPUTS",
        "semantic-input declaration differs",
    )
    expected_firewall = {
        "enabled": True,
        "accepted_input_schemas": [
            "oph.a5_alternating_bracket_reynolds_basis.v1",
            "oph.a5_alternating_bracket_space_stage1.receipt.v1",
        ],
        "exact_input_keysets_enforced": True,
        "arbitrary_input_path_allowed": False,
        "coefficient_target_input_allowed": False,
        "measurement_input_allowed": False,
        "environment_target_read": False,
        "network_target_read": False,
        "nonstdlib_imports": [],
        "audited_import_roots_by_file": dict(imports),
        "exact_import_allowlists_by_file": {
            filename: sorted(roots) for filename, roots in sorted(EXACT_IMPORT_ROOTS_BY_FILE.items())
        },
        "runtime_forbidden_import_roots_by_file": {
            filename: sorted(roots) for filename, roots in sorted(RUNTIME_FORBIDDEN_IMPORT_ROOTS.items())
        },
    }
    require(receipt.get("target_firewall") == expected_firewall, "RECEIPT_FIREWALL", "firewall declaration differs")
    expected_attained = {
        "general_parameter_count": 14,
        "raw_jacobi_coordinate_count": 2640,
        "equivariant_target_orbit_count": 44,
        "independent_jacobi_quadrics": 38,
        "fixed_line_product_quadrics": 11,
        "residual_independent_quadrics": 27,
        "exact_channel_coordinate_change": True,
        "derivation_weight_relative_open_flat_count": 28,
        "generic_weight_stratum_forces_all_non_derivation_channels_zero": True,
        "derivation_only_three_parameter_family_certified": True,
    }
    require(receipt.get("attained") == expected_attained, "RECEIPT_COUNTS", "attained count or rank declaration differs")
    require(receipt.get("not_attained") == {key: False for key in NOT_ATTAINED_KEYS}, "RECEIPT_BOUNDARY", "non-attained boundary differs")
    require(
        receipt.get("system_artifact")
        == {
            "path": "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_system_reduction.json",
            "canonical_json_sha256": object_sha256(system),
        },
        "RECEIPT_SYSTEM_ARTIFACT",
        "system artifact path or hash differs",
    )
    pins = receipt["implementation_pins"]
    require(
        pins.get("producer_sha256") == file_sha256(Path(__file__).resolve())
        and pins.get("independent_verifier_sha256") == file_sha256(VERIFIER_PATH.resolve())
        and pins.get("test_sha256") == file_sha256(TEST_PATH.resolve()),
        "RECEIPT_IMPLEMENTATION_PIN",
        "implementation pin differs",
    )
    require(all(row.get("passed") is True for row in receipt["mutation_tests"]), "RECEIPT_MUTATION", "producer mutation failed")


def validate_target_covariance(
    raw: Mapping[Coordinate3, QRow], orbits: Sequence[Mapping[Coordinate3, int]]
) -> None:
    for orbit in orbits:
        seed = min(orbit)
        seed_row = raw[seed]
        for coordinate, sign in orbit.items():
            require(raw[coordinate] == scalar_row(seed_row, Fraction(sign)), "JACOBI_COVARIANCE", "raw Jacobi tensor is not equivariant")


def validate_product_rowspace(fixed_basis: Sequence[Sequence[Q5]], products: Sequence[Sequence[Q5]]) -> None:
    require(q5_rank(fixed_basis) == 11 and q5_rank(products) == 11, "PRODUCT_RANK", "fixed-line or product rank differs")
    require(q5_rank([*fixed_basis, *products]) == 11, "PRODUCT_ROWSPACE", "product equations do not equal the fixed-line rowspace")


def validate_full_decomposition(
    full_rows: Sequence[Sequence[Q5]], products: Sequence[Sequence[Q5]], residual: Sequence[Sequence[Q5]]
) -> None:
    require(len(residual) == 27, "RESIDUAL_COUNT", "residual complement must contain twenty-seven rows")
    require(q5_rank(full_rows) == 38, "FULL_RANK", "full transformed system rank differs")
    require(q5_rank([*products, *residual]) == 38, "RESIDUAL_SPAN", "product plus residual rows do not span the full system")
    require(q5_rank([*full_rows, *products, *residual]) == 38, "RESIDUAL_SPAN", "decomposed system changes the full rowspace")


def run_mutation_tests(
    basis_raw: Mapping[str, Any],
    receipt_raw: Mapping[str, Any],
    raw: Mapping[Coordinate3, QRow],
    orbits: Sequence[Mapping[Coordinate3, int]],
    transform: Sequence[Q5Row],
    fixed_y: Sequence[Q5Row],
    products_y: Sequence[Q5Row],
    full_y: Sequence[Q5Row],
    residual_y: Sequence[Q5Row],
) -> list[dict[str, Any]]:
    tests: list[dict[str, Any]] = []

    injected_basis = copy.deepcopy(basis_raw)
    injected_basis["target_payload"] = {"accepted": True}
    tests.append(
        expect_error(
            "inject_unregistered_upstream_target",
            "FIREWALL_BASIS_KEYS",
            lambda: validate_upstream_keysets(injected_basis, receipt_raw),
        )
    )

    tampered_raw = dict(raw)
    orbit = orbits[0]
    coordinate = next(item for item in sorted(orbit) if item != min(orbit))
    row = list(tampered_raw[coordinate])
    row[0] += 1
    tampered_raw[coordinate] = tuple(row)
    tests.append(
        expect_error(
            "tamper_raw_jacobi_coordinate",
            "JACOBI_COVARIANCE",
            lambda: validate_target_covariance(tampered_raw, orbits),
        )
    )

    singular = [tuple(row) for row in transform]
    singular[-1] = singular[-2]
    tests.append(
        expect_error(
            "duplicate_channel_coordinate",
            "CHANNEL_TRANSFORM_SINGULAR",
            lambda: require(bool(q5_determinant(singular)), "CHANNEL_TRANSFORM_SINGULAR", "mutated transform is singular"),
        )
    )

    bad_products = [tuple(row) for row in products_y]
    changed = list(bad_products[0])
    changed[MONOMIAL_INDEX[(3, 3)]] = changed[MONOMIAL_INDEX[(3, 3)]] + Q5_ONE
    bad_products[0] = tuple(changed)
    tests.append(
        expect_error(
            "alter_derivation_product",
            "PRODUCT_ROWSPACE",
            lambda: validate_product_rowspace(fixed_y, bad_products),
        )
    )

    tests.append(
        expect_error(
            "drop_residual_equation",
            "RESIDUAL_COUNT",
            lambda: validate_full_decomposition(full_y, products_y, residual_y[:-1]),
        )
    )

    duplicated_residual = list(residual_y)
    duplicated_residual[-1] = duplicated_residual[-2]
    tests.append(
        expect_error(
            "duplicate_residual_equation",
            "RESIDUAL_SPAN",
            lambda: validate_full_decomposition(full_y, products_y, duplicated_residual),
        )
    )

    promoted_verifier_imports = set(EXACT_IMPORT_ROOTS_BY_FILE["verify.py"])
    promoted_verifier_imports.add("subprocess")
    tests.append(
        expect_error(
            "promote_verifier_runtime_import_boundary",
            "FIREWALL_RUNTIME_IMPORT",
            lambda: validate_import_contract("verify.py", promoted_verifier_imports),
        )
    )

    require(all(test["passed"] for test in tests), "MUTATION_TEST", "a stage-two mutation control failed")
    return tests


def rehash_receipt(receipt: dict[str, Any]) -> None:
    unhashed = dict(receipt)
    unhashed.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = object_sha256(unhashed)


def run_artifact_mutation_tests(
    system: Mapping[str, Any], receipt: Mapping[str, Any], imports: Mapping[str, Sequence[str]]
) -> list[dict[str, Any]]:
    tests: list[dict[str, Any]] = []

    def receipt_attack(name: str, expected: str, mutate: Any) -> None:
        changed = copy.deepcopy(receipt)
        mutate(changed)
        rehash_receipt(changed)
        tests.append(
            expect_error(
                name,
                expected,
                lambda: (
                    validate_system_declarations(system),
                    validate_receipt_declarations(system, changed, imports),
                ),
            )
        )

    def system_attack(name: str, expected: str, mutate: Any) -> None:
        changed_system = copy.deepcopy(system)
        changed_receipt = copy.deepcopy(receipt)
        mutate(changed_system)
        changed_receipt["system_artifact"]["canonical_json_sha256"] = object_sha256(changed_system)
        changed_receipt["semantic_inputs"]["paths"] = copy.deepcopy(changed_system["upstream"])
        rehash_receipt(changed_receipt)
        tests.append(
            expect_error(
                name,
                expected,
                lambda: (
                    validate_system_declarations(changed_system),
                    validate_receipt_declarations(changed_system, changed_receipt, imports),
                ),
            )
        )

    receipt_attack(
        "rehashed_false_full_classification_claim",
        "RECEIPT_BOUNDARY",
        lambda value: value["not_attained"].__setitem__("full_jacobi_solution_variety_classified", True),
    )
    receipt_attack(
        "rehashed_false_preferred_algebra_claim",
        "RECEIPT_ATTAINED_KEYS",
        lambda value: value["attained"].__setitem__("standard_model_identified", True),
    )
    receipt_attack(
        "rehashed_false_source_selection_claim",
        "RECEIPT_BOUNDARY",
        lambda value: value["not_attained"].__setitem__("source_selection", True),
    )
    receipt_attack(
        "rehashed_false_compactness_claim",
        "RECEIPT_BOUNDARY",
        lambda value: value["not_attained"].__setitem__("compactness_established", True),
    )
    receipt_attack("rehashed_receipt_issue_change", "RECEIPT_ISSUE", lambda value: value.__setitem__("issue", 567))
    receipt_attack(
        "rehashed_receipt_count_change",
        "RECEIPT_COUNTS",
        lambda value: value["attained"].__setitem__("raw_jacobi_coordinate_count", 2639),
    )
    receipt_attack(
        "rehashed_semantic_path_change",
        "RECEIPT_SEMANTIC_INPUTS",
        lambda value: value["semantic_inputs"]["paths"].__setitem__("basis_path", "unregistered.json"),
    )
    receipt_attack(
        "rehashed_firewall_change",
        "RECEIPT_FIREWALL",
        lambda value: value["target_firewall"].__setitem__("measurement_input_allowed", True),
    )
    receipt_attack("rehashed_receipt_extra_key", "RECEIPT_KEYS", lambda value: value.__setitem__("extra_claim", True))

    system_attack("rehashed_system_issue_change", "SYSTEM_ISSUE", lambda value: value.__setitem__("issue", 567))
    system_attack("rehashed_system_count_change", "SYSTEM_JACOBI_CLAIMS", lambda value: value["jacobi_tensor"].__setitem__("coordinate_equation_count", 2639))
    system_attack("rehashed_system_rank_change", "SYSTEM_RANKS", lambda value: value["fixed_line_reduction"].__setitem__("residual_rank_over_Qsqrt5", 26))
    system_attack("rehashed_system_split_change", "SYSTEM_SPLIT", lambda value: value["fixed_line_reduction"].__setitem__("full_rank_split", [11, 26, 37]))
    system_attack("rehashed_system_field_change", "SYSTEM_FIELD", lambda value: value.__setitem__("field", "R"))
    system_attack(
        "rehashed_minimal_polynomial_change",
        "SYSTEM_MINIMAL_POLYNOMIAL",
        lambda value: value["fixed_line_reduction"]["channel_decomposition"]["spectral_decomposition"].__setitem__("minimal_polynomial", "(t-5)(t+1)"),
    )
    system_attack(
        "rehashed_transform_determinant_change",
        "SYSTEM_DETERMINANT",
        lambda value: value["fixed_line_reduction"]["channel_decomposition"].__setitem__("transform_determinant", [1, 1, 0, 1]),
    )
    system_attack(
        "rehashed_sector_dimension_change",
        "SYSTEM_SECTOR_DIMENSIONS",
        lambda value: value["fixed_line_reduction"]["channel_decomposition"]["spectral_decomposition"]["sectors"][1].__setitem__("dimension", 4),
    )
    system_attack("rehashed_system_extra_key", "SYSTEM_KEYS", lambda value: value.__setitem__("extra_claim", True))
    require(all(row["passed"] for row in tests), "ARTIFACT_MUTATION_TEST", "a rehashed whole-artifact mutation survived")
    return tests


def build_artifacts() -> tuple[dict[str, Any], dict[str, Any]]:
    basis_raw, stage1_receipt, vectors, group = load_stage1()
    constants = full_structure_constants(vectors)
    raw = build_raw_jacobi(constants)
    orbits = target_orbits(group)
    validate_target_covariance(raw, orbits)
    target_character = character_target(group)

    representatives = []
    representative_rows: list[QRow] = []
    for orbit in orbits:
        seed = min(orbit)
        row = raw[seed]
        representative_rows.append(row)
        representatives.append(
            {
                "target_coordinate": list(seed),
                "signed_orbit_size": len(orbit),
                "equation": encode_rational_equation(row),
            }
        )
    full_rref, full_pivots = rref(representative_rows)
    require(len(full_rref) == 38, "JACOBI_RANK", "Jacobi representative rank is not thirty-eight")
    reduced_integer = [primitive_integer_row(row) for row in full_rref]
    require(rational_rank([[Fraction(value) for value in row] for row in reduced_integer]) == 38, "JACOBI_REDUCED_RANK", "integer reduction lost rank")

    fixed_components = invariant_line_rows(raw)
    fixed_rref, fixed_pivots = rref(fixed_components)
    require(len(fixed_rref) == 11, "FIXED_LINE_RANK", "fixed-line rref rank differs")
    fixed_integer = [primitive_integer_row(row) for row in fixed_rref]
    residual_x = choose_complement(fixed_rref, full_rref)
    residual_integer = [primitive_integer_row(row) for row in residual_x]
    require(rational_rank([*fixed_rref, *residual_x]) == 38, "RESIDUAL_SPAN", "rational residual complement is incomplete")

    transform, inverse_transform, channel_metadata = channel_transform(vectors, constants, group)
    product_x = product_equations_in_x(transform)
    fixed_q5 = rational_rows_as_q5(fixed_rref)
    validate_product_rowspace(fixed_q5, product_x)

    full_y = [transform_quadratic(row, inverse_transform) for row in full_rref]
    fixed_y = [transform_quadratic(row, inverse_transform) for row in fixed_rref]
    residual_y = [
        transform_quadratic([Fraction(value) for value in row], inverse_transform)
        for row in residual_integer
    ]
    products_y = simple_product_equations()
    validate_product_rowspace(fixed_y, products_y)
    validate_full_decomposition(full_y, products_y, residual_y)

    pure_derivation_indices = [
        MONOMIAL_INDEX[(a, b)] for a in range(3) for b in range(a, 3)
    ]
    require(
        all(not row[index] for row in residual_y for index in pure_derivation_indices),
        "DERIVATION_ONLY_FAMILY",
        "residual equations do not vanish on the derivation-only subspace",
    )
    central_columns = [
        MONOMIAL_INDEX[(a, b)] for a in range(3, PARAMETERS) for b in range(a, PARAMETERS)
    ]
    central_rows = [tuple(row[index] for index in central_columns) for row in residual_y]
    central_rank = q5_rank(central_rows)
    flats = arrangement_flats()

    mutation_tests = run_mutation_tests(
        basis_raw,
        stage1_receipt,
        raw,
        orbits,
        transform,
        fixed_y,
        products_y,
        full_y,
        residual_y,
    )
    imports = audit_imports((Path(__file__).resolve(), VERIFIER_PATH.resolve(), TEST_PATH.resolve()))

    system = {
        "schema": "oph.a5_jacobi_system_reduction.v1",
        "issue": 566,
        "field": "Q",
        "parameter_variables": [f"x{i:02d}" for i in range(PARAMETERS)],
        "monomial_order": [list(monomial) for monomial in MONOMIALS],
        "upstream": {
            "basis_path": "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_reynolds_basis.json",
            "basis_file_sha256": file_sha256(STAGE1_BASIS),
            "basis_canonical_json_sha256": object_sha256(basis_raw),
            "receipt_path": "code/a5_closure/issue_566_bracket_space_stage1/a5_alternating_bracket_space_stage1.receipt.json",
            "receipt_file_sha256": file_sha256(STAGE1_RECEIPT),
            "receipt_self_sha256": stage1_receipt["receipt_sha256"],
            "proper_group_sha256": stage1_receipt["proper_port_action"]["permutation_rows_sha256"],
        },
        "jacobi_tensor": {
            "convention": "J(u,v,w)=B(B(u,v),w)+B(B(v,w),u)+B(B(w,u),v)",
            "coordinate_equation_count": len(raw),
            "all_coordinate_equations_nonzero": True,
            "all_coordinate_equations_sha256": raw_component_hash(raw),
            "target_space_dimension": PORTS * math.comb(PORTS, 3),
            "target_invariant_character_certificate": target_character,
            "signed_target_orbit_count": len(orbits),
            "signed_target_orbit_sizes": [len(orbit) for orbit in orbits],
            "representative_equations": representatives,
            "representative_coefficient_rank_over_Q": len(full_rref),
            "representative_rref_pivot_monomial_indices": full_pivots,
            "reduced_integer_equations": [encode_integer_equation(row) for row in reduced_integer],
            "reduced_integer_equation_count": len(reduced_integer),
            "zero_locus_equivalence": "all 2640 components vanish iff the 44 orbit representatives vanish iff the coefficient-row rank 38 reduced integer quadrics vanish",
        },
        "fixed_line_reduction": {
            "fixed_vector": [1] * PORTS,
            "contracted_coordinate_count": len(fixed_components),
            "contracted_coefficient_rank_over_Q": len(fixed_rref),
            "contracted_rref_pivot_monomial_indices": fixed_pivots,
            "contracted_integer_equations_in_x": [encode_integer_equation(row) for row in fixed_integer],
            "channel_decomposition": channel_metadata,
            "derivation_product_equations_in_channel_variables": [
                {
                    "channel_id": spec[0],
                    "weight_coefficients_on_d_plus_d_minus_d_five": list(spec[4]),
                    "equation": encode_q5_equation(row),
                }
                for spec, row in zip(CHANNEL_SPECS, products_y)
            ],
            "product_equation_rank_over_Qsqrt5": q5_rank(products_y),
            "product_rowspace_equals_contracted_jacobi_rowspace": True,
            "residual_integer_equations_in_x": [encode_integer_equation(row) for row in residual_integer],
            "residual_equations_in_channel_variables": [encode_q5_equation(row) for row in residual_y],
            "residual_rank_over_Qsqrt5": q5_rank(residual_y),
            "full_rank_split": [q5_rank(products_y), q5_rank(residual_y), q5_rank([*products_y, *residual_y])],
            "derivation_only_subspace": {
                "dimension": 3,
                "definition": "all eleven non-derivation channel variables are zero",
                "all_jacobi_equations_vanish_identically": True,
                "scope": DERIVATION_ONLY_SCOPE,
            },
            "central_specialization": {
                "definition": "d_plus=d_minus=d_five=0",
                "remaining_channel_variable_count": 11,
                "residual_coefficient_rank_over_Qsqrt5": central_rank,
                "solution_set_classified": False,
            },
            "derivation_weight_arrangement": {
                "distinct_weight_hyperplane_count": len(DISTINCT_WEIGHTS),
                "flat_count": len(flats),
                "flat_dimension_counts": {str(dimension): sum(row["dimension"] == dimension for row in flats) for dimension in range(4)},
                "relative_open_flats": flats,
                "scope": ARRANGEMENT_SCOPE,
            },
        },
    }
    validate_system_declarations(system)
    system_hash = object_sha256(system)
    receipt: dict[str, Any] = {
        "schema": "oph.a5_jacobi_stage2.receipt.v1",
        "issue": 566,
        "status": STATUS,
        "claim_boundary": CLAIM_BOUNDARY,
        "semantic_inputs": {
            "count": 2,
            "paths": system["upstream"],
            "stage1_later_gates_all_false": not any(stage1_receipt["later_gates"].values()),
        },
        "target_firewall": {
            "enabled": True,
            "accepted_input_schemas": [
                "oph.a5_alternating_bracket_reynolds_basis.v1",
                "oph.a5_alternating_bracket_space_stage1.receipt.v1",
            ],
            "exact_input_keysets_enforced": True,
            "arbitrary_input_path_allowed": False,
            "coefficient_target_input_allowed": False,
            "measurement_input_allowed": False,
            "environment_target_read": False,
            "network_target_read": False,
            "nonstdlib_imports": [],
            "audited_import_roots_by_file": imports,
            "exact_import_allowlists_by_file": {
                filename: sorted(roots)
                for filename, roots in sorted(EXACT_IMPORT_ROOTS_BY_FILE.items())
            },
            "runtime_forbidden_import_roots_by_file": {
                filename: sorted(roots)
                for filename, roots in sorted(RUNTIME_FORBIDDEN_IMPORT_ROOTS.items())
            },
        },
        "attained": {
            "general_parameter_count": PARAMETERS,
            "raw_jacobi_coordinate_count": len(raw),
            "equivariant_target_orbit_count": len(orbits),
            "independent_jacobi_quadrics": len(full_rref),
            "fixed_line_product_quadrics": len(fixed_rref),
            "residual_independent_quadrics": len(residual_y),
            "exact_channel_coordinate_change": True,
            "derivation_weight_relative_open_flat_count": len(flats),
            "generic_weight_stratum_forces_all_non_derivation_channels_zero": True,
            "derivation_only_three_parameter_family_certified": True,
        },
        "not_attained": {
            "full_jacobi_solution_variety_classified": False,
            "special_weight_flat_residual_solutions_classified": False,
            "rational_descent_components_classified": False,
            "preferred_bracket_selected": False,
            "source_selection": False,
            "compactness_established": False,
            "issue_566_closed": False,
        },
        "system_artifact": {
            "path": "code/a5_closure/issue_566_bracket_space_stage2/a5_jacobi_system_reduction.json",
            "canonical_json_sha256": system_hash,
        },
        "mutation_tests": mutation_tests,
        "implementation_pins": {
            "producer_sha256": file_sha256(Path(__file__).resolve()),
            "independent_verifier_sha256": file_sha256(VERIFIER_PATH.resolve()),
            "test_sha256": file_sha256(TEST_PATH.resolve()),
        },
    }
    rehash_receipt(receipt)
    receipt["mutation_tests"].extend(run_artifact_mutation_tests(system, receipt, imports))
    rehash_receipt(receipt)
    validate_receipt_declarations(system, receipt, imports)
    return system, receipt


def serialized(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="check committed artifacts without writing")
    args = parser.parse_args(argv)
    system, receipt = build_artifacts()
    system_text, receipt_text = serialized(system), serialized(receipt)
    if args.check:
        require(SYSTEM_PATH.read_text(encoding="utf-8") == system_text, "ARTIFACT_STALE", "Jacobi system artifact is stale")
        require(RECEIPT_PATH.read_text(encoding="utf-8") == receipt_text, "ARTIFACT_STALE", "Jacobi receipt is stale")
    else:
        SYSTEM_PATH.write_text(system_text, encoding="utf-8")
        RECEIPT_PATH.write_text(receipt_text, encoding="utf-8")
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "raw_jacobi_coordinates": receipt["attained"]["raw_jacobi_coordinate_count"],
                "equivariant_orbits": receipt["attained"]["equivariant_target_orbit_count"],
                "independent_quadrics": receipt["attained"]["independent_jacobi_quadrics"],
                "rank_split": [
                    receipt["attained"]["fixed_line_product_quadrics"],
                    receipt["attained"]["residual_independent_quadrics"],
                ],
                "full_classification_open": not receipt["not_attained"]["full_jacobi_solution_variety_classified"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CertificateError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
