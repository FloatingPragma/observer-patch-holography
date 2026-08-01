#!/usr/bin/env python3
"""Independent semantic verifier for the issue-656 refinement packet.

This verifier deliberately does not import the producer.  It recomputes the
exact arithmetic, group census, path bounds, rotating covariance fixture, and
mesh formulas through a separate implementation and fails closed on semantic
or custody drift.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Sequence


HERE = Path(__file__).resolve().parent
DEFAULT_RECEIPT = HERE / "runtime" / "discrete_refinement_theorem_receipt.json"
REPO = HERE.parent.parent
SCHEMA = "oph.discrete_refinement_theorem_packet.v1"
STATUS = (
    "EXACT_DISCRETE_REFINEMENT_THEOREMS_AND_DIVISIBILITY_MESH__"
    "PHYSICAL_PRODUCER_MISSING"
)
EXIT = (
    "THEOREM_PACKET_AND_MESH_SCAFFOLD_ATTAINED__"
    "SOURCE_NATIVE_PHYSICAL_BIREFINEMENT_OPEN"
)
MUTATION_IDS = (
    "dr1a_lattice_ratio",
    "dr2_nonresonance_difference",
    "dr3_path_domain",
    "dr3_pell_identity",
    "dr4_rotation_action",
    "finite_group_order_census",
    "mesh_vertex_count",
    "mesh_refinement_square",
    "observer_level_scope_promotion",
    "external_numeric_ancestor",
    "lean_digest_pin",
)


class VerificationError(ValueError):
    """The independent verifier rejected the packet."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def tagged_sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def body_digest(payload: dict[str, Any]) -> str:
    body = {key: value for key, value in payload.items() if key != "receipt_sha256"}
    return tagged_sha256(canonical_json_bytes(body))


def parse_fraction(value: Any) -> Fraction:
    require(isinstance(value, (str, int)), f"nonexact rational value: {value!r}")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise VerificationError(f"invalid exact rational: {value!r}") from error


def independent_profile(coordinate: Fraction) -> Fraction:
    floor = math.floor(coordinate)
    cell = coordinate - floor
    return Fraction(1) + min(cell, 1 - cell) / 5


def independent_lattice_value(coordinate: Fraction, multiplier: Fraction) -> Fraction:
    return multiplier ** math.floor(coordinate) * independent_profile(coordinate)


def verify_one_ratio(payload: dict[str, Any]) -> None:
    section = payload["one_ratio"]
    dr1 = section["DR-1"]
    control = dr1["exact_nonconstant_control"]
    multiplier = parse_fraction(control["multiplier"])
    require(multiplier > 0, "DR-1 multiplier is not positive")
    require(control["profile_positive"] is True, "DR-1 positivity flag drift")
    require(control["profile_nonconstant"] is True, "DR-1 nonconstant flag drift")
    rows = control["scale_equation_rows"]
    require(len(rows) == 31, "DR-1 grid row count drift")
    for row in rows:
        coordinate = parse_fraction(row["log_coordinate"])
        at_x = parse_fraction(row["F_at_x"])
        at_next = parse_fraction(row["F_at_x_plus_one"])
        require(at_x == independent_lattice_value(coordinate, multiplier), "DR-1 F row drift")
        require(
            at_next == independent_lattice_value(coordinate + 1, multiplier),
            "DR-1 shifted F row drift",
        )
        require(at_next == multiplier * at_x, "DR-1 scale equation failure")

    dr1a = section["DR-1A"]
    require("every integer" in dr1a["statement"], "DR-1A integer scope missing")
    require("every scale move" in dr1a["path_condition"], "DR-1A path boundary missing")
    lattice_rows = dr1a["exact_lattice_rows"]
    require(len(lattice_rows) == 27, "DR-1A lattice row count drift")
    for row in lattice_rows:
        anchor = parse_fraction(row["anchor"])
        exponent = row["integer_exponent"]
        require(isinstance(exponent, int) and -4 <= exponent <= 4, "DR-1A exponent drift")
        ratio = parse_fraction(row["ratio"])
        independent = independent_lattice_value(anchor + exponent, multiplier) / independent_lattice_value(
            anchor, multiplier
        )
        require(ratio == independent == multiplier**exponent, "DR-1A ratio failure")


def permutation_parity(permutation: Sequence[int]) -> int:
    return sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    ) % 2


def permutation_order(permutation: Sequence[int]) -> int:
    visited: set[int] = set()
    result = 1
    for start in range(len(permutation)):
        if start in visited:
            continue
        cursor = start
        cycle = 0
        while cursor not in visited:
            visited.add(cursor)
            cycle += 1
            cursor = permutation[cursor]
        result = math.lcm(result, cycle)
    return result


def verify_finite_group_no_go(payload: dict[str, Any]) -> None:
    section = payload["finite_group_scale_no_go"]
    members = [
        permutation
        for permutation in itertools.permutations(range(5))
        if permutation_parity(permutation) == 0
    ]
    census: dict[str, int] = {}
    for member in members:
        key = str(permutation_order(member))
        census[key] = census.get(key, 0) + 1
    require(len(members) == section["exact_group_order"] == 60, "A5 cardinality failure")
    require(
        census == {"1": 1, "2": 15, "3": 20, "5": 6 * 4},
        "unexpected A5 census",
    )
    enumeration = [
        {"permutation": list(member), "order": permutation_order(member)}
        for member in members
    ]
    require(section["element_order_support"] == [1, 2, 3, 5], "A5 order support drift")
    require(
        section["enumeration_sha256"] == tagged_sha256(canonical_json_bytes(enumeration)),
        "A5 enumeration digest failure",
    )
    require("strictly positive reals" in section["statement"], "scale no-go codomain drift")
    require("separate scale law" in section["consequence"], "scale no-go boundary drift")


def verify_bi_refinement(payload: dict[str, Any]) -> None:
    dr2 = payload["bi_refinement"]["DR-2"]
    require("same physical one-dimensional" in dr2["required_same_object"], "DR-2 scalar-ray premise missing")
    require(len(dr2["proof_contract"]) == 4, "DR-2 proof contract drift")
    require(
        "explicit premises" in dr2["lean_formalization_boundary"]
        and "continuity" in dr2["lean_formalization_boundary"],
        "DR-2 Lean boundary drift",
    )
    incommensurate = dr2["binary_ternary_incommensurability"]
    require(incommensurate["statement"] == "log(2)/log(3) is irrational", "DR-2 ratio drift")
    rows = incommensurate["finite_executable_rows"]
    require(len(rows) == 144, "DR-2 nonresonance census drift")
    seen: set[tuple[int, int]] = set()
    for row in rows:
        exponent_two = row["binary_exponent"]
        exponent_three = row["ternary_exponent"]
        require(1 <= exponent_two <= 12 and 1 <= exponent_three <= 12, "DR-2 exponent range drift")
        difference = 2**exponent_two - 3**exponent_three
        require(difference != 0 and difference == row["difference"], "DR-2 resonance row failure")
        seen.add((exponent_two, exponent_three))
    require(len(seen) == 144, "DR-2 duplicate nonresonance rows")

    control = dr2["exact_pure_power_control"]
    exponent = control["theta"]
    lambda_two = parse_fraction(control["lambda_at_ratio_two"])
    lambda_three = parse_fraction(control["lambda_at_ratio_three"])
    require(exponent == 2 and lambda_two == Fraction(1, 4), "DR-2 binary control drift")
    require(lambda_three == Fraction(1, 9), "DR-2 ternary control drift")
    require(len(control["rows"]) == 4, "DR-2 pure-power row count drift")
    for row in control["rows"]:
        scale = parse_fraction(row["scale"])
        value = parse_fraction(row["F"])
        twice = parse_fraction(row["F_at_twice_scale"])
        thrice = parse_fraction(row["F_at_thrice_scale"])
        require(value == scale ** (-exponent), "DR-2 pure-power value drift")
        require(twice == lambda_two * value, "DR-2 binary eigenvalue drift")
        require(thrice == lambda_three * value, "DR-2 ternary eigenvalue drift")
    require("global rigidity" in dr2["finite_window_boundary"], "DR-2 domain boundary missing")


def verify_translation_path(record: dict[str, Any]) -> None:
    domain = record["certified_domain"]
    horizontal_domain = domain["horizontal"]
    vertical_domain = domain["vertical"]
    require(horizontal_domain == [-3, 3] and vertical_domain == [-3, 3], "path domain drift")
    nodes = record["complete_path"]
    edges = record["edges"]
    require(len(nodes) == len(edges) + 1 and nodes[0] == [0, 0], "path chain shape drift")
    require(nodes[-1] == record["endpoint"], "path endpoint drift")
    error_sum = Fraction(0)
    bound_sum = Fraction(0)
    for index, edge in enumerate(edges):
        require(edge["from"] == nodes[index] and edge["to"] == nodes[index + 1], "disconnected path edge")
        start, finish = edge["from"], edge["to"]
        delta = (finish[0] - start[0], finish[1] - start[1])
        require(delta in {(1, 0), (-1, 0), (0, 1), (0, -1)}, "nonunit graph translation")
        observed = parse_fraction(edge["observed_increment"])
        nominal = parse_fraction(edge["nominal_increment"])
        error = parse_fraction(edge["edge_error"])
        limit = parse_fraction(edge["edge_error_limit"])
        require(observed - nominal == error, "edge residual arithmetic drift")
        require(abs(error) <= limit, "edge residual exceeds its certificate")
        error_sum += error
        bound_sum += limit
    require(
        all(
            horizontal_domain[0] <= node[0] <= horizontal_domain[1]
            and vertical_domain[0] <= node[1] <= vertical_domain[1]
            for node in nodes
        ),
        "translation path leaves certified domain",
    )
    require(record["path_stays_inside_domain"] is True, "path condition flag drift")
    require(parse_fraction(record["actual_minus_nominal"]) == error_sum, "path residual drift")
    require(parse_fraction(record["triangle_bound"]) == bound_sum, "path triangle bound drift")
    require(abs(error_sum) <= bound_sum, "path theorem violated")


def verify_approximate_rigidity(payload: dict[str, Any]) -> None:
    dr3 = payload["approximate_rigidity"]["DR-3"]
    require("epsilon/" in dr3["fourier_mode_bound"], "DR-3 Fourier denominator missing")
    require("tail" in dr3["stability_boundary"], "DR-3 smoothness boundary missing")
    graph = dr3["finite_translation_graph"]
    require("every intermediate" in graph["mandatory_path_condition"], "DR-3 path condition missing")
    records = graph["exact_path_records"]
    require([record["endpoint"] for record in records] == [[2, 1], [-2, 3], [3, -2]], "path endpoints drift")
    for record in records:
        verify_translation_path(record)

    control = dr3["small_divisor_negative_control"]
    require(control["irrational_rotation"] == "sqrt(2)", "small-divisor rotation drift")
    rows = control["rows"]
    require(len(rows) == 9, "Pell row count drift")
    previous_pair: tuple[int, int] | None = None
    prior_pair: tuple[int, int] | None = None
    prior_bound: Fraction | None = None
    fixed_amplitude: Fraction | None = None
    for index, row in enumerate(rows):
        require(row["index"] == index, "Pell index drift")
        pair = (row["p"], row["q"])
        pell = pair[0] ** 2 - 2 * pair[1] ** 2
        require(pell == row["p_squared_minus_two_q_squared"] and abs(pell) == 1, "Pell identity failure")
        if prior_pair is not None and previous_pair is not None:
            require(
                pair == (2 * previous_pair[0] + prior_pair[0], 2 * previous_pair[1] + prior_pair[1]),
                "Pell recurrence drift",
            )
        bound = parse_fraction(row["abs_q_sqrt_two_minus_p_upper"])
        require(bound == Fraction(1, pair[0] + pair[1]), "Pell upper bound drift")
        if prior_bound is not None:
            require(bound < prior_bound, "small-divisor bounds fail to shrink")
        peak = parse_fraction(row["fixed_log_profile_peak_to_trough"])
        if fixed_amplitude is None:
            fixed_amplitude = peak
        require(peak == fixed_amplitude == Fraction(2, 5), "negative-control amplitude drift")
        require(
            parse_fraction(row["log_profile_residual_over_two_pi_upper"])
            == Fraction(1, 5) * bound,
            "negative-control residual bound drift",
        )
        prior_pair, previous_pair = previous_pair, pair
        prior_bound = bound


Matrix = tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]
IDENTITY: Matrix = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)))


def parse_matrix(value: Any) -> Matrix:
    require(isinstance(value, list) and len(value) == 2, "matrix row count drift")
    require(all(isinstance(row, list) and len(row) == 2 for row in value), "matrix column count drift")
    return tuple(tuple(parse_fraction(entry) for entry in row) for row in value)  # type: ignore[return-value]


def multiply(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(sum((left[i][k] * right[k][j] for k in range(2)), Fraction(0)) for j in range(2))
        for i in range(2)
    )  # type: ignore[return-value]


def transpose(matrix: Matrix) -> Matrix:
    return ((matrix[0][0], matrix[1][0]), (matrix[0][1], matrix[1][1]))


def power(matrix: Matrix, exponent: int) -> Matrix:
    result = IDENTITY
    for _ in range(exponent):
        result = multiply(result, matrix)
    return result


def verify_scalar_ray(payload: dict[str, Any]) -> None:
    dr4 = payload["scalar_ray"]["DR-4"]
    require(len(dr4["producer_checklist"]) == 5, "DR-4 checklist length drift")
    counter = dr4["exact_rotating_countermodel"]
    first = parse_matrix(counter["first_action"])
    second = parse_matrix(counter["second_action"])
    base = parse_matrix(counter["base_covariance"])
    require(multiply(first, second) == multiply(second, first), "DR-4 actions fail to commute")
    require(multiply(first, transpose(first)) == IDENTITY, "DR-4 first action not orthogonal")
    require(multiply(second, transpose(second)) == IDENTITY, "DR-4 second action not orthogonal")
    require(multiply(first, first) == ((-Fraction(1), Fraction(0)), (Fraction(0), -Fraction(1))), "DR-4 no-ray witness drift")
    require(counter["actions_commute"] is True and counter["actions_orthogonal"] is True, "DR-4 flags drift")
    require(counter["positive_spectrum"] == ["1", "2"], "DR-4 spectrum drift")
    rows = counter["rows"]
    require(len(rows) == 9, "DR-4 covariance row count drift")
    for row in rows:
        m, n = row["binary_level"], row["ternary_level"]
        rotation = multiply(power(first, m), power(second, n))
        covariance = multiply(multiply(rotation, base), transpose(rotation))
        require(parse_matrix(row["covariance"]) == covariance, "DR-4 covariance row drift")
        trace = covariance[0][0] + covariance[1][1]
        determinant = covariance[0][0] * covariance[1][1] - covariance[0][1] * covariance[1][0]
        require(parse_fraction(row["trace"]) == trace == 3, "DR-4 trace drift")
        require(parse_fraction(row["determinant"]) == determinant == 2, "DR-4 determinant drift")


def verify_mesh(payload: dict[str, Any]) -> None:
    tower = payload["icosahedral_divisibility_tower"]
    require("mesh combinatorics only" in tower["scope"], "mesh/physics type boundary missing")
    require(tower["physics_promotion_allowed"] is False, "mesh illegally promotes physics")
    require(tower["base_automorphism_counts"] == {"full": 120, "proper": 60}, "mesh action census drift")
    records = tower["mesh_records"]
    require([record["n"] for record in records] == [1, 2, 3, 6], "mesh frequency rows drift")
    for record in records:
        denominator = record["n"]
        faces = 20 * denominator * denominator
        edges = 30 * denominator * denominator
        vertices = 12 + 30 * (denominator - 1) + 20 * (denominator - 1) * (denominator - 2) // 2
        require(vertices == 10 * denominator * denominator + 2, "independent vertex formula mismatch")
        require(record["vertices"] == vertices, "mesh vertex count drift")
        require(record["edges"] == edges, "mesh edge count drift")
        require(record["faces"] == faces, "mesh face count drift")
        require(vertices - edges + faces == record["euler_characteristic"] == 2, "mesh Euler drift")
        require(record["proper_action_order"] == 60, "mesh proper action drift")
        require(record["proper_action_preserves_vertices_and_faces"] is True, "mesh equivariance flag drift")
        require(
            isinstance(record["mesh_sha256"], str)
            and record["mesh_sha256"].startswith("sha256:")
            and len(record["mesh_sha256"]) == 71,
            "mesh digest shape drift",
        )

    morphisms = tower["refinement_morphisms"]
    require(morphisms["commuting_square"] == "R2 R3 = R3 R2 = R6", "mesh square statement drift")
    rows = []
    for denominator in range(1, 7):
        for i in range(denominator + 1):
            for j in range(denominator + 1 - i):
                k = denominator - i - j
                expected = [6 * i, 6 * j, 6 * k]
                rows.append(
                    {
                        "source_denominator": denominator,
                        "source_coordinate": [i, j, k],
                        "binary_then_ternary": expected,
                        "ternary_then_binary": expected,
                    }
                )
    expected_rows = sum((n + 1) * (n + 2) // 2 for n in range(1, 7))
    require(len(rows) == expected_rows, "mesh square row census drift")
    require(
        morphisms["verified_source_denominators"] == [1, 2, 3, 4, 5, 6],
        "mesh source-denominator inventory drift",
    )
    require(morphisms["exact_square_row_count"] == expected_rows, "mesh square count drift")
    require(
        morphisms["exact_square_sha256"] == tagged_sha256(canonical_json_bytes(rows)),
        "mesh square digest drift",
    )
    seen: set[tuple[int, int, int, int]] = set()
    for row in rows:
        denominator = row["source_denominator"]
        coordinate = row["source_coordinate"]
        require(len(coordinate) == 3 and sum(coordinate) == denominator, "mesh barycentric row drift")
        expected = [6 * value for value in coordinate]
        require(row["binary_then_ternary"] == expected, "R2R3 coordinate drift")
        require(row["ternary_then_binary"] == expected, "R3R2 coordinate drift")
        seen.add((denominator, *coordinate))
    require(len(seen) == expected_rows, "duplicate mesh square rows")


def independent_numeric_leaf_paths(value: Any, target: int, prefix: str = "") -> list[str]:
    if isinstance(value, bool):
        return []
    if isinstance(value, int):
        return [prefix] if value == target else []
    if isinstance(value, dict):
        return [
            path
            for key, child in value.items()
            for path in independent_numeric_leaf_paths(
                child, target, f"{prefix}.{key}" if prefix else key
            )
        ]
    if isinstance(value, list):
        return [
            path
            for index, child in enumerate(value)
            for path in independent_numeric_leaf_paths(child, target, f"{prefix}[{index}]")
        ]
    return []


def verify_numeric_ancestry(payload: dict[str, Any]) -> None:
    scientific_sections = {
        key: payload[key]
        for key in (
            "one_ratio",
            "bi_refinement",
            "approximate_rigidity",
            "scalar_ray",
            "finite_group_scale_no_go",
            "icosahedral_divisibility_tower",
        )
    }
    quarantined_integer = int("2" + "4")
    expected = independent_numeric_leaf_paths(scientific_sections, quarantined_integer)
    require(
        expected == [],
        "target-like decimal leaked into serialized scientific payload",
    )
    audit = payload["numeric_ancestry_audit"]
    require("semantic numeric ancestry" in audit["audit_scope"], "ancestry scope drift")
    require(audit["lexical_absence_claimed"] is False, "false lexical absence claim introduced")
    require(audit["external_numeric_ancestors"] == [], "external numeric ancestor introduced")
    require(
        audit["comparison_or_calibration_ancestors"] == [],
        "comparison/calibration ancestor introduced",
    )
    collision = audit["derived_decimal_collision"]
    require(collision["value_name"] == "twenty-four", "derived collision label drift")
    require(collision["numeric_leaf_paths"] == expected, "derived collision path drift")
    require(collision["numeric_leaf_count"] == len(expected), "derived collision count drift")
    require(collision["imported_external_constant"] is False, "derived collision retyped as import")
    require(collision["serialized_as_decimal_leaf"] is False, "derived collision serialized")
    require(
        "order-five" in collision["derivation"] and "barycentric" in collision["derivation"],
        "derived collision provenance missing",
    )


def verify_lean_binding(payload: dict[str, Any]) -> None:
    binding = payload["lean_kernel_binding"]
    require(binding["path"] == "Lean/Screen/DiscreteRefinement.lean", "Lean path drift")
    path = REPO / binding["path"]
    raw = path.read_bytes()
    require(binding["bytes"] == len(raw), "Lean byte count drift")
    require(binding["sha256"] == tagged_sha256(raw), "Lean source digest drift")
    text = raw.decode("utf-8")
    expected = [
        "lattice_law_nat",
        "periodic_normal_form_scales",
        "positive_bounded_shift_multiplier_eq_one",
        "bi_refinement_periodic_factor_constant",
        "bi_refinement_multiplier_and_shape",
        "bi_refinement_normal_form_is_pure_power",
        "two_pow_ne_three_pow",
        "finite_group_to_positive_reals_trivial",
        "tower_counts_six",
        "binary_ternary_square",
    ]
    require(binding["theorems"] == expected, "Lean theorem inventory drift")
    for name in expected:
        require(f"theorem {name}" in text, f"Lean theorem absent: {name}")
    require("sorry" not in text and binding["contains_sorry"] is False, "Lean sorry boundary failed")
    require(
        "multiplier to one" in binding["dr2_formal_scope"]
        and "two exact periods" in binding["dr2_formal_scope"],
        "formal DR-2 scope missing",
    )
    require("does not construct" in binding["physical_boundary"], "Lean physical boundary missing")


def verify_independent_binding(payload: dict[str, Any]) -> None:
    binding = payload["independent_verifier_binding"]
    require(
        binding["path"] == "code/refinement/verify_discrete_refinement_independent.py",
        "independent verifier path drift",
    )
    raw = Path(__file__).resolve().read_bytes()
    require(binding["bytes"] == len(raw), "independent verifier byte count drift")
    require(binding["sha256"] == tagged_sha256(raw), "independent verifier digest drift")
    require(binding["mutation_ids"] == list(MUTATION_IDS), "mutation inventory drift")


def verify_payload(payload: dict[str, Any], raw: bytes | None = None) -> None:
    require(isinstance(payload, dict), "receipt must be an object")
    if raw is not None:
        require(raw == canonical_json_bytes(payload), "receipt bytes are not canonical")
    require(payload.get("schema") == SCHEMA, "schema drift")
    require(payload.get("status") == STATUS, "status drift")
    require(payload.get("issue") == 656, "issue drift")
    require(payload.get("receipt_sha256") == body_digest(payload), "receipt self-digest drift")
    scope = payload["theorem_scope"]
    require(scope["mathematical_packet"] is True, "mathematical packet flag drift")
    for field in (
        "physical_refinement_ratio_emitted",
        "physical_covariance_eigenvalue_emitted",
        "observer_level_prediction_emitted",
        "comparison_data_read",
    ):
        require(scope[field] is False, f"scope field promoted: {field}")
    verify_one_ratio(payload)
    verify_bi_refinement(payload)
    verify_approximate_rigidity(payload)
    verify_scalar_ray(payload)
    verify_finite_group_no_go(payload)
    verify_mesh(payload)
    verify_numeric_ancestry(payload)
    verify_lean_binding(payload)
    verify_independent_binding(payload)
    require(payload.get("certified_exit") == EXIT, "certified exit drift")


def _mutate_dr1a(payload: dict[str, Any]) -> None:
    payload["one_ratio"]["DR-1A"]["exact_lattice_rows"][0]["ratio"] = "7/9"


def _mutate_dr2(payload: dict[str, Any]) -> None:
    payload["bi_refinement"]["DR-2"]["binary_ternary_incommensurability"][
        "finite_executable_rows"
    ][0]["difference"] += 1


def _mutate_path(payload: dict[str, Any]) -> None:
    payload["approximate_rigidity"]["DR-3"]["finite_translation_graph"][
        "exact_path_records"
    ][0]["complete_path"][1] = [9, 9]


def _mutate_pell(payload: dict[str, Any]) -> None:
    payload["approximate_rigidity"]["DR-3"]["small_divisor_negative_control"][
        "rows"
    ][2]["p"] += 1


def _mutate_matrix(payload: dict[str, Any]) -> None:
    payload["scalar_ray"]["DR-4"]["exact_rotating_countermodel"]["second_action"][0][0] = "2/5"


def _mutate_group(payload: dict[str, Any]) -> None:
    payload["finite_group_scale_no_go"]["enumeration_sha256"] = "sha256:" + "0" * 64


def _mutate_mesh_count(payload: dict[str, Any]) -> None:
    payload["icosahedral_divisibility_tower"]["mesh_records"][3]["vertices"] += 1


def _mutate_square(payload: dict[str, Any]) -> None:
    payload["icosahedral_divisibility_tower"]["refinement_morphisms"][
        "exact_square_sha256"
    ] = "sha256:" + "0" * 64


def _mutate_scope(payload: dict[str, Any]) -> None:
    payload["theorem_scope"]["observer_level_prediction_emitted"] = True


def _mutate_ancestry(payload: dict[str, Any]) -> None:
    payload["numeric_ancestry_audit"]["external_numeric_ancestors"].append("external")


def _mutate_lean_pin(payload: dict[str, Any]) -> None:
    payload["lean_kernel_binding"]["sha256"] = "sha256:" + "0" * 64


def run_mutation_suite(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Resign and reject semantic mutations through the sibling implementation."""

    mutations: tuple[tuple[str, Callable[[dict[str, Any]], None]], ...] = (
        (MUTATION_IDS[0], _mutate_dr1a),
        (MUTATION_IDS[1], _mutate_dr2),
        (MUTATION_IDS[2], _mutate_path),
        (MUTATION_IDS[3], _mutate_pell),
        (MUTATION_IDS[4], _mutate_matrix),
        (MUTATION_IDS[5], _mutate_group),
        (MUTATION_IDS[6], _mutate_mesh_count),
        (MUTATION_IDS[7], _mutate_square),
        (MUTATION_IDS[8], _mutate_scope),
        (MUTATION_IDS[9], _mutate_ancestry),
        (MUTATION_IDS[10], _mutate_lean_pin),
    )
    results = []
    for mutation_id, mutate in mutations:
        tampered = copy.deepcopy(payload)
        mutate(tampered)
        tampered["receipt_sha256"] = body_digest(tampered)
        try:
            verify_payload(tampered)
        except VerificationError as error:
            results.append(
                {
                    "mutation_id": mutation_id,
                    "rejected": True,
                    "reason": str(error),
                }
            )
        else:
            raise VerificationError(f"semantic mutation was accepted: {mutation_id}")
    require([row["mutation_id"] for row in results] == list(MUTATION_IDS), "mutation run drift")
    return results


def verify_path(path: Path) -> None:
    raw = path.read_bytes()
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError("receipt is invalid JSON") from error
    verify_payload(payload, raw)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", nargs="?", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--mutations", action="store_true")
    args = parser.parse_args(argv)
    raw = args.receipt.read_bytes()
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError("receipt is invalid JSON") from error
    verify_payload(payload, raw)
    if args.mutations:
        results = run_mutation_suite(payload)
        print(json.dumps(results, indent=2))
    print("DISCRETE_REFINEMENT_INDEPENDENT_VALID")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
