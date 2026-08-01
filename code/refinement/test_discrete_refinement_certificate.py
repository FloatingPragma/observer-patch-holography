from __future__ import annotations

import copy
import re
from fractions import Fraction

import pytest

import discrete_refinement_certificate as producer
import verify_discrete_refinement_independent as sibling


@pytest.fixture(scope="module")
def receipt() -> dict:
    return producer.build_receipt()


def resign(payload: dict) -> None:
    payload["receipt_sha256"] = sibling.body_digest(payload)


def test_committed_receipt_replays_byte_exactly(receipt: dict) -> None:
    assert producer.RECEIPT_PATH.read_bytes() == producer.canonical_json_bytes(receipt)
    producer.verify_committed()


def test_independent_verifier_accepts_committed_packet(receipt: dict) -> None:
    sibling.verify_payload(receipt)
    sibling.verify_path(producer.RECEIPT_PATH)


def test_independent_resigned_mutation_suite_fails_closed(receipt: dict) -> None:
    binding = receipt["independent_verifier_binding"]
    assert binding["implementation_independent_of_producer_import"] is True
    assert binding["mutations_are_resigned_before_semantic_verification"] is True
    results = sibling.run_mutation_suite(receipt)
    assert [row["mutation_id"] for row in results] == binding["mutation_ids"]
    assert all(row["rejected"] is True for row in results)


def test_packet_does_not_promote_mesh_or_theorems_to_physics(receipt: dict) -> None:
    scope = receipt["theorem_scope"]
    assert scope["mathematical_packet"] is True
    assert scope["physical_refinement_ratio_emitted"] is False
    assert scope["physical_covariance_eigenvalue_emitted"] is False
    assert scope["observer_level_prediction_emitted"] is False
    assert scope["comparison_data_read"] is False
    tower = receipt["icosahedral_divisibility_tower"]
    assert tower["physics_promotion_allowed"] is False
    assert "mesh combinatorics only" in tower["scope"]


def test_one_ratio_control_is_nonconstant_and_exact(receipt: dict) -> None:
    control = receipt["one_ratio"]["DR-1"]["exact_nonconstant_control"]
    multiplier = Fraction(control["multiplier"])
    values = set()
    for row in control["scale_equation_rows"]:
        coordinate = Fraction(row["log_coordinate"])
        at_x = Fraction(row["F_at_x"])
        at_next = Fraction(row["F_at_x_plus_one"])
        assert at_next == multiplier * at_x
        values.add(producer.periodic_profile(coordinate))
    assert len(values) > 1


def test_one_ratio_lattice_law_includes_negative_exponents(receipt: dict) -> None:
    rows = receipt["one_ratio"]["DR-1A"]["exact_lattice_rows"]
    assert {row["integer_exponent"] for row in rows} == set(range(-4, 5))
    multiplier = Fraction(4, 5)
    for row in rows:
        assert Fraction(row["ratio"]) == multiplier ** row["integer_exponent"]


def test_binary_ternary_control_and_incommensurability(receipt: dict) -> None:
    dr2 = receipt["bi_refinement"]["DR-2"]
    assert dr2["binary_ternary_incommensurability"]["statement"] == (
        "log(2)/log(3) is irrational"
    )
    control = dr2["exact_pure_power_control"]
    assert control["theta"] == 2
    assert Fraction(control["lambda_at_ratio_two"]) == Fraction(1, 4)
    assert Fraction(control["lambda_at_ratio_three"]) == Fraction(1, 9)


def test_dr2_is_bound_to_kernel_checked_normalized_rigidity(receipt: dict) -> None:
    binding = receipt["lean_kernel_binding"]
    assert "positive_bounded_shift_multiplier_eq_one" in binding["theorems"]
    assert "bi_refinement_multiplier_and_shape" in binding["theorems"]
    assert "bi_refinement_periodic_factor_constant" in binding["theorems"]
    assert "bi_refinement_normal_form_is_pure_power" in binding["theorems"]
    assert binding["contains_sorry"] is False
    sibling.verify_lean_binding(receipt)


def test_dr3_paths_print_and_satisfy_domain_condition(receipt: dict) -> None:
    graph = receipt["approximate_rigidity"]["DR-3"]["finite_translation_graph"]
    assert "every intermediate node" in graph["mandatory_path_condition"]
    for record in graph["exact_path_records"]:
        sibling.verify_translation_path(record)
        assert record["path_stays_inside_domain"] is True
        assert abs(Fraction(record["actual_minus_nominal"])) <= Fraction(
            record["triangle_bound"]
        )


def test_small_divisor_control_keeps_amplitude_fixed(receipt: dict) -> None:
    rows = receipt["approximate_rigidity"]["DR-3"]["small_divisor_negative_control"]["rows"]
    bounds = [Fraction(row["abs_q_sqrt_two_minus_p_upper"]) for row in rows]
    assert all(left > right for left, right in zip(bounds, bounds[1:]))
    assert {Fraction(row["fixed_log_profile_peak_to_trough"]) for row in rows} == {
        Fraction(2, 5)
    }
    assert all(abs(row["p_squared_minus_two_q_squared"]) == 1 for row in rows)


def test_dr4_countermodel_has_commuting_actions_without_scalar_ray(receipt: dict) -> None:
    counter = receipt["scalar_ray"]["DR-4"]["exact_rotating_countermodel"]
    first = sibling.parse_matrix(counter["first_action"])
    second = sibling.parse_matrix(counter["second_action"])
    assert sibling.multiply(first, second) == sibling.multiply(second, first)
    assert sibling.multiply(first, first) == (
        (Fraction(-1), Fraction(0)),
        (Fraction(0), Fraction(-1)),
    )
    assert all(row["trace"] == "3" and row["determinant"] == "2" for row in counter["rows"])


def test_a5_positive_scale_no_go_has_exact_order_census(receipt: dict) -> None:
    no_go = receipt["finite_group_scale_no_go"]
    assert no_go["exact_group_order"] == 60
    assert no_go["element_order_support"] == [1, 2, 3, 5]
    assert no_go["enumeration_sha256"].startswith("sha256:")
    assert "only when it equals one" in no_go["proof"]


def test_colliding_integer_has_derived_group_ancestry_only(receipt: dict) -> None:
    audit = receipt["numeric_ancestry_audit"]
    assert audit["lexical_absence_claimed"] is False
    assert "semantic numeric ancestry" in audit["audit_scope"]
    assert audit["external_numeric_ancestors"] == []
    assert audit["comparison_or_calibration_ancestors"] == []
    collision = audit["derived_decimal_collision"]
    assert collision["value_name"] == "twenty-four"
    assert collision["numeric_leaf_count"] == 0
    assert collision["numeric_leaf_paths"] == []
    assert collision["imported_external_constant"] is False
    assert collision["serialized_as_decimal_leaf"] is False


@pytest.mark.parametrize(
    ("denominator", "expected"),
    [
        (1, (12, 30, 20)),
        (2, (42, 120, 80)),
        (3, (92, 270, 180)),
        (6, (362, 1080, 720)),
    ],
)
def test_mesh_counts_from_constructed_complex(
    denominator: int, expected: tuple[int, int, int]
) -> None:
    vertices, faces, edges = producer.build_mesh(denominator)
    assert (len(vertices), len(edges), len(faces)) == expected
    assert len(vertices) - len(edges) + len(faces) == 2


def test_mesh_action_and_refinement_square_are_exact(receipt: dict) -> None:
    assert len(producer.graph_automorphisms()) == 120
    assert len(producer.proper_automorphisms()) == 60
    vertices_one, _, _ = producer.build_mesh(1)
    vertices_two, _, _ = producer.build_mesh(2)
    vertices_three, _, _ = producer.build_mesh(3)
    vertices_six, _, _ = producer.build_mesh(6)
    assert vertices_one <= vertices_two <= vertices_six
    assert vertices_one <= vertices_three <= vertices_six
    rows = receipt["icosahedral_divisibility_tower"]["refinement_morphisms"][
        "exact_square_row_count"
    ]
    assert rows == len(producer.exact_refinement_square_rows())


def mutation_one_ratio(payload: dict) -> None:
    payload["one_ratio"]["DR-1A"]["exact_lattice_rows"][0]["ratio"] = "7/9"


def mutation_dr2(payload: dict) -> None:
    payload["bi_refinement"]["DR-2"]["binary_ternary_incommensurability"][
        "finite_executable_rows"
    ][0]["difference"] += 1


def mutation_path(payload: dict) -> None:
    payload["approximate_rigidity"]["DR-3"]["finite_translation_graph"][
        "exact_path_records"
    ][0]["complete_path"][1] = [9, 9]


def mutation_pell(payload: dict) -> None:
    payload["approximate_rigidity"]["DR-3"]["small_divisor_negative_control"][
        "rows"
    ][2]["p"] += 1


def mutation_matrix(payload: dict) -> None:
    payload["scalar_ray"]["DR-4"]["exact_rotating_countermodel"]["second_action"][0][0] = "2/5"


def mutation_group(payload: dict) -> None:
    payload["finite_group_scale_no_go"]["enumeration_sha256"] = "sha256:" + "0" * 64


def mutation_mesh_count(payload: dict) -> None:
    payload["icosahedral_divisibility_tower"]["mesh_records"][3]["vertices"] += 1


def mutation_square(payload: dict) -> None:
    payload["icosahedral_divisibility_tower"]["refinement_morphisms"][
        "exact_square_sha256"
    ] = "sha256:" + "0" * 64


def mutation_scope(payload: dict) -> None:
    payload["theorem_scope"]["observer_level_prediction_emitted"] = True


def mutation_ancestry(payload: dict) -> None:
    payload["numeric_ancestry_audit"]["external_numeric_ancestors"].append("external")


def mutation_lean_binding(payload: dict) -> None:
    payload["lean_kernel_binding"]["sha256"] = "sha256:" + "0" * 64


@pytest.mark.parametrize(
    "mutate",
    [
        mutation_one_ratio,
        mutation_dr2,
        mutation_path,
        mutation_pell,
        mutation_matrix,
        mutation_group,
        mutation_mesh_count,
        mutation_square,
        mutation_scope,
        mutation_ancestry,
        mutation_lean_binding,
    ],
    ids=lambda function: function.__name__,
)
def test_independent_verifier_rejects_semantic_mutations(receipt: dict, mutate) -> None:
    tampered = copy.deepcopy(receipt)
    mutate(tampered)
    resign(tampered)
    with pytest.raises(sibling.VerificationError):
        sibling.verify_payload(tampered)


def test_independent_verifier_rejects_stale_self_digest(receipt: dict) -> None:
    tampered = copy.deepcopy(receipt)
    tampered["status"] = "tampered"
    with pytest.raises(sibling.VerificationError, match="status drift|self-digest drift"):
        sibling.verify_payload(tampered)


def test_receipt_contains_no_comparison_payload(receipt: dict) -> None:
    def semantic_strings(value, key: str = "") -> list[str]:
        if "sha256" in key:
            return []
        if isinstance(value, str):
            return [value]
        if isinstance(value, dict):
            return [
                item
                for child_key, child in value.items()
                for item in semantic_strings(child, child_key)
            ]
        if isinstance(value, list):
            return [item for child in value for item in semantic_strings(child, key)]
        return []

    tokens = ["n" + "_" + "s", chr(80), "2" + "4", "4" + "8"]
    patterns = [
        re.compile(rf"(?<![A-Za-z0-9_]){re.escape(token)}(?![A-Za-z0-9_])")
        for token in tokens
    ]
    for text in semantic_strings(receipt):
        assert all(pattern.search(text) is None for pattern in patterns)
    assert receipt["theorem_scope"]["comparison_data_read"] is False
