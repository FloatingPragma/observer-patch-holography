from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

import a5_multipole_fixed_point_certificate as base
import a5_multipole_fixed_point_hardening_certificate as hardening
import a5_multipole_persistence_certificate as persistence


def test_through_eighth_order_template_is_exact_and_typed() -> None:
    row = persistence.through_eighth_order_i6_template_certificate()
    assert row["declared_range"] == "0 < x <= 1"
    assert row["kernel_status"] == "DECLARED_EQUAL_WEIGHT_COSINE_BRANCH"
    assert row["declared_branch_premise"] is True
    assert row["source_selected"] is False
    assert row["architecture_forced"] is False
    assert row["physical_source_selection_owner"] == 655
    assert row["domain_assumptions"] == {
        "angular_argument": "n in S^2 with ||n|| = 1",
        "carrier_directions": "u_i in S^2 with ||u_i|| = 1, i = 1,...,12",
        "carrier_unit_norm_verified": True,
        "normalization_method": (
            "the exact Cartesian vertices have common squared norm "
            "5/2+1/2*sqrt5; each direction is divided by the square root "
            "of that common squared norm"
        ),
    }
    assert row["moment_decomposition"] == {
        "sum_i (u_i.n)^6": "12/7 + (64/175) I6(n)",
        "sum_i (u_i.n)^8": "4/3 + (256/375) I6(n)",
    }
    assert row["x6_coefficient"] == "4/7875"
    assert row["x8_subtracted_coefficient"] == "2/118125"
    assert row["positivity_certificate"]["strictly_positive"] is True


def test_normalized_tail_bounds_are_exact_at_canonical_endpoint() -> None:
    row = persistence.cosine_tail_bounds()
    assert row["x_max"] == "1"
    assert row["range"] == "0 < x <= 1"
    assert row["normalization_lower_bound"] == (
        "A(x) >= 58/118125 x^6 for 0 < x <= 1"
    )
    assert row["normalized_C1_gradient_bound"] == "6875/101152"
    assert row["normalized_C2_intrinsic_hessian_bound"] == "383125/562658"
    assert Fraction(row["normalized_C1_gradient_bound"]) < Fraction(7, 100)
    assert Fraction(row["normalized_C2_intrinsic_hessian_bound"]) < Fraction(7, 10)
    assert row["derivation_parameters"]["gradient_geometric_sum_factor"] == "110/109"
    assert row["derivation_parameters"]["euclidean_hessian_geometric_sum_factor"] == "90/89"
    assert row["derivation_parameters"]["radial_hessian_correction_included"] is True


def test_tail_bounds_support_every_declared_smaller_endpoint() -> None:
    row = persistence.cosine_tail_bounds(x_max=Fraction(1, 2))
    assert row["x_max"] == "1/2"
    assert row["range"] == "0 < x <= 1/2"
    assert row["normalization_lower_bound"] == (
        "A(x) >= 17/33750 x^6 for 0 < x <= 1/2"
    )
    assert row["normalized_C1_gradient_bound"] == "6875/1660288"
    assert row["normalized_C2_intrinsic_hessian_bound"] == "383125/9235352"


def test_critical_axis_separation_is_exact() -> None:
    row = persistence.critical_axis_separation()
    assert row["unoriented_axis_count"] == 31
    assert row["oriented_critical_point_count"] == 62
    assert row["axis_counts_by_orbit"] == {
        "vertex": 6,
        "face": 10,
        "edge": 15,
    }
    assert row["maximum_squared_axis_cosine"] == "1/2+1/6*sqrt5"
    assert row["minimum_squared_axis_sine"] == "1/2-sqrt5/6"
    assert row["maximizing_pair_orbit_types"] == ["edge", "face"]
    assert row["declared_local_neighborhood"].endswith("1/64")
    assert row["promotion_from_separation_permitted"] is False


def test_v3_parent_is_append_only_and_byte_exact() -> None:
    receipt = persistence.build_receipt()
    parent_bytes = persistence.PARENT_PATH.read_bytes()
    parent = json.loads(parent_bytes)
    assert receipt["parent_pin"]["sha256"] == base.tagged_sha256(parent_bytes)
    assert receipt["parent_pin"]["receipt_sha256"] == parent["receipt_sha256"]
    assert receipt["extends"]["schema"] == hardening.SCHEMA
    assert "supersedes" not in receipt
    assert "v1 and v2 remain immutable" in receipt["extends"]["reason"]


def test_semantically_rehashed_parent_tamper_fails(tmp_path: Path) -> None:
    parent = json.loads(persistence.PARENT_PATH.read_bytes())
    parent["critical_points"]["census"]["maxima"] = 13
    parent["receipt_sha256"] = persistence.artifact_self_hash(parent)
    path = tmp_path / "v2.json"
    path.write_bytes(base.canonical_json_bytes(parent))
    with pytest.raises(base.FingerprintError, match="byte-exact producer replay"):
        persistence.build_receipt(parent_path=path)


def test_range_port_count_and_moment_mutations_fail_closed() -> None:
    with pytest.raises(base.FingerprintError, match="x <= 1"):
        persistence.cosine_tail_bounds(x_max=Fraction(2))
    with pytest.raises(base.FingerprintError, match="twelve-port"):
        persistence.cosine_tail_bounds(port_count=13)
    with pytest.raises(base.FingerprintError, match="sixth-moment amplitude"):
        persistence.through_eighth_order_i6_template_certificate(
            t6_amplitude_override=Fraction(65, 175)
        )
    with pytest.raises(base.FingerprintError, match="eighth-moment amplitude"):
        persistence.through_eighth_order_i6_template_certificate(
            t8_amplitude_override=Fraction(257, 375)
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    (
        ({"gradient_geometric_factor_override": Fraction(111, 109)}, "gradient"),
        ({"hessian_geometric_factor_override": Fraction(91, 89)}, "Hessian"),
        (
            {"normalization_lower_coefficient_override": Fraction(59, 118125)},
            "normalization",
        ),
        ({"include_radial_hessian_correction": False}, "radial-gradient"),
    ),
)
def test_tail_derivation_mutations_fail_closed(kwargs: dict, message: str) -> None:
    with pytest.raises(base.FingerprintError, match=message):
        persistence.cosine_tail_bounds(**kwargs)


def test_on_sphere_geometry_mutation_breaks_axis_certificate() -> None:
    vertices = base.cartesian_vertices()
    antipode = hardening.derive_antipode(vertices)
    mutated = list(vertices)

    def rotate_x(vector: tuple[base.Q5, base.Q5, base.Q5]):
        x, y, z = vector
        return (
            x,
            base.q5_add(
                base.q5_scale(y, Fraction(3, 5)),
                base.q5_scale(z, Fraction(-4, 5)),
            ),
            base.q5_add(
                base.q5_scale(y, Fraction(4, 5)),
                base.q5_scale(z, Fraction(3, 5)),
            ),
        )

    mutated[0] = rotate_x(vertices[0])
    mutated[antipode[0]] = rotate_x(vertices[antipode[0]])
    with pytest.raises(base.FingerprintError):
        persistence.critical_axis_separation(mutated)


def test_neighborhood_and_promotion_mutations_fail_closed() -> None:
    with pytest.raises(base.FingerprintError, match="neighborhoods are not separated"):
        persistence.critical_axis_separation(
            neighborhood_sine_squared=base.q5(Fraction(1, 4))
        )
    with pytest.raises(base.FingerprintError, match="cannot promote"):
        persistence.critical_axis_separation(
            promote_separation_to_global_persistence=True
        )


def test_receipt_types_unfinished_persistence_fail_closed() -> None:
    receipt = persistence.build_receipt()
    assert persistence.STATUS.startswith("EXACT_THROUGH_EIGHTH_ORDER_I6_TEMPLATE")
    assert "full_cosine_coefficient" not in receipt
    row = receipt["quantitative_persistence"]
    assert row["declared_full_cosine_kernel_for_tail_bounds"] is True
    assert row["through_eighth_order_i6_template"] is True
    assert row["exact_C1_C2_tail_bounds"] is True
    assert row["exact_critical_axis_separation"] is True
    assert row["global_interval_gradient_cover"] is False
    assert row["local_interval_newton_uniqueness_boxes"] is False
    assert row["finite_exactly_62_persistence_range"] is False
    assert "does not infer" in row["claim_boundary"]


def test_committed_v3_receipt_is_byte_exact() -> None:
    assert persistence.RECEIPT_PATH.read_bytes() == base.canonical_json_bytes(
        persistence.build_receipt()
    )
