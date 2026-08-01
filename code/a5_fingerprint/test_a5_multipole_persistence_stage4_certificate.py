from __future__ import annotations

import json

import a5_multipole_fixed_point_certificate as base
import a5_multipole_persistence_stage4_certificate as stage4


def test_stage4_receipt_is_byte_exact_and_promotes_only_the_math() -> None:
    receipt = stage4.build_receipt()
    assert stage4.RECEIPT_PATH.read_bytes() == base.canonical_json_bytes(receipt)
    assert receipt["status"] == stage4.STATUS
    theorem = receipt["persistence_theorem"]
    assert theorem["finite_exactly_62_persistence_range"] is True
    assert theorem["range"] == "0 < |a k| <= 1"
    assert theorem["oriented_stationary_direction_count"] == 62
    assert theorem["source_selection_proved"] is False
    assert theorem["physical_attachment_proved"] is False
    assert receipt["comparison_boundary"] == {
        "public_measurement_read": False,
        "comparison_permitted": False,
        "physical_map_open": True,
    }


def test_stage4_partition_and_local_margins_are_strict() -> None:
    receipt = json.loads(stage4.RECEIPT_PATH.read_text())
    cover = receipt["global_off_neighborhood_cover"]
    assert cover["leaf_count"] == 2624
    assert max(row["maximum_depth"] for row in cover["chart_rows"]) == 19
    assert cover["all_off_neighborhood_boxes_gradient_certified"] is True
    local = receipt["local_uniqueness_boxes"]
    assert local["every_local_neighborhood_has_exactly_one_stationary_direction"] is True
    assert [row["orbit"] for row in local["orbit_representatives"]] == [
        "vertex",
        "face",
        "edge",
    ]
    assert all(int(row["strict_margin"].split("/")[0]) > 0 for row in local["orbit_representatives"])


def test_stage4_producer_controls_and_stabilizers() -> None:
    receipt = json.loads(stage4.RECEIPT_PATH.read_text())
    controls = receipt["fail_closed_controls"]["controls"]
    assert len(controls) == 4
    assert all(row["detector_fired"] is True for row in controls)
    stabilizers = receipt["exact_axis_stationarity"]["representative_rotations"]
    assert [(row["orbit"], row["rotation_order"]) for row in stabilizers] == [
        ("vertex", 5),
        ("face", 3),
        ("edge", 2),
    ]
    assert all(row["tangent_fixed_space_dimension"] == 0 for row in stabilizers)
    transport = receipt["exact_axis_stationarity"]["full_group_transport"]
    assert transport["proper_rotation_count"] == 60
    assert len(transport["port_permutations"]) == 60
    assert [
        (row["orbit"], row["projective_orbit_size"])
        for row in transport["projective_axis_orbits"]
    ] == [("vertex", 6), ("face", 10), ("edge", 15)]
    assert transport["projective_orbit_union_axis_count"] == 31
    assert (
        transport["minimum_distinct_projective_axis_sine_squared"]
        == "1/2+-1/6*sqrt5"
    )
    assert transport["minimum_separation_exceeds_twice_local_radius"] is True
    join = receipt["local_uniqueness_boxes"]["symmetry_transport_join"]
    assert join["port_permutation_sha256"] == transport["port_permutation_sha256"]
    assert join["transported_local_box_count"] == 31
    assert join["neighborhoods_sine_squared_1_over_4096_pairwise_disjoint"] is True
