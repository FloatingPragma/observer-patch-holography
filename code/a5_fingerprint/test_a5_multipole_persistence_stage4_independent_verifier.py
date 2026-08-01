from __future__ import annotations

import copy

import pytest

import a5_multipole_fixed_point_certificate as base
import a5_multipole_persistence_stage4_independent_verifier as independent


def test_independent_replay_and_resigned_semantic_mutations() -> None:
    parent_raw, parent = independent.load(independent.PARENT_PATH)
    _, receipt = independent.load(independent.RECEIPT_PATH)
    cover = independent.global_cover()
    local = independent.local_rows()
    result = independent.verify_object(
        receipt,
        parent,
        parent_raw,
        expected_cover=cover,
        expected_local=local,
    )
    assert result["status"] == "PASS"
    assert result["independent_global_leaf_replay"] is True
    assert result["independent_local_margin_replay"] is True
    assert result["independent_full_group_transport_replay"] is True
    assert result["canonical_parent_bytes_verified"] is True
    assert independent.mutations(receipt, parent, parent_raw, cover, local) == [
        "leaf_digest",
        "local_margin",
        "range",
        "source_promotion",
        "comparison",
        "parent_pin",
        "morse_type",
        "orbit_count",
        "extra_field",
        "global_boolean",
        "local_boolean",
        "implementation_pin",
        "operator_bound",
        "y_chart_bound",
        "symmetry_group_order",
        "symmetry_port_permutation",
        "symmetry_projective_orbit_size",
        "symmetry_axis_union_count",
        "symmetry_minimum_separation",
        "symmetry_neighborhood_disjointness",
        "parent_source_promotion",
        "parent_architecture_promotion",
        "parent_physical_map_promotion",
        "parent_comparison_promotion",
        "parent_measured_target_injection",
    ]


def test_rehashed_parent_promotion_is_not_accepted_as_v3() -> None:
    _, parent = independent.load(independent.PARENT_PATH)
    candidate = copy.deepcopy(parent)
    candidate["through_eighth_order_i6_template"]["source_selected"] = True
    candidate["receipt_sha256"] = independent.self_hash(candidate)
    candidate_raw = base.canonical_json_bytes(candidate)
    with pytest.raises(independent.VerificationError, match="canonical parent"):
        independent.verify_canonical_parent(candidate, candidate_raw)
