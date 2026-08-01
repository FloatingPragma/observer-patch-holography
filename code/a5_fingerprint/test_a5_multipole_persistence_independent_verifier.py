from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any, Callable

import pytest

import a5_multipole_fixed_point_certificate as base
import a5_multipole_persistence_independent_verifier as verifier


Mutation = Callable[[dict[str, Any]], None]


def mutated_receipt(
    tmp_path: Path,
    mutation: Mutation,
    *,
    rehash: bool = True,
) -> Path:
    receipt = json.loads(verifier.RECEIPT_PATH.read_bytes())
    mutation(receipt)
    if rehash:
        receipt["receipt_sha256"] = verifier.self_hash(receipt)
    path = tmp_path / "mutated-v3.json"
    path.write_bytes(base.canonical_json_bytes(receipt))
    return path


def test_independent_verifier_recomputes_every_v3_surface() -> None:
    result = verifier.verify_receipt()
    assert result["status"] == "PASS"
    assert result["independent_template_recomputed"] is True
    assert result["independent_tail_bounds_recomputed"] is True
    assert result["independent_axis_separation_recomputed"] is True
    assert result["promotion_boundaries_fail_closed"] is True


def test_self_hash_mutation_fails_before_semantic_replay(tmp_path: Path) -> None:
    path = mutated_receipt(
        tmp_path,
        lambda row: row["through_eighth_order_i6_template"].__setitem__(
            "x8_subtracted_coefficient", "3/118125"
        ),
        rehash=False,
    )
    with pytest.raises(verifier.PersistenceVerificationError, match="self-digest"):
        verifier.verify_receipt(path)


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        (
            lambda row: row["through_eighth_order_i6_template"].__setitem__(
                "x8_subtracted_coefficient", "3/118125"
            ),
            "template mismatch",
        ),
        (
            lambda row: row["normalized_tail_bounds"]["derivation_parameters"].__setitem__(
                "gradient_geometric_sum_factor", "111/109"
            ),
            "tail-bound mismatch",
        ),
        (
            lambda row: row["normalized_tail_bounds"]["derivation_parameters"].__setitem__(
                "euclidean_hessian_geometric_sum_factor", "91/89"
            ),
            "tail-bound mismatch",
        ),
        (
            lambda row: row["normalized_tail_bounds"]["derivation_parameters"].__setitem__(
                "normalization_lower_coefficient", "59/118125"
            ),
            "tail-bound mismatch",
        ),
        (
            lambda row: row["normalized_tail_bounds"]["derivation_parameters"].__setitem__(
                "radial_hessian_correction_included", False
            ),
            "tail-bound mismatch",
        ),
        (
            lambda row: row["critical_axis_separation"].__setitem__(
                "declared_local_neighborhood",
                "sin^2(angle to an axis) <= 1/4",
            ),
            "axis separation mismatch",
        ),
        (
            lambda row: row["critical_axis_separation"].__setitem__(
                "promotion_from_separation_permitted", True
            ),
            "axis separation mismatch",
        ),
        (
            lambda row: row["quantitative_persistence"].__setitem__(
                "finite_exactly_62_persistence_range", True
            ),
            "persistence boundary drift",
        ),
        (
            lambda row: row["through_eighth_order_i6_template"].__setitem__(
                "source_selected", True
            ),
            "template mismatch",
        ),
        (
            lambda row: row["comparison_boundary"].__setitem__(
                "comparison_permitted", True
            ),
            "comparison or physical-map boundary drift",
        ),
    ),
)
def test_semantically_rehashed_mutations_fail_independent_replay(
    tmp_path: Path,
    mutation: Mutation,
    message: str,
) -> None:
    path = mutated_receipt(tmp_path, mutation)
    with pytest.raises(verifier.PersistenceVerificationError, match=message):
        verifier.verify_receipt(path)


def test_parent_pin_is_recomputed_from_actual_parent_bytes(tmp_path: Path) -> None:
    parent = json.loads(verifier.PARENT_PATH.read_bytes())
    parent["critical_points"] = copy.deepcopy(parent["critical_points"])
    parent["critical_points"]["census"]["maxima"] = 13
    parent["receipt_sha256"] = verifier.self_hash(parent)
    parent_path = tmp_path / "mutated-v2.json"
    parent_path.write_bytes(base.canonical_json_bytes(parent))
    with pytest.raises(verifier.PersistenceVerificationError, match="parent pin drift"):
        verifier.verify_receipt(parent_path=parent_path)
