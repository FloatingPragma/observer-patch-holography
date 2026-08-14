from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

import named_law_horizon_bridge_verdict as verdict


def test_horizon_lane_short_circuits_without_overclaim() -> None:
    value = verdict.build()
    assert value["status"] == verdict.STATUS
    assert value["short_circuit"]["input_available"] is False
    assert value["short_circuit"]["horizon_identification_attempted"] is False
    scope = value["scope"]
    assert scope["declared_named_law_branch_only"] is True
    assert all(
        scope[key] is False
        for key in (
            "inhabited_gravity_carrier_rejected",
            "all_future_capacity_laws_excluded",
            "finite_de_sitter_identities_affected",
            "cosmological_value_emitted",
            "comparison_permitted",
        )
    )


def test_runtime_is_byte_exact() -> None:
    verdict.verify()
    assert verdict.OUTPUT_PATH.read_bytes() == verdict.canonical(verdict.build())


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("comparison_boundary", "named_law_comparison_promotable"), True),
        (("comparison_boundary", "cosmological_prediction"), True),
        (("exit_scope", "bounded_to_declared_attachment_class"), False),
        (("exit_scope", "all_possible_global_capacity_laws_excluded"), True),
    ],
)
def test_parent_promotion_mutations_fail_closed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    path: tuple[str, str],
    value: bool,
) -> None:
    parent = copy.deepcopy(verdict.load(verdict.CAPACITY_VERDICT))
    parent[path[0]][path[1]] = value
    mutated = tmp_path / "capacity.json"
    mutated.write_text(json.dumps(parent), encoding="utf-8")
    monkeypatch.setattr(verdict, "CAPACITY_VERDICT", mutated)
    with pytest.raises(ValueError):
        verdict.build()
