from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

import named_law_n_closure_verdict as verdict


def _write(tmp_path: Path, name: str, value: dict) -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def test_verdict_is_bounded_and_does_not_promote_a_prediction() -> None:
    value = verdict.build()
    assert value["status"] == verdict.STATUS
    assert value["exit_scope"] == {
        "bounded_to_declared_attachment_class": True,
        "all_possible_global_capacity_laws_excluded": False,
        "new_source_law_forbidden": False,
        "reason": value["exit_scope"]["reason"],
    }
    boundary = value["comparison_boundary"]
    assert boundary["retrospective_arithmetic_retained"] is True
    assert all(
        boundary[key] is False
        for key in (
            "selective_weight",
            "cosmological_prediction",
            "prospective_forecast",
            "horizon_interpretation",
            "named_law_comparison_promotable",
        )
    )
    assert all(selected is False for selected in value["branch_status"].values())


def test_runtime_is_byte_exact() -> None:
    verdict.verify()
    assert verdict.OUTPUT_PATH.read_bytes() == verdict.canonical_json_bytes(verdict.build())


@pytest.mark.parametrize(
    ("parent", "mutation"),
    [
        ("branch", lambda p: p["scope"].__setitem__("branch_selected", True)),
        (
            "branch",
            lambda p: p["scope"].__setitem__("global_capacity_derived", True),
        ),
        (
            "attachment",
            lambda p: p["comparison_boundary"].__setitem__(
                "numeric_cosmic_capacity_emitted", True
            ),
        ),
        (
            "attachment",
            lambda p: p["comparison_boundary"].__setitem__(
                "forecast_or_comparison_permitted", True
            ),
        ),
    ],
)
def test_parent_promotion_mutations_fail_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, parent: str, mutation
) -> None:
    if parent == "branch":
        changed = copy.deepcopy(verdict._load(verdict.BRANCH_PACKET))
        mutation(changed)
        monkeypatch.setattr(
            verdict, "BRANCH_PACKET", _write(tmp_path, "branch.json", changed)
        )
    else:
        changed = copy.deepcopy(verdict._load(verdict.ATTACHMENT_PACKET))
        mutation(changed)
        monkeypatch.setattr(
            verdict, "ATTACHMENT_PACKET", _write(tmp_path, "attachment.json", changed)
        )
    with pytest.raises((ValueError, SystemExit)):
        verdict.build()
