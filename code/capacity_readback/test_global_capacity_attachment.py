from __future__ import annotations

import copy
import json

import pytest

import global_capacity_attachment as cert


@pytest.fixture(scope="module")
def payload() -> dict:
    return cert.build()


def test_two_lawful_global_completions_disagree(payload: dict) -> None:
    assert payload["status"] == cert.STATUS
    for row in payload["exact_witnesses"]:
        assert row["different_at_one_cut"] is True
        completions = {item["completion"]: item for item in row["completions"]}
        assert set(completions) == {"neutral", "multiplicative"}
        for item in completions.values():
            assert item["identity_at_zero_cuts"] is True
            assert item["positive"] is True
            assert item["disconnected_cut_composition"] is True
            assert item["refinement_regrouping_invariant"] is True
            assert all(
                check["commutes"] is True
                and check["coarse_value"] == check["iterated_refined_value"]
                for check in item["refinement_partition_checks"]
            )


def test_blocked_event_semantics_are_not_selected(payload: dict) -> None:
    verdict = payload["blocked_event_verdict"]
    assert verdict["one_class_selected"] is False
    assert verdict["six_class_total_selected"] is False
    assert verdict["no_capacity_action_selected"] is False
    assert all(
        row["blocked_event_factors_pairwise_distinct"] is True
        and row["blocked_event_actions_pairwise_distinct_at_one_cut"] is True
        for row in payload["exact_witnesses"]
    )
    for row in payload["exact_witnesses"]:
        completions = {
            item["completion"]: item
            for item in row["blocked_event_completions"]
        }
        assert set(completions) == {
            "no_capacity_action",
            "one_class_projector",
            "six_class_total_projector",
        }
        assert len(
            {item["values_m_0_to_3"][1] for item in completions.values()}
        ) == 3
        for item in completions.values():
            assert item["positive"] is True
            assert item["disconnected_cut_composition"] is True
            assert item["refinement_regrouping_invariant"] is True


def test_no_physical_or_comparison_promotion(payload: dict) -> None:
    assert all(value is False for value in payload["comparison_boundary"].values())
    assert payload["composition_boundary"]["per_edge_multiplication_authorized"] is False


@pytest.mark.parametrize(
    ("mutation", "code"),
    [
        (lambda p: p.__setitem__("status", "SOURCE_POSITIVE_ONE_FINITE_CUT_FACTOR"), "STATUS"),
        (
            lambda p: p["exact_witnesses"][0].__setitem__("different_at_one_cut", False),
            "ARITHMETIC",
        ),
        (
            lambda p: p["exact_witnesses"][0]["completions"][1].__setitem__(
                "disconnected_cut_composition", False
            ),
            "ARITHMETIC",
        ),
        (
            lambda p: p["exact_witnesses"][0]["completions"][1][
                "refinement_partition_checks"
            ][1].__setitem__("commutes", False),
            "ARITHMETIC",
        ),
        (
            lambda p: p["blocked_event_verdict"].__setitem__("one_class_selected", True),
            "SELECTION",
        ),
        (
            lambda p: p["exact_witnesses"][0]["blocked_event_completions"][2].__setitem__(
                "disconnected_cut_composition", False
            ),
            "ARITHMETIC",
        ),
        (
            lambda p: p["comparison_boundary"].__setitem__(
                "forecast_or_comparison_permitted", True
            ),
            "PROMOTION",
        ),
        (
            lambda p: p["producer_cone"].__setitem__(
                "cosmological_target_payload_used_in_countermodel", True
            ),
            "PROMOTION",
        ),
        (
            lambda p: p["composition_boundary"].__setitem__(
                "per_edge_multiplication_authorized", True
            ),
            "COMPOSITION",
        ),
        (
            lambda p: p["source_pin"].__setitem__("sha256", "sha256:" + "0" * 64),
            "SOURCE_PIN",
        ),
    ],
)
def test_mutations_fail_closed(payload: dict, mutation, code: str) -> None:
    changed = copy.deepcopy(payload)
    mutation(changed)
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.validate(changed)
    assert excinfo.value.code == code


def test_runtime_is_byte_exact(payload: dict) -> None:
    assert cert.OUTPUT_PATH.read_bytes() == cert.canonical_json_bytes(payload)
    cert.validate(payload)


@pytest.mark.parametrize(
    ("mutation", "code"),
    [
        (
            lambda source: source["reserve_trace_branch"].__setitem__(
                "binding_theorem_present", True
            ),
            "SOURCE_SEMANTICS",
        ),
        (
            lambda source: source["P_certified_enclosure"].__setitem__("hi", "5"),
            "DOMAIN",
        ),
    ],
)
def test_finite_source_boundary_mutations_fail_closed(
    tmp_path, monkeypatch: pytest.MonkeyPatch, mutation, code: str
) -> None:
    source = json.loads(cert.SOURCE_PROJECTION_PATH.read_text(encoding="utf-8"))
    mutation(source)
    changed = tmp_path / "finite_source.json"
    changed.write_text(json.dumps(source), encoding="utf-8")
    monkeypatch.setattr(cert, "SOURCE_PROJECTION_PATH", changed)
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.build()
    assert excinfo.value.code == code
