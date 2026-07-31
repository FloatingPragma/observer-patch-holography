from __future__ import annotations

import copy
import json

import pytest

import n_closure_branch_certificate as cert


@pytest.fixture(scope="module")
def payload() -> dict:
    return cert.build()


def test_two_retrospective_branches_are_explicitly_unselected(payload: dict) -> None:
    assert payload["status"] == cert.STATUS
    assert payload["scope"] == {
        "retrospective": True,
        "target_blind_forecast": False,
        "branch_selected": False,
        "comparison_data_consumed": False,
        "global_capacity_derived": False,
        "horizon_attachment_derived": False,
        "prediction_promoted": False,
    }
    branches = {row["branch_id"]: row for row in payload["branches"]}
    assert set(branches) == {cert.FINITE_BRANCH, cert.POISSON_BRANCH}
    assert all(row["selected"] is False for row in branches.values())


def test_source_typed_finite_branch_is_separate_from_poisson(payload: dict) -> None:
    branches = {row["branch_id"]: row for row in payload["branches"]}
    finite = branches[cert.FINITE_BRANCH]
    poisson = branches[cert.POISSON_BRANCH]
    assert finite["source_type"] == "finite_one_step_presence"
    assert finite["mean_count_or_projective_limit_carrier_required"] is False
    assert poisson["source_type"] == "poisson_or_projective_limit"
    assert poisson["mean_count_or_projective_limit_carrier_required"] is True
    assert poisson["additional_open_gate"] == {
        "mean_count_or_projective_limit_carrier": "open"
    }
    assert float(finite["factor"]) < float(poisson["factor"]) < 1


def test_all_physical_and_promotion_gates_fail_closed(payload: dict) -> None:
    assert all(value is False for value in payload["promotion_controls"].values())
    for row in payload["branches"]:
        assert all(value == "open" for value in row["open_gates"].values())
        assert row["open_gates"]["one_class_reserve_attachment"] == "open"
        assert row["open_gates"]["scalar_weighted_reserve_receipt"] == "open"
    assert payload["self_reference_boundary"]["same_typed_quantity_identified"] is False
    assert payload["self_reference_boundary"][
        "source_return_map_physically_attached"
    ] is False


@pytest.mark.parametrize(
    ("mutation", "code"),
    [
        (
            lambda value: value["branches"][0].__setitem__(
                "source_type", "poisson_or_projective_limit"
            ),
            "TYPE_CONFUSION",
        ),
        (
            lambda value: value["branches"][1].__setitem__(
                "mean_count_or_projective_limit_carrier_required", False
            ),
            "TYPE_CONFUSION",
        ),
        (
            lambda value: value["branches"][1].__setitem__(
                "additional_open_gate", {}
            ),
            "TYPE_CONFUSION",
        ),
        (
            lambda value: value["scope"].__setitem__("prediction_promoted", True),
            "PROMOTION",
        ),
        (
            lambda value: value["branches"][0].__setitem__("selected", True),
            "SELECTION",
        ),
        (
            lambda value: value["branches"][0].__setitem__("N", "1.0"),
            "ARITHMETIC",
        ),
        (
            lambda value: value["source_pins"][0].__setitem__(
                "sha256", "sha256:" + "0" * 64
            ),
            "SOURCE_PINS",
        ),
    ],
)
def test_scientific_mutations_fail_closed(
    payload: dict, mutation, code: str
) -> None:
    changed = copy.deepcopy(payload)
    mutation(changed)
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.validate(changed)
    assert excinfo.value.code == code


def test_source_status_mutation_fails_closed(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = json.loads(cert.P_SOURCE_PATH.read_text(encoding="utf-8"))
    source["status"] = "physical_endpoint_theorem"
    changed = tmp_path / "p_source.json"
    changed.write_text(json.dumps(source), encoding="utf-8")
    monkeypatch.setattr(cert, "P_SOURCE_PATH", changed)
    with pytest.raises(cert.CertificateError) as excinfo:
        cert.build()
    assert excinfo.value.code == "SOURCE_STATUS"


def test_manifest_is_byte_exact_and_validates(payload: dict) -> None:
    assert cert.OUTPUT_PATH.read_bytes() == cert.canonical_json_bytes(payload)
    cert.validate(payload)
