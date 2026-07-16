#!/usr/bin/env python3
"""Protocol tests for source-only emission and fail-closed scoring."""

from __future__ import annotations

import copy
import hashlib
import json
import sys
from decimal import Decimal, localcontext
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
TARGET_PATH = (
    REPO
    / "falsification"
    / "frozen_targets"
    / "hadronic_closure_target_2026-07-16_v3.json"
)
V2_PATH = (
    REPO
    / "falsification"
    / "frozen_targets"
    / "hadronic_closure_target_2026-07-14_v2.json"
)
HISTORICAL_V1_PATH = HERE / "runtime" / "ward_projected_payload_bracket_current.json"

if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import payload_harness as ph  # noqa: E402
import run_bracket  # noqa: E402
import score_bracket  # noqa: E402


@pytest.fixture(scope="module")
def target() -> dict:
    return json.loads(TARGET_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def source_artifact() -> dict:
    # The emitter is deliberately restricted to the source-derived P.  It is
    # not permitted to read the measurement-located P from the target.
    return run_bracket.build_bracket(ph.build_evaluation_point(), fast=True)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rehash(artifact: dict) -> None:
    artifact["content_sha256"] = score_bracket.artifact_content_sha256(artifact)


def _source_contract(source_artifact: dict) -> dict:
    return {
        "payload_object": source_artifact["payload_object"],
        "source_family_id": source_artifact["source_family_id"],
        "scheme_id": source_artifact["scheme"]["scheme_id"],
        "current_definition_id": source_artifact["current_definition_id"],
        "p_domain": copy.deepcopy(source_artifact["p_domain"]),
        "coordinate_schema": copy.deepcopy(source_artifact["coordinate_schema"]),
    }


def test_historical_records_are_hash_locked() -> None:
    assert (
        _sha256(V2_PATH)
        == "8dcefe93f124e21295b8e7cd85f2524ab17c65ad003c4fdac6205213d03ac6b2"
    )
    assert (
        _sha256(HISTORICAL_V1_PATH)
        == "395b1b6bc3d53ec40222e03004cced920236da3ed3b3e46cd407e71ab4852a15"
    )


def test_corrective_target_is_not_retroactively_activated(target: dict) -> None:
    assert target["registration_status"].endswith("not_scorable")
    assert target["frozen_utc"] is None
    assert target["promotion_or_falsification_allowed"] is False
    assert target["supersedes_for_future_evaluation"].endswith("2026-07-14_v2")


def test_corrective_target_has_independent_map_algebra(target: dict) -> None:
    score_bracket._validate_corrective_target(target)
    measurement = target["measurement_coordinate"]
    cl1 = target["map_targets"]["CL-1"]["point_diagnostics_only"]
    cl2 = target["map_targets"]["CL-2"]["point_diagnostics_only"]
    alpha_u = Decimal(measurement["alpha_U_at_P_target_point"])
    q_naive = Decimal(measurement["Delta_quark_naive_at_P_target_point"])

    with localcontext() as context:
        context.prec = 80
        assert (
            Decimal(cl1["Delta_source_total_target"])
            - Decimal(cl2["Delta_source_total_target"])
            == alpha_u
        )
        assert abs(
            Decimal(cl1["S_QEW_effective_target"])
            - Decimal(cl2["S_QEW_effective_target"])
            - alpha_u / q_naive
        ) < Decimal("1e-36")
    assert cl1["S_hadronic_target"] is None
    assert cl2["S_hadronic_target"] is None


def test_emitter_source_contains_no_target_or_tolerance_constants() -> None:
    source = (HERE / "run_bracket.py").read_text(encoding="utf-8")
    forbidden = (
        "CANON_NON_BLIND",
        "PASS_TOLERANCE",
        "S_required",
        "Delta_required",
        "Delta_missing",
        "pass_tolerance",
        "width_over_tolerance",
        "0.895400132647",
        "0.041465861005",
        "2.115394964782e-8",
        "4.286673538686e-9",
    )
    for token in forbidden:
        assert token not in source


def test_emitter_artifact_is_source_only_and_unscored(source_artifact: dict) -> None:
    serialized = json.dumps(source_artifact, sort_keys=True)
    target_constant = "0.895400132647658797805800284922791290"
    residual_constant = "0.041465861005223389053448723971273398"
    assert source_artifact["artifact"] == score_bracket.EMISSION_ARTIFACT
    assert source_artifact["schema_version"] == score_bracket.EMISSION_SCHEMA_VERSION
    assert set(source_artifact["coordinate_schema"]) == set(score_bracket.COORDINATES)
    assert "delta_source_total_alpha_inv" in source_artifact["bracket"]
    assert target_constant not in serialized
    assert residual_constant not in serialized
    assert "pass_tolerance" not in serialized.lower()
    assert source_artifact["target_or_measurement_inputs_used_in_computation"] is False
    assert source_artifact["promotion_allowed"] is False
    for coordinate, expected in score_bracket.REQUIRED_COORDINATE_TYPES.items():
        actual = source_artifact["coordinate_schema"][coordinate]
        assert (actual["kind"], actual["units"], actual["scoring_role"]) == expected
    assert source_artifact["p_domain"]["kind"] == "singleton_source_diagnostic"
    assert source_artifact["certification"]["status"].startswith("uncertified")
    assert (
        score_bracket.artifact_content_sha256(source_artifact)
        == source_artifact["content_sha256"]
    )


def test_measurement_located_p_cannot_enter_emitter(target: dict) -> None:
    target_p = target["measurement_coordinate"]["P_target_point_diagnostic"]
    ep = ph.build_evaluation_point(
        p_value=target_p,
        p_provenance="test_only_measurement_coordinate",
    )
    with pytest.raises(ValueError, match="measurement-located"):
        run_bracket.build_bracket(ep, fast=True)


def test_unactivated_target_fails_closed(source_artifact: dict, target: dict) -> None:
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket.score_artifact(source_artifact, target)
    assert excinfo.value.code == "target_not_activated"


def test_target_coordinate_tamper_fails_before_activation(
    source_artifact: dict, target: dict
) -> None:
    invalid_target = copy.deepcopy(target)
    invalid_target["map_targets"]["CL-2"]["point_diagnostics_only"][
        "Delta_source_total_target"
    ] = invalid_target["map_targets"]["CL-1"]["point_diagnostics_only"][
        "Delta_source_total_target"
    ]
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket.score_artifact(source_artifact, invalid_target)
    assert excinfo.value.code == "target_coordinate_mismatch"


@pytest.mark.parametrize("mutation", ["A", "joint_totals", "joint_s", "formulas"])
def test_corrective_target_rejects_algebra_preserving_looking_tamper(
    target: dict, mutation: str
) -> None:
    invalid = copy.deepcopy(target)
    if mutation == "A":
        invalid["measurement_coordinate"]["A_target_inverse_alpha"] = "999"
    elif mutation == "joint_totals":
        with localcontext() as context:
            context.prec = 100
            for map_name in ("CL-1", "CL-2"):
                point = invalid["map_targets"][map_name]["point_diagnostics_only"]
                point["Delta_source_total_target"] = str(
                    Decimal(point["Delta_source_total_target"]) + 1
                )
                point["Delta_source_residual_vs_implemented"] = str(
                    Decimal(point["Delta_source_residual_vs_implemented"]) + 1
                )
    elif mutation == "joint_s":
        with localcontext() as context:
            context.prec = 100
            for map_name in ("CL-1", "CL-2"):
                point = invalid["map_targets"][map_name]["point_diagnostics_only"]
                point["S_QEW_effective_target"] = str(
                    Decimal(point["S_QEW_effective_target"]) + Decimal("0.1")
                )
    else:
        invalid["map_targets"]["CL-1"]["map_formula"] = "garbage"
        invalid["map_targets"]["CL-2"]["map_formula"] = "garbage alpha_U"

    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket._validate_corrective_target(invalid)
    assert excinfo.value.code == "target_coordinate_mismatch"


def test_current_v3_cannot_be_activated_in_place(target: dict) -> None:
    invalid = copy.deepcopy(target)
    invalid["registration_status"] = "active"
    invalid["frozen_utc"] = "2026-07-16T00:00:00Z"
    invalid["promotion_or_falsification_allowed"] = True
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket._require_activated_target(invalid)
    assert excinfo.value.code == "target_not_activated"


def test_certified_interval_width_must_be_exact_decimal() -> None:
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket._validate_interval(
            {
                "lo": "1.0000000000000001",
                "hi": "1.0000000000000002",
                "width": "0",
            },
            "test.interval",
        )
    assert excinfo.value.code == "coordinate_schema_mismatch"


def test_all_zero_receipt_digest_is_rejected() -> None:
    with pytest.raises(score_bracket.ScoringError):
        score_bracket._verified_receipt(
            {"status": "PASS", "verified": True, "sha256": "0" * 64},
            "test.receipt",
        )


def test_false_same_subtraction_attestation_is_rejected(
    source_artifact: dict,
) -> None:
    invalid = copy.deepcopy(source_artifact)
    invalid["scheme"]["same_subtraction_as_a0"] = False
    _rehash(invalid)
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket._validate_artifact(invalid, _source_contract(invalid))
    assert excinfo.value.code == "artifact_schema_mismatch"


def test_scorer_rejects_p_domain_mismatch(source_artifact: dict) -> None:
    contract = _source_contract(source_artifact)
    contract["p_domain"]["hi"] = "1.7"
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket._validate_artifact(source_artifact, contract)
    assert excinfo.value.code == "evaluation_point_mismatch"


def test_current_grid_fails_certification_gate(source_artifact: dict) -> None:
    contract = _source_contract(source_artifact)
    contract["p_domain"].update({"kind": "registered_P_basin", "hi": "1.7"})
    uncertified = copy.deepcopy(source_artifact)
    uncertified["p_domain"] = copy.deepcopy(contract["p_domain"])
    _rehash(uncertified)
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket._validate_artifact(uncertified, contract)
    assert excinfo.value.code == "artifact_not_certified"


@pytest.mark.parametrize(
    ("coordinate", "wrong_kind"),
    [
        (score_bracket.TOTAL_COORDINATE, "residual"),
        (score_bracket.RESIDUAL_COORDINATE, "total"),
        (score_bracket.S_QEW_COORDINATE, "residual"),
        (score_bracket.S_HADRONIC_COORDINATE, "screening_ratio_qew"),
    ],
)
def test_scorer_rejects_coordinate_kind_mismatch(
    source_artifact: dict, coordinate: str, wrong_kind: str
) -> None:
    contract = _source_contract(source_artifact)
    contract["p_domain"]["kind"] = "registered_P_basin"
    contract["p_domain"]["hi"] = "1.7"
    invalid = copy.deepcopy(source_artifact)
    invalid["p_domain"] = copy.deepcopy(contract["p_domain"])
    invalid["coordinate_schema"][coordinate]["kind"] = wrong_kind
    _rehash(invalid)
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket._validate_artifact(invalid, contract)
    assert excinfo.value.code == "coordinate_schema_mismatch"


def test_hash_tamper_fails_artifact_validation(source_artifact: dict) -> None:
    tampered = copy.deepcopy(source_artifact)
    tampered["bracket"]["s_hadronic"]["lo"] += 0.01
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket._validate_artifact(tampered, {})
    assert excinfo.value.code == "artifact_hash_mismatch"


def test_schema_version_mismatch_fails_artifact_validation(
    source_artifact: dict,
) -> None:
    mismatched = copy.deepcopy(source_artifact)
    mismatched["schema_version"] -= 1
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket._validate_artifact(mismatched, {})
    assert excinfo.value.code == "artifact_schema_mismatch"


def test_historical_v1_is_not_a_v3_source_artifact() -> None:
    historical = json.loads(HISTORICAL_V1_PATH.read_text(encoding="utf-8"))
    with pytest.raises(score_bracket.ScoringError) as excinfo:
        score_bracket._validate_artifact(historical, {})
    assert excinfo.value.code == "artifact_schema_mismatch"
