from __future__ import annotations

import json
from decimal import Decimal, getcontext
from pathlib import Path


TARGET = (
    Path(__file__).resolve().parents[2]
    / "falsification"
    / "frozen_targets"
    / "hadronic_closure_target_2026-07-16_v3.json"
)


def _load() -> dict:
    return json.loads(TARGET.read_text(encoding="utf-8"))


def test_v3_fails_closed_until_a_new_external_registration() -> None:
    target = _load()

    assert target["registration_status"] == (
        "corrective_contract_not_yet_externally_timestamped_not_scorable"
    )
    assert target["frozen_utc"] is None
    assert target["promotion_or_falsification_allowed"] is False
    assert target["activation_requirements"]["external_timestamp_required"] is True
    assert target["activation_requirements"]["source_method_selection_status"] == (
        "post_target_access_protocol_scaffold_not_eligible"
    )
    assert "Do not populate this file in place" in target["activation_requirements"][
        "note"
    ]


def test_v3_separates_total_and_residual_coordinates() -> None:
    getcontext().prec = 80
    target = _load()
    measurement = target["measurement_coordinate"]
    maps = target["map_targets"]
    implemented = Decimal(measurement["Delta_source_implemented_at_P_target_point"])

    for row in ("CL-1", "CL-2"):
        point = maps[row]["point_diagnostics_only"]
        total = Decimal(point["Delta_source_total_target"])
        residual = Decimal(point["Delta_source_residual_vs_implemented"])
        assert total - implemented == residual
        assert "residual" not in maps[row]["pass_rule"].lower()

    type_rule = target["coordinate_definitions"]["type_rule"]
    assert "total payload" in type_rule
    assert "total target" in type_rule


def test_v3_has_distinct_map_targets_separated_by_alpha_u() -> None:
    getcontext().prec = 80
    target = _load()
    maps = target["map_targets"]
    alpha_u = Decimal(target["measurement_coordinate"]["alpha_U_at_P_target_point"])
    t1 = Decimal(maps["CL-1"]["point_diagnostics_only"]["Delta_source_total_target"])
    t2 = Decimal(maps["CL-2"]["point_diagnostics_only"]["Delta_source_total_target"])

    assert maps["CL-1"]["closure_rows"] == ["CL-1"]
    assert maps["CL-2"]["closure_rows"] == ["CL-2"]
    assert maps["CL-1"]["target_id"] != maps["CL-2"]["target_id"]
    assert t1 != t2
    assert t1 - t2 == alpha_u
    assert target["cross_map_invariants"]["one_scalar_closes_both_maps"] is False


def test_zero_ew_s_values_are_diagnostic_and_map_specific() -> None:
    getcontext().prec = 80
    target = _load()
    measurement = target["measurement_coordinate"]
    maps = target["map_targets"]
    alpha_u = Decimal(measurement["alpha_U_at_P_target_point"])
    q_naive = Decimal(measurement["Delta_quark_naive_at_P_target_point"])
    s1 = Decimal(maps["CL-1"]["point_diagnostics_only"]["S_QEW_effective_target"])
    s2 = Decimal(maps["CL-2"]["point_diagnostics_only"]["S_QEW_effective_target"])

    assert abs((s1 - s2) - alpha_u / q_naive) < Decimal("1e-36")
    assert maps["CL-1"]["point_diagnostics_only"]["S_hadronic_target"] is None
    assert maps["CL-2"]["point_diagnostics_only"]["S_hadronic_target"] is None
    assert target["shared_payload_contract"]["Delta_EW_zero_status"] == (
        "declared_zero_branch_unproven"
    )


def test_primary_contract_requires_a_function_and_independent_map_solves() -> None:
    target = _load()
    payload = target["shared_payload_contract"]
    scoring = target["primary_scoring_rule"]

    assert "function or certified interval enclosure" in payload["required_object"]
    assert payload["singleton_P_payload_eligible"] is False
    assert payload["sampled_grid_extrema_interval_certificate"] is False
    assert scoring["map_verdicts_independent"] is True
    assert scoring["missing_scheme_or_error_field_result"] == "NOT_EVALUABLE"
    assert scoring["open_Delta_EW_gate_result"] == "NOT_EVALUABLE"
    assert scoring["point_diagnostics_can_decide_verdict"] is False
    assert scoring["decision_policy_status"] == "not_yet_registered"
    assert len(scoring["required_future_decision_fields"]) == 5
    assert target["declared_outcomes"]["status"] == (
        "inactive_placeholders_for_a_future_successor"
    )
