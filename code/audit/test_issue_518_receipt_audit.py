"""Adversarial tests for the issue #518 receipt promotion firewall."""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
AUDITOR_PATH = ROOT / "tools" / "audit_issue_518_receipts.py"
REGISTRY_PATH = Path(__file__).with_name("issue_518_receipt_registry.json")

spec = importlib.util.spec_from_file_location("audit_issue_518_receipts", AUDITOR_PATH)
auditor = importlib.util.module_from_spec(spec)
sys.modules["audit_issue_518_receipts"] = auditor
spec.loader.exec_module(auditor)


def _registry() -> dict:
    return json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))


def _artifact(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def _row(registry: dict, receipt_id: str) -> dict:
    return next(
        row for row in registry["rows"]
        if row["receipt_id"] == receipt_id
    )


def _audit(registry: dict, overrides: dict[str, dict] | None = None) -> dict:
    return auditor.audit_registry(registry, root=ROOT, overrides=overrides)


def test_live_registry_passes_and_promotes_only_literal_recomputations() -> None:
    report = _audit(_registry())

    assert report["pass"] is True
    assert report["class_source"] == "independent_registry_only"
    assert report["subject_declared_classes_ignored"] is True
    assert report["promoted_receipts"] == [
        "hierarchy-ru-independent-interval",
        "petz-classical-recovery",
        "hierarchy-capacity-antecedent-only-nogo",
        "higgs-naturality-antecedent-only-nogo",
    ]


def test_registry_uses_complete_controlled_vocabulary() -> None:
    registry = _registry()

    assert set(registry["receipt_classes"]) == {
        "identity",
        "schema",
        "producer",
        "recomputation",
        "prediction",
        "comparison",
    }
    assert registry["promotion_policy"]["never_sufficient_classes"] == [
        "identity",
        "schema",
        "comparison",
    ]


def test_subject_cannot_select_its_own_receipt_class() -> None:
    registry = _registry()
    path = "code/consensus/runs/verified_tree_packet_net_domain.json"
    artifact = _artifact(path)
    artifact["petz_domain"]["receipt_class"] = "comparison"

    report = _audit(registry, {path: artifact})

    assert report["pass"] is True
    petz = next(
        row for row in report["rows"]
        if row["receipt_id"] == "petz-classical-recovery"
    )
    assert petz["receipt_class"] == "recomputation"
    assert petz["promoted"] is True


def test_unknown_receipt_class_fails_closed() -> None:
    registry = _registry()
    _row(registry, "petz-classical-recovery")["receipt_class"] = "theorem_grade"

    report = _audit(registry)

    assert report["pass"] is False
    assert any("unknown_receipt_class" in failure for failure in report["failures"])


def test_identity_receipt_cannot_be_promoted() -> None:
    registry = _registry()
    row = _row(registry, "hierarchy-ew-capacity-identity")
    row["promoted"] = True
    row["open_gates"] = []

    report = _audit(registry)

    assert report["pass"] is False
    assert any(
        "receipt_class_not_promotion_eligible" in failure
        for failure in report["failures"]
    )


def test_backsolved_identity_cannot_reenter_hierarchy_closed_bundle() -> None:
    registry = _registry()
    path = "code/particles/hierarchy/manifest.json"
    artifact = _artifact(path)
    artifact["claim_boundary"]["closed_by_bundle"].append(
        "bridge-defined fixed-point identity N_CRC^EW(P)=target"
    )

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert "backsolved_capacity_identity_present_in_closed_bundle" in report[
        "global_failures"
    ]


def test_nogo_results_cannot_disappear_from_hierarchy_closed_bundle() -> None:
    registry = _registry()
    path = "code/particles/hierarchy/manifest.json"
    artifact = _artifact(path)
    artifact["claim_boundary"]["closed_by_bundle"] = [
        item
        for item in artifact["claim_boundary"]["closed_by_bundle"]
        if "non-identifiability theorem" not in item
    ]

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert "capacity_nonidentifiability_missing_from_closed_bundle" in report[
        "global_failures"
    ]
    assert "naturality_nonidentifiability_missing_from_closed_bundle" in report[
        "global_failures"
    ]


def test_promoted_row_cannot_consume_target_artifact() -> None:
    registry = _registry()
    row = _row(registry, "petz-classical-recovery")
    row["antecedent_artifacts"].append("fixture/target.json")
    row["target_artifacts"].append("fixture/target.json")

    report = _audit(registry)

    assert report["pass"] is False
    assert any(
        "promoted_row_consumes_target_artifact" in failure
        for failure in report["failures"]
    )


def test_petz_broken_column_normalization_fails_literal_recomputation() -> None:
    registry = _registry()
    path = "code/consensus/runs/verified_tree_packet_net_domain.json"
    artifact = _artifact(path)
    artifact["petz_domain"]["recovery_matrix"][0][0] = "3/5"

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any("petz_literal_mismatch" in failure for failure in report["failures"])


def test_petz_negative_choi_entry_fails_literal_recomputation() -> None:
    registry = _registry()
    path = "code/consensus/runs/verified_tree_packet_net_domain.json"
    artifact = _artifact(path)
    artifact["petz_domain"]["recovery_matrix"][1][0] = "-1/4"

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any("petz_literal_mismatch" in failure for failure in report["failures"])


def test_petz_zero_reference_support_fails_literal_recomputation() -> None:
    registry = _registry()
    path = "code/consensus/runs/verified_tree_packet_net_domain.json"
    artifact = _artifact(path)
    artifact["petz_domain"]["reference_sigma_B"] = {
        "0": "0",
        "1": "1/2",
        "2": "1/2",
    }

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any("petz_literal_mismatch:full_support" in failure for failure in report["failures"])


def test_hierarchy_defining_antecedent_perturbation_fails_closed() -> None:
    registry = _registry()
    path = (
        "code/particles/hierarchy/certificates/"
        "independent_interval_witness_receipt.json"
    )
    artifact = _artifact(path)
    artifact["phi_at_lower"] = "[-0.1, -0.1]"

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "hierarchy_lower_endpoint_not_strictly_positive" in failure
        for failure in report["failures"]
    )


def test_hierarchy_mutation_control_boolean_is_not_trusted_if_false() -> None:
    registry = _registry()
    path = (
        "code/particles/hierarchy/certificates/"
        "independent_interval_witness_receipt.json"
    )
    artifact = _artifact(path)
    artifact["perturbation_controls"]["pixel_perturbed_bracket_broken"] = False

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "hierarchy_mutation_controls_not_all_fail_closed" in failure
        for failure in report["failures"]
    )


def test_capacity_nogo_requires_same_antecedent_fingerprint() -> None:
    registry = _registry()
    path = (
        "code/particles/hierarchy/certificates/"
        "antecedent_only_nonidentifiability_receipt.json"
    )
    artifact = _artifact(path)
    artifact["capacity_nonidentifiability"]["countermodels"][1][
        "antecedent_fingerprint"
    ] = "0" * 64

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "capacity_nogo_antecedent_mismatch" in failure
        for failure in report["failures"]
    )


def test_capacity_nogo_requires_distinct_fixed_points() -> None:
    registry = _registry()
    path = (
        "code/particles/hierarchy/certificates/"
        "antecedent_only_nonidentifiability_receipt.json"
    )
    artifact = _artifact(path)
    artifact["capacity_nonidentifiability"]["countermodels"][1][
        "fixed_log_capacity"
    ] = "1"

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "capacity_nogo_outputs_not_distinct" in failure
        for failure in report["failures"]
    )


def test_capacity_nogo_recomputes_contraction_instead_of_trusting_boolean() -> None:
    registry = _registry()
    path = (
        "code/particles/hierarchy/certificates/"
        "antecedent_only_nonidentifiability_receipt.json"
    )
    artifact = _artifact(path)
    artifact["capacity_nonidentifiability"]["countermodels"][0][
        "strict_contraction"
    ] = False

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "capacity_nogo_literal_mismatch" in failure
        for failure in report["failures"]
    )


def test_antecedent_only_nogo_rejects_target_consumption() -> None:
    registry = _registry()
    path = (
        "code/particles/hierarchy/certificates/"
        "antecedent_only_nonidentifiability_receipt.json"
    )
    artifact = _artifact(path)
    artifact["target_artifacts_consumed"] = [
        "R_EW_global_capacity_certificate.json"
    ]

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "nogo_consumes_target_artifact" in failure
        for failure in report["failures"]
    )


def test_naturality_nogo_recomputes_defects_from_finite_maps() -> None:
    registry = _registry()
    path = (
        "code/particles/hierarchy/certificates/"
        "antecedent_only_nonidentifiability_receipt.json"
    )
    artifact = _artifact(path)
    artifact["naturality_nonidentifiability"]["countermodels"][1][
        "evaluated_defects"
    ] = {"epsilon_n": 0, "epsilon_h": 0, "epsilon_H": 0}

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "naturality_nogo_defect_recomputation_failed" in failure
        for failure in report["failures"]
    )


def test_naturality_nogo_requires_commuting_and_noncommuting_completions() -> None:
    registry = _registry()
    path = (
        "code/particles/hierarchy/certificates/"
        "antecedent_only_nonidentifiability_receipt.json"
    )
    artifact = _artifact(path)
    block = artifact["naturality_nonidentifiability"]
    commuting_maps = copy.deepcopy(block["countermodels"][0]["maps"])
    for model in block["countermodels"]:
        model["maps"] = copy.deepcopy(commuting_maps)
        model["evaluated_defects"] = {
            "epsilon_n": 0,
            "epsilon_h": 0,
            "epsilon_H": 0,
        }
    block["evaluated_epsilon_H_values"] = [0, 0, 0]

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "naturality_nogo_does_not_exhibit_distinct_defects" in failure
        for failure in report["failures"]
    )
    assert any(
        "naturality_nogo_normal_square_countermodel_missing" in failure
        for failure in report["failures"]
    )
    assert any(
        "naturality_nogo_obstruction_square_countermodel_missing" in failure
        for failure in report["failures"]
    )


def test_nogo_source_packet_mutation_breaks_every_antecedent_binding() -> None:
    registry = _registry()
    path = (
        "code/particles/hierarchy/certificates/"
        "antecedent_only_nonidentifiability_receipt.json"
    )
    artifact = _artifact(path)
    artifact["source_packet"]["P"] = "2"

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert "capacity_nogo_source_packet_hash_mismatch" in " ".join(
        report["failures"]
    )
    assert "naturality_nogo_source_packet_hash_mismatch" in " ".join(
        report["failures"]
    )


def test_nogo_source_schema_rejects_hidden_target_even_after_rehash() -> None:
    registry = _registry()
    path = (
        "code/particles/hierarchy/certificates/"
        "antecedent_only_nonidentifiability_receipt.json"
    )
    artifact = _artifact(path)
    artifact["source_packet"]["capacity_target"] = "281.0326"
    forged_hash = auditor._canonical_hash(artifact["source_packet"])
    artifact["source_packet_sha256"] = forged_hash
    for block_name in (
        "capacity_nonidentifiability",
        "naturality_nonidentifiability",
    ):
        for model in artifact[block_name]["countermodels"]:
            model["antecedent_fingerprint"] = forged_hash

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "nogo_source_packet_schema_mismatch" in failure
        for failure in report["failures"]
    )


def test_thomson_theorem_grade_boolean_cannot_change_verdict() -> None:
    registry = _registry()
    path = "code/P_derivation/runtime/thomson_endpoint_package_current.json"
    artifact = _artifact(path)
    artifact["issue_223_acceptance"]["theorem_grade_object_defined"] = False

    report = _audit(registry, {path: artifact})

    assert report["pass"] is True
    thomson = next(
        row for row in report["rows"]
        if row["receipt_id"] == "thomson-endpoint-blocker"
    )
    assert thomson["receipt_class"] == "comparison"
    assert thomson["promoted"] is False


def test_thomson_component_perturbation_fails_literal_recomputation() -> None:
    registry = _registry()
    path = "code/P_derivation/runtime/thomson_endpoint_package_current.json"
    artifact = _artifact(path)
    exact = artifact["codata_mapped_endpoint_packet"]["exact_one_loop_package"]
    exact["quark_delta_alpha_inv_screened"] = "4.5"

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "thomson_transport_components_do_not_sum" in failure
        for failure in report["failures"]
    )


def test_thomson_subject_promotion_boolean_fails_closed() -> None:
    registry = _registry()
    path = "code/P_derivation/runtime/thomson_endpoint_package_current.json"
    artifact = _artifact(path)
    artifact["promotion_allowed"] = True

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "thomson_comparison_improperly_promoted" in failure
        for failure in report["failures"]
    )


def test_planck_input_output_overlap_fails_closed() -> None:
    registry = _registry()
    row = _row(registry, "planck-posterior-capacity-comparison")
    row["planck_output_pointer"] = row["planck_input_pointer"]

    report = _audit(registry)

    assert report["pass"] is False
    assert any(
        "planck_input_output_pointers_overlap" in failure
        for failure in report["failures"]
    )


def test_planck_derived_gap_perturbation_fails_recomputation() -> None:
    registry = _registry()
    path = (
        "code/capacity_readback/planck_posterior/"
        "planck_lambda_to_N_propagation.json"
    )
    artifact = _artifact(path)
    artifact["combos"]["TTTEEE_lowE_lensing"]["gap_ratio"] = 1.0

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any(
        "planck_gap_ratio_recomputation_failed" in failure
        for failure in report["failures"]
    )


def test_dark_sector_calibration_fixture_cannot_become_active() -> None:
    registry = _registry()
    path = "code/audit/fixtures/issue_518/dark_sector_calibration_negative.json"
    artifact = _artifact(path)
    artifact["active_theorem_bundle"] = True

    report = _audit(registry, {path: artifact})

    assert report["pass"] is False
    assert any("dark_fixture_marked_active" in failure for failure in report["failures"])


def test_a5_layers_have_four_separate_receipt_classes() -> None:
    registry = _registry()
    expected = {
        "coefficient_algebra": "recomputation",
        "physical_currents": "producer",
        "global_descent": "identity",
        "matter_realization": "schema",
    }

    assert registry["a5_role_classes"] == expected
    assert len(set(expected.values())) == 4
    report = _audit(registry)
    assert report["pass"] is True


def test_a5_role_class_alias_fails_closed() -> None:
    registry = _registry()
    _row(registry, "a5-global-descent")["receipt_class"] = "producer"

    report = _audit(registry)

    assert report["pass"] is False
    assert any(
        "a5_role_receipt_class_mismatch" in failure
        for failure in report["failures"]
    )


def test_cli_exits_zero_for_live_registry() -> None:
    result = subprocess.run(
        [sys.executable, str(AUDITOR_PATH)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(result.stdout)["pass"] is True
