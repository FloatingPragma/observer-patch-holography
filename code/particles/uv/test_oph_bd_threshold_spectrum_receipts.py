#!/usr/bin/env python3
"""Regression tests for the issue-368/369 fail-closed receipt bundle."""

from __future__ import annotations

import copy
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest
from jsonschema import Draft202012Validator


UV_DIR = Path(__file__).resolve().parent
REPO_ROOT = UV_DIR.parents[2]
PAPER_PATH = (
    REPO_ROOT
    / "extra/observer_patch_holography_as_string_vacuum_selector.tex"
)
if str(UV_DIR) not in sys.path:
    sys.path.insert(0, str(UV_DIR))

import build_oph_bd_threshold_spectrum_receipts as receipt_builder  # noqa: E402
from build_oph_bd_threshold_spectrum_receipts import (  # noqa: E402
    DEFAULT_OUT_DIR,
    DEFAULT_PACKET,
    PacketError,
    _canonical_bytes,
    _derive_one_higgs_geometry,
    _moduli_receipt_presence,
    _resolve_repo_path,
    _validate_uv_receipts,
    build_manifest,
    build_receipts,
)
from oph_bd_open_computation import compute_report  # noqa: E402


def _load(name: str) -> dict:
    return json.loads((DEFAULT_OUT_DIR / name).read_text(encoding="utf-8"))


@pytest.fixture(scope="module", autouse=True)
def _synchronize_bundle_after_upstream_artifact_tests() -> None:
    """Earlier particle tests rebuild hashed dependencies in the committed tree."""
    subprocess.run(
        [sys.executable, str(UV_DIR / "build_oph_bd_threshold_spectrum_receipts.py")],
        cwd=REPO_ROOT,
        check=True,
        text=True,
        capture_output=True,
    )


def test_committed_bundle_matches_in_memory_rebuild() -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(UV_DIR / "build_oph_bd_threshold_spectrum_receipts.py"),
            "--check",
        ],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    summary = json.loads(result.stdout)
    assert summary["mismatches"] == []
    assert summary["status"] == "OPEN_SOURCE_DATA_INSUFFICIENT"
    assert summary["compatibility_evaluated"] is False
    assert summary["promotion_allowed"] is False


def test_source_packet_defines_every_required_uv_input_and_leaves_it_missing() -> None:
    packet = json.loads(DEFAULT_PACKET.read_text(encoding="utf-8"))
    assert packet["schema_version"] == 3
    required = {
        "augmented_interval_isolation_receipt",
        "branch_global_coverage_receipt",
        "bundle_moduli_point_and_mass_matrix",
        "bulk_five_brane_sector_and_moduli_or_absence",
        "complex_structure_moduli_point_and_mass_matrix",
        "comparative_catalogue_coverage_receipt",
        "completed_constraint_system_and_stability_receipt",
        "conventional_soft_terms_boundary_conditions",
        "conventional_susy_breaking_and_mediation_data",
        "dilaton_point_and_mass",
        "heavy_nonzero_mode_spectrum",
        "hidden_sector_completion_and_spectrum",
        "kahler_moduli_point_metric_and_mass_matrix",
        "low_energy_spectrum_solver_receipt",
        "normalized_yukawa_boundary_data",
        "bd_equivalent_non_susy_uv_deformation_receipt",
        "physical_moduli_jacobian",
        "physical_quotient_descent_receipt",
        "string_compactification_and_matching_scales",
        "supersymmetric_sector_decoupling_map",
        "target_registry_precommitment_receipt",
        "threshold_matching_receipt",
    }
    assert set(packet["bd_uv_inputs"]) == required
    assert all(packet["bd_uv_inputs"][name] is None for name in required)
    assert packet["literature_inventory"]["available"]["kahler_moduli_count"] == 11
    assert packet["literature_inventory"]["available"]["complex_structure_moduli_count"] == 11
    assert packet["literature_inventory"]["available"]["bundle_moduli_count"] == 51
    assert packet["continuation_issues"] == [369]
    assert packet["comparison_registry"] == {
        "common_threshold_scheme_fixed": False,
        "complete_oph_target": False,
        "status": "five_coordinate_mixed_status_comparison_contract",
    }
    assert (
        packet["literature_inventory"]["available"][
            "one_higgs_locus_codimension_complex"
        ]
        == 2
    )


def test_one_higgs_geometry_is_derived_and_cross_checked() -> None:
    packet = json.loads(DEFAULT_PACKET.read_text(encoding="utf-8"))
    geometry = _derive_one_higgs_geometry(packet)
    assert geometry == {
        "bundle_tangent_complex": 49,
        "cokernel_complex": 2,
        "constituent_tangent_complex": 7,
        "domain_complex": 8,
        "codomain_complex": 9,
        "extension_tangent_complex": 42,
        "kernel_complex": 1,
        "normal_complex": 2,
        "rank": 7,
    }

    wrong_codimension = copy.deepcopy(packet)
    wrong_codimension["literature_inventory"]["available"][
        "one_higgs_locus_codimension_complex"
    ] = 3
    with pytest.raises(PacketError, match="determinantal normal dimension"):
        _derive_one_higgs_geometry(wrong_codimension)

    wrong_decomposition = copy.deepcopy(packet)
    wrong_decomposition["literature_inventory"]["available"][
        "one_higgs_determinantal_data"
    ]["extension_tangent_dimension_complex"] = 41
    with pytest.raises(PacketError, match="tangent decomposition"):
        _derive_one_higgs_geometry(wrong_decomposition)


def test_moduli_receipt_presence_is_derived_from_packet_slots() -> None:
    packet = json.loads(DEFAULT_PACKET.read_text(encoding="utf-8"))
    assert _moduli_receipt_presence(packet) == {
        "all_completion_receipt_slots_populated": False,
        "physical_jacobian_receipt_slot_populated": False,
        "selected_moduli_point_declared": False,
    }
    populated = copy.deepcopy(packet)
    populated["bd_branch"]["selected_moduli_point"] = "declared-point"
    populated["bd_uv_inputs"]["physical_moduli_jacobian"] = {
        "test_only": "populated"
    }
    presence = _moduli_receipt_presence(populated)
    assert presence["selected_moduli_point_declared"] is True
    assert presence["physical_jacobian_receipt_slot_populated"] is True

    route_complete = copy.deepcopy(packet)
    for name in receipt_builder.BASE_REQUIRED_UV_INPUTS:
        route_complete["bd_uv_inputs"][name] = {"test_only": "populated"}
    route_complete["bd_uv_inputs"]["supersymmetric_sector_decoupling_map"] = {
        "route": "conventional_susy_breaking"
    }
    route_complete["bd_uv_inputs"]["conventional_susy_breaking_and_mediation_data"] = {
        "test_only": "populated"
    }
    route_complete["bd_uv_inputs"]["conventional_soft_terms_boundary_conditions"] = {
        "test_only": "populated"
    }
    assert _moduli_receipt_presence(route_complete)[
        "all_completion_receipt_slots_populated"
    ] is True


def test_source_packet_schema_types_the_full_frozen_surface() -> None:
    packet = json.loads(DEFAULT_PACKET.read_text(encoding="utf-8"))
    schema = json.loads(
        (UV_DIR / "oph_bd_threshold_spectrum_inputs.schema.json").read_text(
            encoding="utf-8"
        )
    )
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)
    assert list(validator.iter_errors(packet)) == []

    malformed = copy.deepcopy(packet)
    malformed["external_sources"] = "not-a-source-list"
    malformed["target_coordinates"] = {}
    malformed["numerical_policy"] = "not-a-policy"
    errors = list(validator.iter_errors(malformed))
    assert len(errors) >= 3

    drifted_proxy = copy.deepcopy(packet)
    assumptions = drifted_proxy["proxy_assumptions"]
    assumptions["alpha1_normalization"] = "arbitrary"
    assumptions["gauge_loop_order"] = 2
    assumptions["m_susy_scan_gev"] = ["500"]
    assumptions["mssm_beta_coefficients"] = ["33/5"]
    assumptions["stop_mixing_scan"] = ["garbage"]
    drifted_proxy["numerical_policy"]["decimal_precision_digits"] = 60
    errors = list(validator.iter_errors(drifted_proxy))
    assert len(errors) >= 6

    malformed_coordinates = copy.deepcopy(packet)
    malformed_coordinates["target_coordinates"].pop("v_gev")
    malformed_coordinates["target_coordinates"]["mH_gev"]["value"] = "not-a-number"
    malformed_coordinates["reference_inputs"]["alpha3_reference_mz"].pop(
        "uncertainty"
    )
    errors = list(validator.iter_errors(malformed_coordinates))
    assert len(errors) >= 3


def test_loader_rejects_packet_selected_permissive_schema(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    packet = json.loads(DEFAULT_PACKET.read_text(encoding="utf-8"))
    canonical_schema = tmp_path / "canonical_schema.json"
    canonical_schema.write_text(
        (UV_DIR / "oph_bd_threshold_spectrum_inputs.schema.json").read_text(
            encoding="utf-8"
        ),
        encoding="utf-8",
    )
    permissive_schema = tmp_path / "permissive_schema.json"
    permissive_schema.write_text("{}\n", encoding="utf-8")
    packet["$schema"] = permissive_schema.name
    packet_path = tmp_path / "packet.json"
    packet_path.write_text(
        json.dumps(packet, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(receipt_builder, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(receipt_builder, "SCHEMA_PATH", canonical_schema)
    with pytest.raises(PacketError, match="must use the canonical issue-368 schema"):
        receipt_builder._load_packet(packet_path)


def test_proxy_reproduces_low_energy_arithmetic_without_using_bd_values() -> None:
    proxy = _load("proxy_targets.json")
    low = proxy["low_energy_tree_coordinates"]
    assert proxy["bd_branch_values_used"] is False
    assert proxy["status"] == "PROXY_TARGETS_REPRODUCED_NOT_BD_COMPATIBILITY"
    assert Decimal(low["y_t_tree_coordinate"]).quantize(Decimal("1e-12")) == Decimal(
        "0.987745211164"
    )
    assert Decimal(low["lambda_H_tree_coordinate"]).quantize(
        Decimal("1e-12")
    ) == Decimal("0.128706603202")
    assert Decimal(low["lambda_MSSM_tree_ceiling_coordinate"]).quantize(
        Decimal("1e-12")
    ) == Decimal("0.068973725409")
    assert Decimal(low["Delta_lambda_large_tan_beta_coordinate"]).quantize(
        Decimal("1e-12")
    ) == Decimal("0.059732877792")
    assert Decimal(low["m_h_tree_ceiling_GeV"]).quantize(
        Decimal("1e-9")
    ) == Decimal("91.652460286")
    assert len(proxy["stop_threshold_proxy"]) == 8
    assert Decimal(proxy["stop_threshold_proxy"][0]["required_M_S_GeV"]).quantize(
        Decimal("1e-3")
    ) == Decimal("3046.266")
    assert all(
        row["proxy_validity_warning"]
        == "no_BD_soft_terms_and_no_controlled_running_scheme"
        for row in proxy["stop_threshold_proxy"]
    )


def test_gauge_proxy_uses_current_sourced_reference_and_exposes_scale_sensitivity() -> None:
    gauge = _load("proxy_targets.json")["gauge_unification_proxy"]
    assert gauge["alpha3_reference_mz"] == "0.118"
    assert gauge["alpha3_reference_uncertainty"] == "0.0009"
    assert [row["M_SUSY_GeV"] for row in gauge["rows"]] == [
        "500",
        "1000",
        "2000",
        "5000",
        "10000",
    ]
    one_tev = next(row for row in gauge["rows"] if row["M_SUSY_GeV"] == "1000")
    assert Decimal(one_tev["alpha3_pred_mZ"]).quantize(
        Decimal("1e-12")
    ) == Decimal("0.111511308019")
    interval = one_tev["required_combined_inverse_coupling_correction"][
        "interval_from_alpha3_reference_only"
    ]
    central = Decimal(
        one_tev["required_combined_inverse_coupling_correction"]["central"]
    )
    assert Decimal(interval[0]) < central < Decimal(interval[1])


def test_threshold_and_moduli_certificates_fail_closed() -> None:
    threshold = _load("threshold_spectrum_certificate.json")
    moduli = _load("moduli_locking_certificate.json")
    assert threshold["status"] == "OPEN_SOURCE_DATA_INSUFFICIENT"
    assert threshold["compatibility_evaluated"] is False
    assert threshold["bd_threshold_spectrum_certificate_receipt"] is False
    assert threshold["promotion_allowed"] is False
    assert threshold["decoupling_route_policy"]["oph_incorporates_supersymmetry"] is False
    assert "supersymmetric_sector_decoupling_map" in threshold["missing_source_objects"]
    assert "dilaton_point_and_mass" in threshold["missing_source_objects"]
    assert "normalized_yukawa_boundary_data" in threshold["missing_source_objects"]
    assert (
        "bulk_five_brane_sector_and_moduli_or_absence"
        in threshold["missing_source_objects"]
    )
    assert (
        "conventional_soft_terms_boundary_conditions"
        not in threshold["missing_source_objects"]
    )
    assert threshold["acceptance_criteria"] == {
        "proxy_numbers_are_not_a_heavy_spectrum_certificate": True,
        "selected_string_witness_promoted_to_full_OPH_target": False,
        "threshold_compatibility_backed_by_reproducible_calculation": False,
        "threshold_gate_remains_explicitly_open": True,
    }
    assert moduli["status"] == "FAILED_NO_COMPLETED_RANK_ISOLATION_CERTIFICATE"
    assert moduli["certificate_issue"] == 369
    assert moduli["source_packet_issue"] == 368
    assert moduli["dynamical_stability_verified"] is False
    assert moduli["full_transverse_rank_verified"] is False
    assert moduli["constraint_augmented_rank_verified"] is False
    assert moduli["interval_existence_isolation_verified"] is False
    assert moduli["isolation_verified"] is False
    assert moduli["bd_moduli_locking_certificate_receipt"] is False
    assert moduli["rank_obstruction_certificate_receipt"] is True
    assert moduli["promotion_allowed"] is False
    for required_missing_object in (
        "augmented_interval_isolation_receipt",
        "branch_global_coverage_receipt",
        "comparative_catalogue_coverage_receipt",
        "completed_constraint_system_and_stability_receipt",
        "physical_quotient_descent_receipt",
        "target_registry_precommitment_receipt",
    ):
        assert required_missing_object in moduli["missing_source_objects"]
    assert (
        moduli["physical_source_slice"]["all_completion_receipt_slots_populated"]
        is False
    )
    assert (
        moduli["physical_source_slice"]["completed_constraint_locus_verified"]
        is False
    )
    assert moduli["physical_source_slice"]["completed_dimension_certified"] is False
    assert moduli["physical_source_slice"]["completed_dimension_real"] is None
    published = moduli["physical_source_slice"]["published_precompletion_slice"]
    assert published["dimension_complex"] == 71
    assert published["dimension_real"] == 142
    one_higgs = moduli["one_higgs_stratum_certificate"]
    assert one_higgs["ambient_normal_dimension_complex"] == 2
    assert one_higgs["pullback_transverse_rank_complex"] is None
    assert one_higgs["pullback_transverse_rank_verified"] is False
    assert moduli["rank_test"]["emitted_proxy_jacobian_rank"] == 0
    assert (
        moduli["rank_test"]["physical_jacobian_receipt_slot_populated"]
        is False
    )
    assert moduli["rank_test"]["augmented_constraint_map_evaluable"] is False
    assert moduli["rank_test"]["selected_moduli_point_declared"] is False
    assert (
        moduli["rank_test"]["published_slice_five_coordinate_rank_upper_bound"]
        == 5
    )
    assert (
        moduli["rank_test"][
            "published_slice_five_coordinate_infinitesimal_nullity_lower_bound"
        ]
        == 137
    )
    assert (
        moduli["rank_test"][
            "necessary_net_real_completion_cuts_for_promoted_rows"
        ]
        == 139
    )
    assert (
        moduli["verdict"]["operator_safe_selected_candidate"]
        == "RETRACTED_BY_MODULI_LOCKING_GATE"
    )
    assert moduli["verdict"]["recovered_oph_core"] == "UNAFFECTED"
    assert moduli["target_map"]["common_threshold_scheme_fixed"] is False
    assert moduli["target_map"]["complete_oph_target"] is False
    assert moduli["target_map"]["promoted_coordinate_count"] == 3


def test_manifest_hashes_every_receipt_and_dependency() -> None:
    receipts = build_receipts(DEFAULT_PACKET)
    expected_manifest = build_manifest(receipts, DEFAULT_PACKET)
    committed_manifest = _load("manifest.json")
    assert committed_manifest == expected_manifest
    for filename, payload in receipts.items():
        expected_bytes = _canonical_bytes(payload)
        entry = committed_manifest["artifacts"][filename]
        assert entry["byte_count"] == len(expected_bytes)
        assert entry["sha256"] == hashlib.sha256(expected_bytes).hexdigest()
    source_dag = _load("source_dag.json")
    assert all(
        len(entry["sha256"]) == 64
        for entry in source_dag["local_dependencies"]
    )
    assert all(
        len(entry["sha256"]) == 64
        for entry in source_dag["external_sources"]
    )
    assert all(
        entry["verification_scope"]
        == "recorded_hash; this builder does not fetch external bytes"
        for entry in source_dag["external_sources"]
    )


def test_legacy_open_computation_contract_is_preserved() -> None:
    report = compute_report()
    assert set(report) == {
        "gauge_unification_proxy",
        "inputs",
        "low_energy_targets",
        "open_gates",
        "operator_audit",
        "selected_representative",
        "status",
        "stop_threshold_proxy",
    }
    assert report["inputs"]["alpha3_reference_mz"] == 0.118
    assert report["status"] == "OPEN_SOURCE_DATA_INSUFFICIENT"
    assert "dilaton_point_and_mass" in report["open_gates"]
    committed_legacy = json.loads(
        (
            REPO_ROOT
            / "code/particles/runs/uv/oph_bd_open_computation_results.json"
        ).read_text(encoding="utf-8")
    )
    assert committed_legacy == report


def test_source_paths_reject_absolute_parent_and_pdf_dependencies() -> None:
    with pytest.raises(PacketError):
        _resolve_repo_path("/tmp/injected.json")
    with pytest.raises(PacketError):
        _resolve_repo_path("../injected.json")
    with pytest.raises(PacketError):
        _resolve_repo_path("extra/observer_patch_holography_as_string_vacuum_selector.pdf")


def test_non_null_uv_receipts_are_bound_to_issue_branch_and_slot(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(receipt_builder, "REPO_ROOT", tmp_path)
    packet = json.loads(DEFAULT_PACKET.read_text(encoding="utf-8"))
    receipt_path = "heavy_spectrum_receipt.json"
    envelope = {
        "artifact": "oph_bd_uv_input_receipt_envelope",
        "bd_branch_identifier": packet["bd_branch"]["identifier"],
        "issue": 368,
        "schema_version": 1,
        "uv_input": "heavy_nonzero_mode_spectrum",
    }
    receipt_file = tmp_path / receipt_path
    receipt_file.write_text(
        json.dumps(envelope, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    receipt_bytes = receipt_file.read_bytes()
    valid_receipt = {
        "path": receipt_path,
        "sha256": hashlib.sha256(receipt_bytes).hexdigest(),
        "status": "hash_and_envelope_verified",
    }

    valid_packet = copy.deepcopy(packet)
    valid_packet["bd_uv_inputs"]["heavy_nonzero_mode_spectrum"] = valid_receipt
    ledger = _validate_uv_receipts(valid_packet)
    assert ledger[0]["uv_input"] == "heavy_nonzero_mode_spectrum"
    assert ledger[0]["sha256"] == valid_receipt["sha256"]

    missing_path_packet = copy.deepcopy(packet)
    missing_path_packet["bd_uv_inputs"]["heavy_nonzero_mode_spectrum"] = {
        **valid_receipt,
        "path": "does_not_exist.json",
    }
    with pytest.raises(PacketError, match="does not exist"):
        _validate_uv_receipts(missing_path_packet)

    bad_hash_packet = copy.deepcopy(packet)
    bad_hash_packet["bd_uv_inputs"]["heavy_nonzero_mode_spectrum"] = {
        **valid_receipt,
        "sha256": "0" * 64,
    }
    with pytest.raises(PacketError, match="hash mismatch"):
        _validate_uv_receipts(bad_hash_packet)

    failed_status_packet = copy.deepcopy(packet)
    failed_status_packet["bd_uv_inputs"]["heavy_nonzero_mode_spectrum"] = {
        **valid_receipt,
        "status": "FAILED",
    }
    with pytest.raises(PacketError, match="status=hash_and_envelope_verified"):
        _validate_uv_receipts(failed_status_packet)

    wrong_slot_file = tmp_path / "wrong_slot.json"
    wrong_slot_file.write_text(
        json.dumps(
            {**envelope, "uv_input": "normalized_yukawa_boundary_data"},
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    wrong_slot_packet = copy.deepcopy(packet)
    wrong_slot_packet["bd_uv_inputs"]["heavy_nonzero_mode_spectrum"] = {
        "path": wrong_slot_file.name,
        "sha256": hashlib.sha256(wrong_slot_file.read_bytes()).hexdigest(),
        "status": "hash_and_envelope_verified",
    }
    with pytest.raises(PacketError, match="wrong envelope field uv_input"):
        _validate_uv_receipts(wrong_slot_packet)


def test_paper_tracks_the_frozen_proxy_and_open_gate() -> None:
    paper = PAPER_PATH.read_text(encoding="utf-8")
    assert r"\mathcal C_3^{\rm proxy}=-0.493124115142" in paper
    assert r"\alpha_s(m_Z)=0.1180\pm0.0009" in paper
    assert r"\Delta\lambda_{\rm proxy}=0.059732877792" in paper
    assert r"m_{h,\rm tree}^{\rm max}=91.652460286" in paper
    assert r"3046.266\,\mathrm{GeV}" in paper
    assert "The OPH low-energy target used here has no supersymmetric partner sector." in paper
    assert "The threshold/spectrum certificate is not emitted here" in paper
    assert "0.521754255407" not in paper
    assert "0.118400" not in paper
