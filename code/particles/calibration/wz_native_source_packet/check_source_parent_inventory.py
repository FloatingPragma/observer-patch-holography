#!/usr/bin/env python3
"""Independent resolver for the issue-#594 source-parent frontier.

This checker imports no producer code. It resolves every pinned file, derives
the content digests and source-DAG result, verifies the finite parents with
their native certificate programs, and keeps every physical promotion gate
false while its named interface is open.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path, PurePosixPath
from typing import Any, Mapping

from jsonschema import Draft202012Validator


HERE = Path(__file__).resolve().parent
DEFAULT_REPO_ROOT = HERE.parents[3]
DEFAULT_INVENTORY = HERE / "outputs" / "source_parent_inventory.json"
SCHEMA_PATH = HERE / "schemas" / "source_parent_inventory_v1.schema.json"
POLICY_REL = (
    "code/particles/calibration/wz_native_source_packet/"
    "data/source_parent_policy_v1.json"
)
EXPECTED_POLICY_CANONICAL_SHA256 = (
    "4729fdc68168f31482e24dcf94f4db5e61a4702d053bb88cfc937eae40f4d12f"
)

EXPECTED_POSITIVE = {
    "screen_carrier": {
        "issue": 565,
        "verifier_id": "issue_565_carrier",
        "files": {
            "code/a5_closure/manifests/echosahedral_federation_reference.json":
                ("manifest", "oph.echosahedral_selector_manifest.v1"),
            "code/a5_closure/receipts/echosahedral_federation_reference.receipt.json":
                ("receipt", "oph.echosahedral_selector_receipt.v1"),
        },
    },
    "finite_port_current_algebra": {
        "issue": 566,
        "verifier_id": "issue_566_currents",
        "files": {
            "code/a5_closure/manifests/port_current_response_reference.json":
                ("manifest", "oph.port_current_response_manifest.v5"),
            "code/a5_closure/receipts/port_current_inner_reference.receipt.json":
                ("receipt", "oph.port_current_inner_receipt.v5"),
            "code/a5_closure/manifests/charged_response_semantic_artifact.json":
                ("semantic_artifact", "oph.charged_response_semantic_artifact.v3"),
        },
    },
    "finite_matter_module": {
        "issue": 314,
        "verifier_id": "issue_314_matter",
        "files": {
            "code/a5_closure/manifests/super_tannakian_matter_reference.json":
                ("manifest", "oph.super_tannakian_matter_manifest.v5"),
            "code/a5_closure/receipts/super_tannakian_matter_reference.receipt.json":
                ("receipt", "oph.super_tannakian_matter_receipt.v5"),
            "code/a5_closure/manifests/spin_statistics_semantic_artifact.json":
                ("semantic_artifact", "oph.spin_statistics_semantic_artifact.v1"),
        },
    },
    "finite_global_form": {
        "issue": 567,
        "verifier_id": "issue_567_global_form",
        "files": {
            "code/a5_closure/manifests/axis_center_descent_reference.json":
                ("manifest", "oph.axis_center_descent_manifest.v4"),
            "code/a5_closure/receipts/axis_center_descent_reference.receipt.json":
                ("receipt", "oph.axis_center_descent_receipt.v4"),
            "code/a5_closure/manifests/global_form_semantic_artifact.json":
                ("semantic_artifact", "oph.global_form_semantic_artifact.v1"),
        },
    },
}

EXPECTED_CONDITIONAL = {
    "family_band_candidate": {
        "issue": 569,
        "verifier_id": "issue_569_family_band",
        "files": {
            "code/a5_closure/manifests/family_band_attachment_reference.json":
                ("candidate_certificate", "oph.family_band_attachment_certificate.v6"),
            "code/a5_closure/manifests/charged_response_pole_residue_artifact.json":
                ("candidate_artifact", "oph.charged_response_pole_residue.v2"),
            "code/a5_closure/manifests/matter_attachment_receipt.json":
                (
                    "finite_domain_attachment_receipt",
                    "oph.local-domain-matter-attachment.v1",
                ),
        },
    },
    "matter_completeness_boundary": {
        "issue": 609,
        "verifier_id": "issue_609_matter_boundary",
        "files": {
            "code/a5_closure/manifests/matter_menu_spectral_ledger_reference.json":
                ("boundary_certificate", "oph.matter_menu_spectral_ledger_certificate.v1"),
        },
    },
    "finite_flux_and_local_domain_classical_context": {
        "issue": 311,
        "verifier_id": "issue_311_finite_spectral_classical_control",
        "files": {
            "code/a5_closure/manifests/flux_defect_criterion_reference.json":
                (
                    "boundary_certificate",
                    "oph.flux_defect_criterion_certificate.v3",
                ),
            "code/a5_closure/receipts/flux_defect_criterion_reference.receipt.json":
                ("receipt", "oph.flux_defect_criterion_receipt.v3"),
        },
    },
    "rg_representation_frontier": {
        "issue": 32,
        "verifier_id": "issue_32_rg_frontier",
        "files": {
            "code/P_derivation/source_rg_frontier/outputs/rg_representation_frontier.json":
                ("partial_source_frontier", "oph.rg_representation_frontier.v1"),
        },
    },
    "scalar_yukawa_source_frontier": {
        "issue": 630,
        "verifier_id": "issue_630_scalar_yukawa_frontier",
        "files": {
            "code/particles/hierarchy/higgs_yukawa_source_frontier/outputs/"
            "higgs_yukawa_source_frontier.json":
                ("partial_source_frontier", "oph.higgs_yukawa_source_frontier.v1"),
        },
    },
    "local_ew_order_unit_frontier": {
        "issue": 631,
        "verifier_id": "issue_631_local_carrier_frontier",
        "files": {
            "code/a5_closure/manifests/common_ew_order_unit_carrier_reference.json":
                ("candidate_certificate", "oph.common_ew_order_unit_carrier_frontier.v1"),
            "code/a5_closure/receipts/"
            "common_ew_order_unit_carrier_reference.receipt.json":
                ("receipt", "oph.common_ew_order_unit_carrier_receipt.v1"),
        },
    },
    "source_clock_frontier": {
        "issue": 633,
        "verifier_id": "issue_633_source_clock_frontier",
        "files": {
            "code/particles/hierarchy/source_clock_frontier/outputs/"
            "source_clock_frontier.json":
                ("partial_source_frontier", "oph.source_clock.frontier.v1"),
        },
    },
}

EXPECTED_RESOLVED_BOUNDARIES = {
    "finite_local_domain_boundary": {
        "issue": 634,
        "outcome": "closed_bounded",
        "verdict": "ATTAINED",
        "sentinel": "LOCAL_DOMAIN_INHABITATION_RECEIPT",
        "path": "code/a5_closure/manifests/stage4_receipt.json",
        "schema": "oph.local-domain-stage4.v1",
    },
    "finite_spectral_quantum_discrimination_boundary": {
        "issue": 311,
        "outcome": "closed_no_go",
        "verdict": "CLASSICAL_REALIZATION_MATCHES_DECLARED_FINITE_SPECTRAL_INTERFACE",
        "sentinel": "CLASSICAL_REALIZATION_RECEIPT",
        "path": "code/a5_closure/manifests/classical_realization_receipt.json",
        "schema": "oph.local-domain-classical-realization.v1",
    },
    "physical_unit_boundary": {
        "issue": 633,
        "outcome": "closed_not_evaluable",
        "verdict": "PHYSICAL_UNITS_NOT_EVALUABLE_ON_DECLARED_SERIALIZED_INTERFACE",
        "sentinel": "CLOCK_UNIT_BOUNDED_INTERFACE_AUDIT",
        "path": "code/a5_closure/manifests/clock_unit_verdict.json",
        "schema": "oph.local-domain-clock-unit-verdict.v1",
    },
}

EXPECTED_CONSUMER_SCHEMAS = {
    "canonical_action":
        "code/particles/calibration/wz_upstream_completion/schemas/"
        "sm_eft_action_packet_v1.schema.json",
    "full_yukawa":
        "code/particles/calibration/wz_upstream_completion/schemas/"
        "full_yukawa_packet_v1.schema.json",
    "eft_matching":
        "code/particles/calibration/wz_upstream_completion/schemas/"
        "eft_matching_packet_v1.schema.json",
    "source_law_covariance":
        "code/particles/calibration/wz_upstream_completion/schemas/"
        "source_law_covariance_packet_v1.schema.json",
    "operational_clock":
        "code/particles/calibration/wz_upstream_completion/schemas/"
        "operational_clock_packet_v1.schema.json",
}

EXPECTED_CONSUMER_SCHEMA_STATUS = {
    "canonical_action": "provisional_external_validation_schema",
    "full_yukawa": "provisional_external_validation_schema",
    "eft_matching": "provisional_external_validation_schema",
    "source_law_covariance": "nonpromoting_specification_schema",
    "operational_clock": "nonpromoting_specification_schema",
}

NATIVE_VERIFIER_COMMANDS = {
    "issue_565_carrier": [
        "code/a5_closure/echosahedral_selector_certificate.py",
        "verify",
        "--manifest",
        "code/a5_closure/manifests/echosahedral_federation_reference.json",
        "--receipt",
        "code/a5_closure/receipts/echosahedral_federation_reference.receipt.json",
    ],
    "issue_566_currents": [
        "code/a5_closure/port_current_inner_certificate.py",
        "verify",
        "--manifest",
        "code/a5_closure/manifests/port_current_response_reference.json",
        "--receipt",
        "code/a5_closure/receipts/port_current_inner_reference.receipt.json",
    ],
    "issue_314_matter": [
        "code/a5_closure/super_tannakian_matter_lift_certificate.py",
        "verify",
        "--manifest",
        "code/a5_closure/manifests/super_tannakian_matter_reference.json",
        "--receipt",
        "code/a5_closure/receipts/super_tannakian_matter_reference.receipt.json",
    ],
    "issue_567_global_form": [
        "code/a5_closure/axis_center_descent_certificate.py",
        "verify",
        "--manifest",
        "code/a5_closure/manifests/axis_center_descent_reference.json",
        "--receipt",
        "code/a5_closure/receipts/axis_center_descent_reference.receipt.json",
    ],
    "issue_569_family_band": [
        "code/a5_closure/family_band_attachment_certificate.py",
        "--verify",
    ],
    "issue_609_matter_boundary": [
        "code/a5_closure/matter_menu_spectral_ledger_certificate.py",
        "verify",
        "--manifest",
        "code/a5_closure/manifests/matter_menu_spectral_ledger_reference.json",
    ],
    "issue_311_finite_spectral_classical_control": [
        "code/a5_closure/flux_defect_criterion_certificate.py",
        "verify",
        "--manifest",
        "code/a5_closure/manifests/flux_defect_criterion_reference.json",
        "--receipt",
        "code/a5_closure/receipts/flux_defect_criterion_reference.receipt.json",
    ],
    "issue_32_rg_frontier": [
        "code/P_derivation/source_rg_frontier/check_rg_representation_frontier.py",
    ],
    "issue_630_scalar_yukawa_frontier": [
        "code/particles/hierarchy/higgs_yukawa_source_frontier/"
        "check_higgs_yukawa_source_frontier.py",
    ],
    "issue_631_local_carrier_frontier": [
        "code/a5_closure/check_common_ew_order_unit_carrier.py",
    ],
    "issue_633_source_clock_frontier": [
        "code/particles/hierarchy/source_clock_frontier/"
        "check_source_clock_frontier.py",
    ],
}


class FrontierVerificationError(ValueError):
    """Raised when the source-parent frontier fails closed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FrontierVerificationError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def byte_sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FrontierVerificationError(f"cannot read JSON {path}: {exc}") from exc


def safe_repo_file(repo_root: Path, relative_path: str) -> Path:
    require(isinstance(relative_path, str) and bool(relative_path), "empty repository path")
    pure = PurePosixPath(relative_path)
    require(not pure.is_absolute(), f"absolute path is forbidden: {relative_path}")
    require(".." not in pure.parts, f"path traversal is forbidden: {relative_path}")
    require("\\" not in relative_path, f"non-POSIX path is forbidden: {relative_path}")
    require(":" not in pure.parts[0], f"drive-qualified path is forbidden: {relative_path}")
    root = repo_root.resolve()
    cursor = root
    for part in pure.parts:
        cursor = cursor / part
        require(not cursor.is_symlink(), f"symlinked source input is forbidden: {relative_path}")
    resolved = cursor.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise FrontierVerificationError(
            f"path resolves outside repository: {relative_path}"
        ) from exc
    require(resolved.is_file(), f"referenced source file is missing: {relative_path}")
    return resolved


def validate_schema(inventory: Mapping[str, Any]) -> None:
    schema = load_json(SCHEMA_PATH)
    errors = sorted(
        Draft202012Validator(schema).iter_errors(inventory),
        key=lambda error: list(error.absolute_path),
    )
    if errors:
        details = "; ".join(
            f"{'/'.join(map(str, error.absolute_path)) or '$'}: {error.message}"
            for error in errors
        )
        raise FrontierVerificationError(f"inventory schema validation failed: {details}")


def _verify_file_pin(
    repo_root: Path,
    pin: Mapping[str, Any],
    expected_role: str,
    expected_schema: str,
) -> dict[str, Any]:
    require(pin["role"] == expected_role, f"wrong file role for {pin['path']}")
    require(pin["schema"] == expected_schema, f"wrong declared schema for {pin['path']}")
    path = safe_repo_file(repo_root, pin["path"])
    raw = path.read_bytes()
    payload = load_json(path)
    require(isinstance(payload, dict), f"source file is not an object: {pin['path']}")
    require(payload.get("schema") == expected_schema, f"payload schema drift: {pin['path']}")
    require(pin["bytes"] == len(raw), f"byte count mismatch: {pin['path']}")
    require(pin["byte_sha256"] == byte_sha256(raw), f"byte digest mismatch: {pin['path']}")
    require(
        pin["canonical_json_sha256"] == canonical_sha256(payload),
        f"canonical JSON digest mismatch: {pin['path']}",
    )
    if expected_role in {"semantic_artifact", "candidate_artifact"}:
        declared = payload.get("artifact_sha256")
        body = {key: value for key, value in payload.items() if key != "artifact_sha256"}
        require(
            declared == "sha256:" + canonical_sha256(body),
            f"semantic artifact self-digest mismatch: {pin['path']}",
        )
    if expected_role == "candidate_certificate":
        declared = payload.get("manifest_sha256")
        body = {key: value for key, value in payload.items() if key != "manifest_sha256"}
        require(
            declared == "sha256:" + canonical_sha256(body),
            f"candidate certificate self-digest mismatch: {pin['path']}",
        )
    return payload


def _verify_bindings(
    repo_root: Path,
    bindings: list[dict[str, Any]],
    expected: Mapping[str, Any],
    expected_status: str,
) -> dict[str, dict[str, Any]]:
    require(
        [binding["role"] for binding in bindings] == list(expected),
        "binding order or role set drifted",
    )
    loaded: dict[str, dict[str, Any]] = {}
    for binding in bindings:
        role = binding["role"]
        spec = expected[role]
        require(binding["issue"] == spec["issue"], f"wrong issue owner for {role}")
        require(
            binding["verifier_id"] == spec["verifier_id"],
            f"wrong native verifier for {role}",
        )
        require(binding["status"] == expected_status, f"wrong binding status for {role}")
        expected_files = spec["files"]
        require(
            [item["path"] for item in binding["files"]] == list(expected_files),
            f"file set or order drifted for {role}",
        )
        loaded[role] = {}
        for pin in binding["files"]:
            file_role, schema = expected_files[pin["path"]]
            loaded[role][pin["path"]] = _verify_file_pin(
                repo_root,
                pin,
                file_role,
                schema,
            )
        require(bool(binding["usable_exports"]), f"usable export scope missing for {role}")
        require(
            bool(binding["excluded_promotions"]),
            f"excluded-promotion boundary missing for {role}",
        )
    return loaded


def _verify_resolved_boundaries(
    repo_root: Path,
    boundaries: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    require(
        [item["boundary_id"] for item in boundaries]
        == list(EXPECTED_RESOLVED_BOUNDARIES),
        "resolved boundary order or identifier set drifted",
    )
    loaded: dict[str, dict[str, Any]] = {}
    for boundary in boundaries:
        boundary_id = boundary["boundary_id"]
        expected = EXPECTED_RESOLVED_BOUNDARIES[boundary_id]
        require(
            boundary["owner_issue"] == expected["issue"],
            f"wrong issue owner for resolved boundary {boundary_id}",
        )
        require(
            boundary["outcome"] == expected["outcome"],
            f"outcome drift for resolved boundary {boundary_id}",
        )
        require(
            boundary["verdict"] == expected["verdict"],
            f"verdict drift for resolved boundary {boundary_id}",
        )
        require(
            boundary["resolution_authority"]
            == "pinned_receipt_plus_live_github_issue_tracking",
            f"resolved boundary {boundary_id} lost its receipt authority",
        )
        require(
            len(boundary["files"]) == 1
            and boundary["files"][0]["path"] == expected["path"],
            f"receipt path drift for resolved boundary {boundary_id}",
        )
        payload = _verify_file_pin(
            repo_root,
            boundary["files"][0],
            "boundary_receipt",
            expected["schema"],
        )
        require(
            payload.get("issue") == expected["issue"],
            f"receipt issue drift for resolved boundary {boundary_id}",
        )
        require(
            payload.get("verdict") == expected["verdict"],
            f"receipt verdict drift for resolved boundary {boundary_id}",
        )
        require(
            payload.get("physical_promotion_allowed") is False,
            f"resolved boundary {boundary_id} was physically promoted",
        )
        require(
            payload.get(expected["sentinel"]) is True
            and payload.get("controls_fail_closed") is True
            and payload.get("blockers") == [],
            f"resolved boundary {boundary_id} lacks its attained control receipt",
        )
        loaded[boundary_id] = payload

    stage4_boundary = loaded["finite_local_domain_boundary"].get(
        "claim_boundary", ""
    ).lower()
    require(
        "finite causal and local-operator object" in stage4_boundary
        and "no continuum lorentzian spacetime" in stage4_boundary,
        "#634 boundary lost its finite-only continuum exclusion",
    )
    clock_boundary = loaded["physical_unit_boundary"].get(
        "claim_boundary", ""
    ).lower()
    require(
        "declared serialized interface" in clock_boundary
        and "not a transitive source-closure proof" in clock_boundary
        and "no-go theorem for every extended" in clock_boundary,
        "#633 receipt overstates the bounded serialized-interface verdict",
    )
    classical_boundary = loaded[
        "finite_spectral_quantum_discrimination_boundary"
    ].get("claim_boundary", "").lower()
    require(
        "two-component" in classical_boundary
        and "harmonic" in classical_boundary
        and "local-domain spectral interface" in classical_boundary
        and "no identity bridge" in classical_boundary
        and "does not close a criterion over an extended domain"
        in classical_boundary,
        "#311 receipt lost the bounded vector-spring spectral match",
    )
    classical_identity = loaded[
        "finite_spectral_quantum_discrimination_boundary"
    ].get("spectral_interface_identity", {})
    require(
        classical_identity.get("producer_schema")
        == "oph.local-domain-defect-sector-spectra.v1"
        and classical_identity.get(
            "rer_exact_flux_12_42_vertex_identity_bridge"
        )
        is False
        and classical_identity.get(
            "separate_from_rer_exact_flux_certificate"
        )
        is True
        and classical_identity.get("main_domain", {}).get(
            "visible_node_count"
        )
        == 8662
        and classical_identity.get("ladder_domain", {}).get(
            "visible_node_count"
        )
        == 1052,
        "#311 local-domain receipt was conflated with the exact flux support",
    )
    return loaded


def _verify_transitive_pins(loaded: Mapping[str, Mapping[str, Any]]) -> None:
    carrier_manifest_path = (
        "code/a5_closure/manifests/echosahedral_federation_reference.json"
    )
    carrier_receipt_path = (
        "code/a5_closure/receipts/echosahedral_federation_reference.receipt.json"
    )
    current_manifest_path = "code/a5_closure/manifests/port_current_response_reference.json"
    current_receipt_path = (
        "code/a5_closure/receipts/port_current_inner_reference.receipt.json"
    )
    response_artifact_path = (
        "code/a5_closure/manifests/charged_response_semantic_artifact.json"
    )
    matter_manifest_path = (
        "code/a5_closure/manifests/super_tannakian_matter_reference.json"
    )
    matter_receipt_path = (
        "code/a5_closure/receipts/super_tannakian_matter_reference.receipt.json"
    )
    spin_artifact_path = (
        "code/a5_closure/manifests/spin_statistics_semantic_artifact.json"
    )
    global_manifest_path = (
        "code/a5_closure/manifests/axis_center_descent_reference.json"
    )
    global_receipt_path = (
        "code/a5_closure/receipts/axis_center_descent_reference.receipt.json"
    )
    global_artifact_path = (
        "code/a5_closure/manifests/global_form_semantic_artifact.json"
    )

    carrier_manifest = loaded["screen_carrier"][carrier_manifest_path]
    carrier_receipt = loaded["screen_carrier"][carrier_receipt_path]
    carrier_digest = canonical_sha256(carrier_manifest)
    require(
        carrier_receipt.get("manifest_sha256") == carrier_digest,
        "#565 receipt does not bind its manifest",
    )
    require(carrier_receipt.get("issue") == 565, "#565 receipt issue tag drifted")
    require(
        carrier_receipt.get("source_firewall", {}).get("forbidden_dependency_hits") == [],
        "#565 source firewall is not clear",
    )

    current_manifest = loaded["finite_port_current_algebra"][current_manifest_path]
    current_receipt = loaded["finite_port_current_algebra"][current_receipt_path]
    response_artifact = loaded["finite_port_current_algebra"][response_artifact_path]
    current_digest = canonical_sha256(current_manifest)
    require(
        current_manifest.get("carrier_manifest_sha256") == carrier_digest,
        "#566 manifest does not bind #565",
    )
    require(
        current_receipt.get("manifest_sha256") == current_digest,
        "#566 receipt does not bind its manifest",
    )
    require(
        current_receipt.get("carrier_manifest_sha256") == carrier_digest,
        "#566 receipt carrier pin drifted",
    )
    require(current_receipt.get("issue") == 566, "#566 receipt issue tag drifted")
    require(
        current_manifest.get("semantic_response_artifact", {}).get("artifact_sha256")
        == response_artifact.get("artifact_sha256"),
        "#566 response artifact pin drifted",
    )
    require(
        current_receipt.get("claim_boundary", {}).get("does_not_close"),
        "#566 physical boundary is missing",
    )

    matter_manifest = loaded["finite_matter_module"][matter_manifest_path]
    matter_receipt = loaded["finite_matter_module"][matter_receipt_path]
    spin_artifact = loaded["finite_matter_module"][spin_artifact_path]
    matter_digest = canonical_sha256(matter_manifest)
    require(
        matter_manifest.get("current_manifest_sha256") == current_digest,
        "#314 manifest does not bind #566 manifest",
    )
    require(
        matter_manifest.get("current_receipt_sha256") == canonical_sha256(current_receipt),
        "#314 manifest does not bind #566 receipt",
    )
    require(
        matter_manifest.get("spin_statistics_artifact_sha256")
        == spin_artifact.get("artifact_sha256"),
        "#314 spin artifact pin drifted",
    )
    require(
        matter_receipt.get("manifest_sha256") == matter_digest,
        "#314 receipt does not bind its manifest",
    )
    require(matter_receipt.get("issue") == 314, "#314 receipt issue tag drifted")
    source_gate = matter_receipt.get("physical_source_gate", {})
    require(
        source_gate.get("declared_scalar_content_source_bound") is False,
        "#314 must not promote scalar existence",
    )
    require(
        source_gate.get("scalar_economy_source_bound") is False,
        "#314 must not promote scalar economy",
    )

    global_manifest = loaded["finite_global_form"][global_manifest_path]
    global_receipt = loaded["finite_global_form"][global_receipt_path]
    global_artifact = loaded["finite_global_form"][global_artifact_path]
    global_digest = canonical_sha256(global_manifest)
    require(
        global_manifest.get("matter_receipt_sha256") == canonical_sha256(matter_receipt),
        "#567 manifest does not bind #314 receipt",
    )
    require(
        global_manifest.get("global_form_artifact_sha256")
        == global_artifact.get("artifact_sha256"),
        "#567 global-form artifact pin drifted",
    )
    require(
        global_receipt.get("manifest_sha256") == global_digest,
        "#567 receipt does not bind its manifest",
    )
    require(
        global_receipt.get("matter_receipt_sha256") == canonical_sha256(matter_receipt),
        "#567 receipt matter pin drifted",
    )
    require(global_receipt.get("issue") == 567, "#567 receipt issue tag drifted")
    physical_gate = global_receipt.get("physical_global_form_gate", {})
    require(
        physical_gate.get("laboratory_global_form_attachment") is False,
        "#567 must not promote laboratory global-form attachment",
    )
    require(
        physical_gate.get("four_dimensional_instanton_action_normalization") is False,
        "#567 must not promote four-dimensional instanton normalization",
    )


def _verify_conditional_family(
    loaded: Mapping[str, Mapping[str, Any]],
) -> None:
    certificate_path = (
        "code/a5_closure/manifests/family_band_attachment_reference.json"
    )
    artifact_path = (
        "code/a5_closure/manifests/charged_response_pole_residue_artifact.json"
    )
    attachment_path = "code/a5_closure/manifests/matter_attachment_receipt.json"
    certificate = loaded["family_band_candidate"][certificate_path]
    artifact = loaded["family_band_candidate"][artifact_path]
    attachment = loaded["family_band_candidate"][attachment_path]
    require(certificate.get("issue") == 569, "#569 certificate issue tag drifted")
    require(
        certificate.get("named_interface", {}).get("class")
        == "conditional_open_interface",
        "#569 candidate was relabeled as a positive physical parent",
    )
    open_receipts = set(certificate.get("named_interface", {}).get("open_receipts", []))
    require(
        {
            "matter-pole identification",
            "physical Spin/locality bridge",
            "physical seam action selection",
            "laboratory current identification",
        }.issubset(open_receipts),
        "#569 physical-family boundary is incomplete",
    )
    require(
        certificate.get("open_gates")
        == [
            "matter_pole_identification",
            "continuum_Spin_locality",
            "physical_seam_action_selection",
            "laboratory_current_identification",
        ],
        "#569 open physical gates drifted",
    )
    require(
        certificate.get("promotion")
        == {
            "matter_pole_identified": False,
            "continuum_spin_locality_derived": False,
            "physical_seam_action_selected": False,
            "laboratory_current_identified": False,
            "promotion_allowed": False,
        },
        "#569 conditional finite attachment was physically promoted",
    )
    require(
        "matter-pole identification"
        in certificate.get("claim_boundary", "").lower(),
        "#569 claim boundary lost the matter-pole gate",
    )
    artifact_boundary = artifact.get("claim_boundary", "")
    require(
        "not a matter-pole measurement" in artifact_boundary
        and "matter-pole identification" in artifact_boundary,
        "#569 pole-residue artifact was overpromoted",
    )
    require(
        attachment.get("issue") == 569
        and attachment.get("verdict") == "ATTAINED"
        and attachment.get("MATTER_ATTACHMENT_RECEIPT") is True
        and attachment.get("controls_fail_closed") is True
        and attachment.get("blockers") == []
        and attachment.get("physical_promotion_allowed") is False,
        "#569 finite-domain attachment status drifted",
    )
    require(
        attachment.get("attachment", {}).get("complex_rank") == 45
        and attachment.get("attachment", {}).get("band_rank_measured") == 3,
        "#569 finite rank-forty-five attachment drifted",
    )
    structure = attachment.get("gap_inheritance_certificate", {}).get("structure", "")
    matter_operator = attachment.get("matter_operator_certificate", {})
    gap_inheritance = attachment.get("gap_inheritance_certificate", {})
    require(
        "scalar operator tensor identity" in structure
        and "multiplicity forty-five" in structure,
        "#569 matter operator lost its tensor-identity limitation",
    )
    require(
        attachment.get("declared_matter_packet", {}).get("source_selected")
        is False
        and matter_operator.get("status") == "declared_tensor_extension"
        and matter_operator.get("source_selected") is False
        and gap_inheritance.get("status")
        == (
            "conditional_algebraic_inheritance_under_declared_"
            "tensor_extension"
        )
        and gap_inheritance.get("matter_action_source_selected") is False,
        "#569 declared tensor extension was relabeled as source-selected",
    )
    spin_layer = attachment.get("spin_layer", {})
    require(
        spin_layer.get("packet_status")
        == "separate_pinned_issue_314_packet"
        and spin_layer.get("spin_to_local_domain_bridge_certified") is False
        and spin_layer.get("same_source_domain_certified") is False
        and spin_layer.get("open_interface")
        == "physical Spin/locality bridge",
        "#569 issue-314 spin packet was silently attached to the issue-634 "
        "local domain",
    )
    attachment_clauses = attachment.get("clause_verdicts", {})
    require(
        attachment_clauses.get(
            "separate_issue_314_spin_packet_resolved"
        )
        is True
        and attachment_clauses.get(
            "local_domain_stage2_context_recorded"
        )
        is True
        and attachment_clauses.get(
            "local_stage2_same_source_domain_binding"
        )
        is True
        and attachment_clauses.get(
            "conditional_gap_inheritance_exact"
        )
        is True
        and "same_source_domain_binding" not in attachment_clauses
        and "spin_gates_consumed" not in attachment_clauses
        and "gap_inherited_exact" not in attachment_clauses,
        "#569 attained clauses conflate the spin and local-domain packets",
    )
    attachment_boundary = attachment.get("claim_boundary", "").lower()
    bounded_scan = attachment.get("bounded_declared_key_scan", {})
    local_parent_pins = attachment.get("upstream_pins", {}).get(
        "local_domain_parent_sha256",
        {},
    )
    require(
        "source does not select a matter action" in attachment_boundary
        and "no source, domain, or transport bridge" in attachment_boundary
        and "physical seam-action selection" in attachment_boundary
        and "matter-pole identification" in attachment_boundary
        and "bounded declared-key scan" in attachment_boundary
        and "not semantic input closure" in attachment_boundary
        and bounded_scan.get("fragments")
        == ["yukawa", "pole_mass", "mass_gev", "mev"]
        and bounded_scan.get("hits") == []
        and "declared mapping keys only"
        in str(bounded_scan.get("scope", "")).lower()
        and "no transitive input-closure claim"
        in str(bounded_scan.get("scope", "")).lower()
        and set(local_parent_pins)
        == {
            "source_gap_receipt.json",
            "stage1_arrays.npz.gz",
            "stage1_receipt.json",
            "stage2_receipt.json",
            "stage3_receipt.json",
        }
        and all(
            isinstance(value, str)
            and value.startswith("sha256:")
            and len(value) == 71
            for value in local_parent_pins.values()
        ),
        "#569 finite attachment was overpromoted to a physical matter action",
    )


def _verify_matter_completeness_boundary(
    loaded: Mapping[str, Mapping[str, Any]],
) -> None:
    path = "code/a5_closure/manifests/matter_menu_spectral_ledger_reference.json"
    boundary = loaded["matter_completeness_boundary"][path]
    require(
        boundary.get("verdicts", {}).get("menu_completeness_inside_declared_algebra")
        == "exact",
        "#609 in-algebra completeness verdict drifted",
    )
    require(
        boundary.get("verdicts", {}).get("beyond_declared_algebra")
        == "independence_limited",
        "#609 beyond-algebra boundary was overpromoted",
    )
    sterile = boundary.get("off_menu_controls", {}).get(
        "neutral_singlet_sterile", {}
    )
    require(
        sterile.get("verdict") == "sterile_countermodel_source_invisible"
        and sterile.get("all_current_observable_couplings_zero") is True,
        "#609 sterile countermodel is missing",
    )
    require(
        boundary.get("light_heavy_threshold", {})
        .get("physical_decoupling_interface", {})
        .get("status")
        == "separate_open_physical_interface",
        "#609 declared threshold was relabeled as physical decoupling",
    )


def _verify_finite_spectral_classical_control(
    loaded: Mapping[str, Mapping[str, Any]],
) -> None:
    path = "code/a5_closure/receipts/flux_defect_criterion_reference.receipt.json"
    receipt = loaded["finite_flux_and_local_domain_classical_context"][path]
    require(receipt.get("issue") == 311, "#311 receipt issue tag drifted")
    exact = receipt.get("exact_support_classical_realification", {})
    require(
        exact.get("domain")
        == {
            "name": "certified_icosahedral_support",
            "vertex_count": 12,
            "seam_count": 30,
            "regular_degree": 5,
            "separate_from_local_domain_seam_complex": True,
        }
        and exact.get("stiffness_rule") == "K_k = 5 I - A_k"
        and exact.get("declared_adjacency_spectral_family_recoverable") is True
        and exact.get("scalar_operator_or_gap_matched") is False
        and exact.get("phase_metric_isometry_checks") == 180
        and exact.get("edge_hessian_identity_entry_checks") == 3456
        and len(exact.get("per_class", [])) == 6,
        "#311 exact-support vector-spring realification is missing",
    )
    require(
        all(
            row.get("classical_hessian_certificate", {}).get(
                "coordinate_metric"
            )
            == [[2, 1], [1, 2]]
            and row.get("classical_hessian_certificate", {}).get(
                "metric_leading_principal_minors"
            )
            == [2, 3]
            and row.get("classical_hessian_certificate", {}).get(
                "energy_hessian_equals_metric_times_stiffness"
            )
            is True
            and row.get("classical_hessian_certificate", {}).get(
                "positive_semidefinite_by_edge_sum_of_squares"
            )
            is True
            for row in exact.get("per_class", [])
        ),
        "#311 classical completion lacks the exact positive-metric "
        "edge-Hessian proof",
    )
    require(
        exact.get("complete_interface_ontology_no_go") is False
        and exact.get("extended_domain_non_identifiability") is False,
        "#311 exact-support realification was expanded into a complete no-go",
    )
    local = receipt.get("local_domain_classical_spectral_context", {})
    require(
        local.get("verdict")
        == "CLASSICAL_REALIZATION_MATCHES_DECLARED_FINITE_SPECTRAL_INTERFACE"
        and "two-component classical harmonic network"
        in str(local.get("classical_model", "")).lower()
        and local.get("sector_payload_identity") is True
        and local.get("scalar_gap_payload_identity") is True
        and local.get("ladder_payload_identity") is True
        and local.get("separate_from_exact_flux_support") is True
        and local.get("exact_flux_identity_bridge") is False,
        "#311 separate local-domain vector-spring context is missing",
    )
    require(
        local.get("complete_interface_ontology_no_go") is False
        and local.get("extended_domain_non_identifiability") is False
        and local.get("bounded_lexical_census", {}).get("completeness_theorem")
        is False,
        "#311 local-domain spectral context was expanded into a complete no-go",
    )
    criteria = receipt.get("acceptance_criteria_status", {})
    require(
        criteria.get("quantum_pole_or_equivalent_physical_spectral_criterion_proved")
        is False
        and criteria.get("mass_invariant_and_target_independent") is False,
        "#311 finite receipt fabricated a pole or mass",
    )
    boundary = receipt.get("claim_boundary", {})
    require(
        "complete-interface" in " ".join(boundary.get("does_not_close", []))
        and "extended source domain"
        in " ".join(boundary.get("does_not_close", []))
        and "continuum quantum pole" in " ".join(boundary.get("does_not_close", [])),
        "#311 bounded spectral control lost its physical exclusions",
    )
    require(
        receipt.get("local_domain_spectral_context_boundary", {}).get(
            "identity_bridge_to_exact_flux_support"
        )
        is False,
        "#311 local-domain and exact-support spectra were silently identified",
    )


def _verify_rg_representation_frontier(
    loaded: Mapping[str, Mapping[str, Any]],
) -> None:
    path = (
        "code/P_derivation/source_rg_frontier/outputs/"
        "rg_representation_frontier.json"
    )
    frontier = loaded["rg_representation_frontier"][path]
    require(
        frontier.get("status")
        == "PARTIAL_EXACT_REPRESENTATION_INDICES__SOURCE_MATCHING_OPEN",
        "#32 source frontier status drifted",
    )
    require(frontier.get("promotion_allowed") is False, "#32 frontier promoted itself")
    require(
        frontier.get("qft_import_boundary", {}).get("oph_native_one_loop_beta_theorem")
        is False,
        "#32 imported QFT functional was relabeled as OPH-native",
    )
    require(
        all(
            item.get("status") == "not_emitted"
            for item in frontier.get("matching_objects", {}).values()
        ),
        "#32 frontier fabricated a matching object",
    )


def _verify_scalar_yukawa_source_frontier(
    loaded: Mapping[str, Mapping[str, Any]],
) -> None:
    path = (
        "code/particles/hierarchy/higgs_yukawa_source_frontier/outputs/"
        "higgs_yukawa_source_frontier.json"
    )
    frontier = loaded["scalar_yukawa_source_frontier"][path]
    require(frontier.get("issue") == 630, "#630 frontier has the wrong owner")
    require(
        frontier.get("status")
        == "BOUNDED_NONPROMOTING_FRONTIER__POSITIVE_SOURCE_ACTION_OPEN",
        "#630 frontier status drifted",
    )
    require(frontier.get("promotion_allowed") is False, "#630 frontier promoted itself")
    require(
        frontier.get("coefficient_assignments_emitted") is False
        and frontier.get("physical_source_action_emitted") is False,
        "#630 frontier fabricated a physical action or coefficient assignment",
    )
    require(
        all(row.get("status") != "complete" for row in frontier.get("acceptance_map", [])),
        "#630 acceptance was overpromoted",
    )


def _verify_local_ew_order_unit_frontier(
    loaded: Mapping[str, Mapping[str, Any]],
) -> None:
    path = "code/a5_closure/manifests/common_ew_order_unit_carrier_reference.json"
    frontier = loaded["local_ew_order_unit_frontier"][path]
    require(frontier.get("issue") == 631, "#631 frontier has the wrong owner")
    require(
        frontier.get("status")
        == "FINITE_ORDER_UNIT_INTERTWINER__PHYSICAL_COMMON_CARRIER_OPEN",
        "#631 frontier status drifted",
    )
    promotion = frontier.get("promotion", {})
    require(
        promotion
        and all(value is False for value in promotion.values()),
        "#631 finite line isomorphism was given a physical promotion",
    )
    require(
        all(row.get("status") == "open" for row in frontier.get("open_gates", [])),
        "#631 physical gate was self-attested closed",
    )


def _verify_source_clock_frontier(
    loaded: Mapping[str, Mapping[str, Any]],
) -> None:
    path = (
        "code/particles/hierarchy/source_clock_frontier/outputs/"
        "source_clock_frontier.json"
    )
    frontier = loaded["source_clock_frontier"][path]
    require(frontier.get("issue") == 633, "#633 frontier has the wrong owner")
    require(
        frontier.get("status")
        == "NONPROMOTING_ROUTE_NEUTRAL_CONTRACT__PHYSICAL_CLOCK_ATTACHMENT_OPEN",
        "#633 frontier status drifted",
    )
    require(
        frontier.get("dependency_semantics")
        == {
            "hard_issue_dependencies": [634],
            "alternative_routes_allowed": True,
            "optional_route_owner_issues": [
                32,
                34,
                317,
                318,
                425,
                522,
                545,
                546,
                569,
                633,
            ],
            "downstream_only_issues": [334],
        },
        "#633 route dependency semantics drifted",
    )
    routes = frontier.get("candidate_routes", [])
    require(
        len(routes) == 1
        and routes[0].get("route_id") == "cesium_133_hyperfine"
        and routes[0].get("required_for_issue_633") is False,
        "#633 optional cesium route became mandatory or malformed",
    )
    require(
        frontier.get("dimensionless_clock_gap_emitted") is False
        and frontier.get("source_energy_interval_emitted") is False
        and frontier.get("source_g_si_interval_emitted") is False
        and frontier.get("physical_promotion_allowed") is False,
        "#633 frontier fabricated a clock, source energy, or gravity output",
    )
    require(
        frontier.get("downstream_consumers")
        == [
            {
                "issue": 334,
                "role": "newton_g_composition",
                "required_input": "source_energy_interval",
                "possible_output": "source_g_si_interval",
                "status": "downstream_open_not_issue_633_acceptance_gate",
            }
        ],
        "#334 was not confined to the downstream gravity composition",
    )


def _verify_policy(repo_root: Path, inventory: Mapping[str, Any]) -> dict[str, Any]:
    path = safe_repo_file(repo_root, POLICY_REL)
    raw = path.read_bytes()
    policy = load_json(path)
    require(
        canonical_sha256(policy) == EXPECTED_POLICY_CANONICAL_SHA256,
        "trusted source-parent policy digest drifted",
    )
    pin = inventory["policy"]
    require(pin["path"] == POLICY_REL, "inventory points to a different policy")
    require(pin["bytes"] == len(raw), "policy byte count mismatch")
    require(pin["byte_sha256"] == byte_sha256(raw), "policy byte digest mismatch")
    require(
        pin["canonical_json_sha256"] == canonical_sha256(policy),
        "policy canonical digest mismatch",
    )
    require(
        [item["role"] for item in policy["positive_parents"]] == list(EXPECTED_POSITIVE),
        "policy positive-parent roles drifted",
    )
    require(
        [item["role"] for item in policy["conditional_context"]]
        == list(EXPECTED_CONDITIONAL),
        "policy conditional-context roles drifted",
    )
    require(
        [item["boundary_id"] for item in policy["resolved_boundaries"]]
        == list(EXPECTED_RESOLVED_BOUNDARIES),
        "policy resolved-boundary identifiers drifted",
    )
    return policy


def _verify_consumer_schemas(
    repo_root: Path,
    schemas: list[dict[str, Any]],
) -> None:
    require(
        [item["slot"] for item in schemas] == list(EXPECTED_CONSUMER_SCHEMAS),
        "consumer schema slot order drifted",
    )
    for item in schemas:
        expected_path = EXPECTED_CONSUMER_SCHEMAS[item["slot"]]
        require(item["path"] == expected_path, f"consumer schema path drift: {item['slot']}")
        require(
            item["status"] == EXPECTED_CONSUMER_SCHEMA_STATUS[item["slot"]],
            f"consumer schema status drift: {item['slot']}",
        )
        path = safe_repo_file(repo_root, expected_path)
        raw = path.read_bytes()
        payload = load_json(path)
        require(
            payload.get("$schema") == "https://json-schema.org/draft/2020-12/schema",
            f"consumer schema is not Draft 2020-12: {expected_path}",
        )
        require(item["bytes"] == len(raw), f"consumer schema bytes drift: {expected_path}")
        require(
            item["byte_sha256"] == byte_sha256(raw),
            f"consumer schema hash drift: {expected_path}",
        )


def _verify_dependency_state(
    inventory: Mapping[str, Any],
    policy: Mapping[str, Any],
) -> None:
    for inventory_key, policy_key in (
        ("positive_parent_bindings", "positive_parents"),
        ("conditional_context", "conditional_context"),
    ):
        for binding, source in zip(
            inventory[inventory_key],
            policy[policy_key],
            strict=True,
        ):
            require(
                binding["role"] == source["role"]
                and binding["usable_exports"] == source["usable_exports"]
                and binding["excluded_promotions"] == source["excluded_promotions"],
                f"declared export scope drifted for {source['role']}",
            )
    require(
        [item.get("context_status") for item in inventory["conditional_context"]]
        == [item["status"] for item in policy["conditional_context"]],
        "conditional-context issue outcomes drifted",
    )
    expected_boundaries = inventory["resolved_boundaries"]
    require(
        len(expected_boundaries) == len(policy["resolved_boundaries"]),
        "resolved Phase-2 boundary count drifted",
    )
    for boundary, item in zip(
        expected_boundaries,
        policy["resolved_boundaries"],
        strict=True,
    ):
        require(
            {
                "boundary_id": boundary["boundary_id"],
                "owner_issue": boundary["owner_issue"],
                "outcome": boundary["outcome"],
                "scope": boundary["scope"],
                "effect": boundary["effect"],
                "scientific_parent_status": boundary["scientific_parent_status"],
                "resolution_authority": boundary["resolution_authority"],
                "verdict": boundary["verdict"],
            }
            == {
                "boundary_id": item["boundary_id"],
                "owner_issue": item["owner_issue"],
                "outcome": item["outcome"],
                "scope": item["scope"],
                "effect": item["effect"],
                "scientific_parent_status": item["scientific_parent_status"],
                "resolution_authority": item["resolution_authority"],
                "verdict": item["expected_verdict"],
            },
            "resolved Phase-2 boundary state drifted",
        )
    expected_interfaces = [
        {
            "gate_id": item["gate_id"],
            "owner_issues": item["owner_issues"],
            "required_output": item["required_output"],
            "prerequisites": item["prerequisites"],
            "terminal_for_dimensionless_output":
                item["terminal_for_dimensionless_output"],
            "status": "open",
            "evidence": [],
        }
        for item in policy["open_interfaces"]
    ]
    require(
        inventory["open_interfaces"] == expected_interfaces,
        "open source-interface dependency state drifted",
    )

    dag = inventory["source_dag"]
    expected_nodes = [
        {
            "id": binding["role"],
            "class": "verified_finite_parent",
            "issue": binding["issue"],
            "status": "verified_finite_parent",
        }
        for binding in inventory["positive_parent_bindings"]
    ]
    expected_nodes.extend(
        {
            "id": binding["role"],
            "class": "conditional_context",
            "issue": binding["issue"],
            "status": binding["context_status"],
        }
        for binding in inventory["conditional_context"]
    )
    expected_nodes.extend(
        {
            "id": item["boundary_id"],
            "class": "resolved_boundary",
            "issue": item["owner_issue"],
            "status": item["outcome"],
        }
        for item in expected_boundaries
    )
    expected_nodes.extend(
        {
            "id": item["gate_id"],
            "class": "open_interface",
            "issue": item["owner_issues"][0],
            "status": "open",
        }
        for item in expected_interfaces
    )
    expected_nodes.append(
        {
            "id": "oph_native_dimensionless_packet",
            "class": "candidate_output",
            "issue": 594,
            "status": "blocked",
        }
    )
    require(dag["nodes"] == expected_nodes, "source DAG node state drifted")

    expected_edges = [
        {"from": "screen_carrier", "to": "finite_port_current_algebra"},
        {"from": "finite_port_current_algebra", "to": "finite_matter_module"},
        {"from": "finite_matter_module", "to": "finite_global_form"},
        {"from": "screen_carrier", "to": "family_band_candidate"},
        {"from": "finite_matter_module", "to": "family_band_candidate"},
    ]
    expected_edges.extend(
        {
            "from": binding["role"],
            "to": "oph_native_dimensionless_packet",
        }
        for binding in inventory["positive_parent_bindings"]
    )
    for interface in expected_interfaces:
        expected_edges.extend(
            {
                "from": prerequisite,
                "to": interface["gate_id"],
            }
            for prerequisite in interface["prerequisites"]
        )
        if interface["terminal_for_dimensionless_output"]:
            expected_edges.append(
                {
                    "from": interface["gate_id"],
                    "to": "oph_native_dimensionless_packet",
                }
            )
    require(dag["edges"] == expected_edges, "source DAG dependency edges drifted")
    require(
        {
            "from": "common_screen_electroweak_carrier",
            "to": "source_complete_field_census",
        }
        not in dag["edges"],
        "#632 incorrectly depends on #631",
    )


def _acyclic_and_forbidden_paths(dag: Mapping[str, Any]) -> tuple[bool, list[list[str]]]:
    nodes = dag["nodes"]
    edges = dag["edges"]
    node_ids = [node["id"] for node in nodes]
    require(len(node_ids) == len(set(node_ids)), "source DAG has duplicate nodes")
    known = set(node_ids)
    require(
        all(edge["from"] in known and edge["to"] in known for edge in edges),
        "source DAG has dangling edges",
    )
    outgoing = {node_id: [] for node_id in node_ids}
    indegree = {node_id: 0 for node_id in node_ids}
    incoming = {node_id: [] for node_id in node_ids}
    for edge in edges:
        outgoing[edge["from"]].append(edge["to"])
        incoming[edge["to"]].append(edge["from"])
        indegree[edge["to"]] += 1
    queue = sorted(node_id for node_id, degree in indegree.items() if degree == 0)
    visited: list[str] = []
    while queue:
        current = queue.pop(0)
        visited.append(current)
        for target in sorted(outgoing[current]):
            indegree[target] -= 1
            if indegree[target] == 0:
                queue.append(target)
                queue.sort()
    acyclic = len(visited) == len(node_ids)

    classes = {node["id"]: node["class"] for node in nodes}
    forbidden = set(dag["forbidden_node_classes"])
    paths: list[list[str]] = []

    def visit(current: str, trail: list[str]) -> None:
        if current in trail:
            paths.append(list(reversed(trail + [current])))
            return
        if classes[current] in forbidden:
            paths.append(list(reversed(trail + [current])))
            return
        for parent in incoming[current]:
            visit(parent, trail + [current])

    for output in dag["protected_outputs"]:
        require(output in known, f"unknown protected output: {output}")
        visit(output, [])
    return acyclic, sorted(paths)


def _structured_token_hits(value: Any, tokens: list[str], at: str = "$") -> list[str]:
    hits: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            key_text = str(key).lower()
            for token in tokens:
                pattern = rf"(?<![a-z0-9_]){re.escape(token)}(?![a-z0-9_])"
                if re.search(pattern, key_text):
                    hits.append(f"{at}.{key}:{token}")
            hits.extend(_structured_token_hits(child, tokens, f"{at}.{key}"))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            hits.extend(_structured_token_hits(child, tokens, f"{at}[{index}]"))
    elif isinstance(value, str):
        text = value.lower()
        for token in tokens:
            pattern = rf"(?<![a-z0-9_]){re.escape(token)}(?![a-z0-9_])"
            if re.search(pattern, text):
                hits.append(f"{at}:{token}")
    return hits


def _under_allowed_root(path: str, roots: list[str]) -> bool:
    pure = PurePosixPath(path)
    return any(
        pure == PurePosixPath(root) or PurePosixPath(root) in pure.parents
        for root in roots
    )


def _verify_firewall(
    inventory: Mapping[str, Any],
    policy: Mapping[str, Any],
    loaded_positive: Mapping[str, Mapping[str, Any]],
    loaded_conditional: Mapping[str, Mapping[str, Any]],
    loaded_boundaries: Mapping[str, Mapping[str, Any]],
) -> None:
    firewall = inventory["target_firewall"]
    expected = policy["firewall_policy"]
    require(
        firewall["allowed_source_roots"] == expected["allowed_source_roots"],
        "firewall allowed roots drifted",
    )
    require(
        firewall["forbidden_source_paths"] == expected["forbidden_source_paths"],
        "firewall forbidden paths drifted",
    )
    require(
        firewall["forbidden_structured_tokens"]
        == expected["forbidden_structured_tokens"],
        "firewall token policy drifted",
    )
    resolved = firewall["resolved_source_paths"]
    require(len(resolved) == len(set(resolved)), "firewall resolved path list has duplicates")
    require(
        all(_under_allowed_root(path, firewall["allowed_source_roots"]) for path in resolved),
        "resolved source path lies outside allowlisted roots",
    )
    require(
        not set(resolved).intersection(firewall["forbidden_source_paths"]),
        "a forbidden target file is mounted as source input",
    )
    expected_resolved = sorted(
        [
            path
            for group in (EXPECTED_POSITIVE, EXPECTED_CONDITIONAL)
            for spec in group.values()
            for path in spec["files"]
        ]
        + [
            spec["path"]
            for spec in EXPECTED_RESOLVED_BOUNDARIES.values()
        ]
        + list(EXPECTED_CONSUMER_SCHEMAS.values())
    )
    require(resolved == expected_resolved, "resolved source-input closure drifted")

    tokens = [token.lower() for token in firewall["forbidden_structured_tokens"]]
    hits: list[str] = []
    for group in (loaded_positive, loaded_conditional):
        for role, files in group.items():
            for path, payload in files.items():
                hits.extend(
                    f"{role}:{path}:{hit}"
                    for hit in _structured_token_hits(payload, tokens)
                )
    for boundary_id, payload in loaded_boundaries.items():
        hits.extend(
            f"{boundary_id}:{hit}"
            for hit in _structured_token_hits(payload, tokens)
        )
    require(not hits, f"structured target content found in source inputs: {hits}")
    require(firewall["comparison_channel_present"] is False, "comparison channel mounted")
    require(firewall["network_required"] is False, "source producer requires network")
    require(
        firewall["runtime_execution_receipt"] == "open",
        "runtime isolation must remain open until receipt-backed",
    )
    require(
        firewall["human_formula_selection_ancestry"] == "open",
        "human formula-selection ancestry was self-attested closed",
    )
    require(firewall["full_gate_satisfied"] is False, "target firewall overpromoted")


def _bundle_seed(bindings: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "role": binding["role"],
            "issue": binding["issue"],
            "files": [
                {
                    "path": item["path"],
                    "canonical_json_sha256": item["canonical_json_sha256"],
                }
                for item in binding["files"]
            ],
        }
        for binding in bindings
    ]


def _run_native_verifiers(repo_root: Path) -> dict[str, str]:
    statuses: dict[str, str] = {}
    for verifier_id, argv in NATIVE_VERIFIER_COMMANDS.items():
        result = subprocess.run(
            [sys.executable, *argv],
            cwd=repo_root,
            check=False,
            capture_output=True,
            text=True,
            env={
                "PATH": str(Path(sys.executable).parent),
                "PYTHONDONTWRITEBYTECODE": "1",
            },
        )
        if result.returncode != 0:
            raise FrontierVerificationError(
                f"native verifier {verifier_id} failed: "
                f"{result.stdout[-1000:]}{result.stderr[-1000:]}"
            )
        statuses[verifier_id] = "PASS"
    return statuses


def verify_inventory(
    inventory: Mapping[str, Any],
    *,
    repo_root: Path = DEFAULT_REPO_ROOT,
    run_native_verifiers: bool = True,
) -> dict[str, Any]:
    repo_root = repo_root.resolve()
    validate_schema(inventory)
    policy = _verify_policy(repo_root, inventory)
    positive = _verify_bindings(
        repo_root,
        inventory["positive_parent_bindings"],
        EXPECTED_POSITIVE,
        "verified_finite_parent",
    )
    conditional = _verify_bindings(
        repo_root,
        inventory["conditional_context"],
        EXPECTED_CONDITIONAL,
        "conditional_context_only",
    )
    _verify_transitive_pins(positive)
    _verify_conditional_family(conditional)
    _verify_matter_completeness_boundary(conditional)
    _verify_finite_spectral_classical_control(conditional)
    _verify_rg_representation_frontier(conditional)
    _verify_scalar_yukawa_source_frontier(conditional)
    _verify_local_ew_order_unit_frontier(conditional)
    _verify_source_clock_frontier(conditional)
    resolved_boundaries = _verify_resolved_boundaries(
        repo_root,
        inventory["resolved_boundaries"],
    )
    _verify_consumer_schemas(repo_root, inventory["consumer_contract"]["schemas"])
    _verify_firewall(
        inventory,
        policy,
        positive,
        conditional,
        resolved_boundaries,
    )
    _verify_dependency_state(inventory, policy)

    bridge = inventory["coordinate_bridge"]
    require(
        bridge == {
            "source_coordinate": "v_chart",
            "renormalized_coordinate": "v_F",
            "equality_receipt": None,
            "status": "open",
            "relabel_allowed": False,
        },
        "v_chart to v_F bridge was relabeled without a proof receipt",
    )
    require(
        inventory["unit_scope"]
        == {
            "native_coordinates": "E_star_normalized_dimensionless",
            "physical_clock_required": True,
            "physical_unit_verdict":
                "PHYSICAL_UNITS_NOT_EVALUABLE_ON_DECLARED_SERIALIZED_INTERFACE",
            "dimensionful_values_present": False,
            "emitted_observables": [],
        },
        "physical-unit row drifted or the frontier emitted a dimensionful pole",
    )
    require(inventory["promotion_allowed"] is False, "frontier promoted itself")
    require(
        inventory["consumer_contract"]["frozen_algorithm_substitution_ready"] is False,
        "#593 substitution was marked ready before validation",
    )
    require(
        inventory["consumer_contract"]["common_subject_digest_ready"] is False,
        "common subject digest was marked ready before integration",
    )

    indices = [row["acceptance_index"] for row in inventory["acceptance_map"]]
    require(indices == list(range(1, 10)), "acceptance map must cover rows 1 through 9")
    expected_acceptance = {
        1: (
            "closed_bounded",
            [],
            "issue 634 supplies a finite causal and local-operator domain; no continuum Lorentzian or quantum-EFT promotion follows",
        ),
        2: (
            "closed_not_evaluable",
            [],
            "issue 633 closes the physical-unit row as PHYSICAL_UNITS_NOT_EVALUABLE_ON_DECLARED_SERIALIZED_INTERFACE under a bounded schema/source scan and named channel-nonuse experiment; no complete-domain clock non-identifiability theorem follows",
        ),
        3: (
            "open",
            [
                "physical_family_and_matter_pole_attachment",
                "scalar_yukawa_fj_integration",
                "common_screen_electroweak_carrier",
                "source_complete_field_census",
            ],
            "the finite rank-forty-five candidate is bound, but no source-selected scalar action, complete Yukawa matrices, physical family action, or coupled-sector census is emitted",
        ),
        4: (
            "open",
            [
                "source_to_fj_coordinate_map",
                "scalar_yukawa_fj_integration",
            ],
            "v_chart and v_F remain distinct typed coordinates; the source-to-FJ map is owned by issue 638 and its integration by issue 630",
        ),
        5: (
            "partial",
            ["target_clean_rg_threshold_matching"],
            "the issue-32 frontier supplies exact per-copy representation indices and a parametric one-loop gauge law under an imported QFT functional; ordered intervals, thresholds, finite maps, Jacobians, masks, and remainders are absent",
        ),
        6: (
            "open",
            ["finite_to_lorentzian_quantum_eft_transfer"],
            "the finite domain has no certified transfer to a Lorentzian Spin quantum EFT accepted by the pole consumer",
        ),
        7: (
            "open",
            ["unique_source_root_and_joint_law", "validated_qft_consumer"],
            "substitution into unchanged validated QFT algorithms waits for the source packet and the production consumer",
        ),
        8: (
            "partial",
            ["runtime_and_human_target_firewall"],
            "source paths and structured ancestry are allowlisted and target-free; hermetic runtime and human-selection receipts remain open",
        ),
        9: (
            "open",
            [
                "unique_source_root_and_joint_law",
                "validated_qft_consumer",
                "runtime_and_human_target_firewall",
            ],
            "independent replay of the complete native source-to-consumer conjunction requires every open producer and firewall gate",
        ),
    }
    for row in inventory["acceptance_map"]:
        require(
            (row["status"], row["blocking_gates"], row["summary"])
            == expected_acceptance[row["acceptance_index"]],
            (
                f"acceptance row {row['acceptance_index']} status, blockers, "
                "or claim scope drifted"
            ),
        )

    dag_acyclic, forbidden_paths = _acyclic_and_forbidden_paths(inventory["source_dag"])
    require(dag_acyclic is True, "source DAG is cyclic")
    require(forbidden_paths == [], f"source DAG has forbidden ancestry: {forbidden_paths}")
    require(inventory["source_dag"]["acyclic"] is True, "displayed DAG result drifted")
    require(
        inventory["source_dag"]["forbidden_ancestry_paths"] == [],
        "displayed forbidden-ancestry result drifted",
    )

    require(
        inventory["finite_parent_bundle_digest"]
        == canonical_sha256(_bundle_seed(inventory["positive_parent_bindings"])),
        "finite parent bundle digest mismatch",
    )
    require(
        inventory["conditional_context_digest"]
        == canonical_sha256(_bundle_seed(inventory["conditional_context"])),
        "conditional context digest mismatch",
    )
    body = {key: value for key, value in inventory.items() if key != "inventory_digest"}
    require(
        inventory["inventory_digest"] == canonical_sha256(body),
        "inventory self-digest mismatch",
    )

    native = _run_native_verifiers(repo_root) if run_native_verifiers else {
        verifier_id: "SKIPPED_BY_CALLER"
        for verifier_id in NATIVE_VERIFIER_COMMANDS
    }
    return {
        "status": "PASS",
        "issue": 594,
        "claim_lane": inventory["claim_lane"],
        "positive_parent_count": len(inventory["positive_parent_bindings"]),
        "conditional_context_count": len(inventory["conditional_context"]),
        "open_interface_count": len(inventory["open_interfaces"]),
        "promotion_allowed": False,
        "native_verifiers": native,
        "inventory_digest": inventory["inventory_digest"],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--skip-native-verifiers", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    inventory = load_json(args.inventory)
    result = verify_inventory(
        inventory,
        repo_root=args.repo_root,
        run_native_verifiers=not args.skip_native_verifiers,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
