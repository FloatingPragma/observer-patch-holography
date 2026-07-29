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
    "d1a824a7fba641b2ff3116e566c2ea120023495f2f34ca4340420e8867b547f5"
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
    certificate = loaded["family_band_candidate"][certificate_path]
    artifact = loaded["family_band_candidate"][artifact_path]
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
            "continuum Spin/locality receipt",
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
    require(not hits, f"structured target content found in source parents: {hits}")
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
    _verify_rg_representation_frontier(conditional)
    _verify_scalar_yukawa_source_frontier(conditional)
    _verify_local_ew_order_unit_frontier(conditional)
    _verify_source_clock_frontier(conditional)
    _verify_consumer_schemas(repo_root, inventory["consumer_contract"]["schemas"])
    _verify_firewall(inventory, policy, positive, conditional)

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
    require(inventory["unit_scope"]["emitted_observables"] == [], "frontier emitted poles")
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
    require(
        all(row["status"] in {"open", "partial"} for row in inventory["acceptance_map"]),
        "an issue-#594 acceptance row was overpromoted",
    )
    require(
        all(row["blocking_gates"] for row in inventory["acceptance_map"]),
        "acceptance row lost its blockers",
    )
    acceptance_nine = next(
        row
        for row in inventory["acceptance_map"]
        if row["acceptance_index"] == 9
    )
    require(
        acceptance_nine["blocking_gates"]
        == [
            "event_and_spacetime_action_parent",
            "physical_family_and_matter_pole_attachment",
            "scalar_higgs_and_fj_coordinate",
            "full_yukawa_operator_and_coefficients",
            "common_screen_electroweak_carrier",
            "source_complete_field_census",
        ],
        "acceptance row 9 lost a live native-action attachment",
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
