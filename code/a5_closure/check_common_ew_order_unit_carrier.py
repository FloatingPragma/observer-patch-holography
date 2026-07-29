#!/usr/bin/env python3
"""Independent replay for the finite #631 order-unit precursor.

This checker imports no code from the producer.  It resolves the pinned
upstream artifacts, recomputes the finite invariant and center dimensions,
checks every non-promotion guard, and emits a compact deterministic receipt.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

MANIFEST_SCHEMA = "oph.common_ew_order_unit_carrier_frontier.v1"
RECEIPT_SCHEMA = "oph.common_ew_order_unit_carrier_receipt.v1"
ISSUE = 631
STATUS = "FINITE_ORDER_UNIT_INTERTWINER__PHYSICAL_COMMON_CARRIER_OPEN"

MANIFEST_PATH = (
    MODULE_DIR / "manifests" / "common_ew_order_unit_carrier_reference.json"
)
RECEIPT_PATH = (
    MODULE_DIR / "receipts" / "common_ew_order_unit_carrier_reference.receipt.json"
)

UPSTREAM_FILES = {
    "carrier_manifest": (
        565,
        "screen_carrier_manifest",
        MODULE_DIR / "manifests" / "echosahedral_federation_reference.json",
    ),
    "carrier_receipt": (
        565,
        "screen_carrier_receipt",
        MODULE_DIR
        / "receipts"
        / "echosahedral_federation_reference.receipt.json",
    ),
    "current_manifest": (
        566,
        "port_current_manifest",
        MODULE_DIR / "manifests" / "port_current_response_reference.json",
    ),
    "current_receipt": (
        566,
        "port_current_receipt",
        MODULE_DIR / "receipts" / "port_current_inner_reference.receipt.json",
    ),
    "matter_manifest": (
        314,
        "matter_manifest",
        MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json",
    ),
    "matter_receipt": (
        314,
        "matter_receipt",
        MODULE_DIR / "receipts" / "super_tannakian_matter_reference.receipt.json",
    ),
    "global_form_manifest": (
        567,
        "global_form_manifest",
        MODULE_DIR / "manifests" / "axis_center_descent_reference.json",
    ),
    "global_form_receipt": (
        567,
        "global_form_receipt",
        MODULE_DIR / "receipts" / "axis_center_descent_reference.receipt.json",
    ),
    "scalar_boundary": (
        616,
        "scalar_nonselection_boundary",
        MODULE_DIR / "manifests" / "multiplicity_window_reference.json",
    ),
}


def _pin(issue: int, role: str, path: Path, payload: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "issue": issue,
        "role": role,
        "path": path.relative_to(MODULE_DIR).as_posix(),
        "sha256": sha256_json(payload),
    }


def load_upstreams() -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    payloads = {
        key: load_json(path)
        for key, (_, _, path) in UPSTREAM_FILES.items()
    }
    pins = {
        key: _pin(issue, role, path, payloads[key])
        for key, (issue, role, path) in UPSTREAM_FILES.items()
    }
    return payloads, pins


def rational_rank(rows: Sequence[Sequence[int | Fraction]]) -> int:
    matrix = [[Fraction(value) for value in row] for row in rows]
    if not matrix:
        return 0
    width = len(matrix[0])
    require(
        all(len(row) == width for row in matrix),
        "RANK_MATRIX",
        "independent rank input has inconsistent row widths",
    )
    rank = 0
    column = 0
    while rank < len(matrix) and column < width:
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column] != 0),
            None,
        )
        if pivot is None:
            column += 1
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = matrix[rank][column]
        matrix[rank] = [entry / scale for entry in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                left - factor * right
                for left, right in zip(matrix[row], matrix[rank])
            ]
        rank += 1
        column += 1
    return rank


def center_constraint_rank(size: int) -> int:
    rows: list[list[int]] = []
    for a in range(size):
        for b in range(size):
            for i in range(size):
                for j in range(size):
                    row = [0] * (size * size)
                    if j == b:
                        row[i * size + a] += 1
                    if i == a:
                        row[b * size + j] -= 1
                    if any(row):
                        rows.append(row)
    return rational_rank(rows)


def verify_upstream_chain(
    upstream: Mapping[str, Mapping[str, Any]],
) -> None:
    carrier_manifest = upstream["carrier_manifest"]
    carrier_receipt = upstream["carrier_receipt"]
    current_manifest = upstream["current_manifest"]
    current_receipt = upstream["current_receipt"]
    matter_manifest = upstream["matter_manifest"]
    matter_receipt = upstream["matter_receipt"]
    global_manifest = upstream["global_form_manifest"]
    global_receipt = upstream["global_form_receipt"]
    scalar_boundary = upstream["scalar_boundary"]

    require(
        carrier_receipt.get("issue") == 565
        and carrier_receipt.get("schema") == "oph.echosahedral_selector_receipt.v1"
        and carrier_receipt.get("manifest_sha256") == sha256_json(carrier_manifest),
        "CARRIER_CHAIN",
        "the carrier manifest/receipt pair does not resolve",
    )
    require(
        current_receipt.get("issue") == 566
        and current_receipt.get("schema") == "oph.port_current_inner_receipt.v5"
        and current_receipt.get("manifest_sha256") == sha256_json(current_manifest)
        and current_manifest.get("carrier_manifest_sha256")
        == sha256_json(carrier_manifest),
        "CURRENT_CHAIN",
        "the current packet does not resolve to the carrier packet",
    )
    require(
        matter_receipt.get("issue") == 314
        and matter_receipt.get("schema")
        == "oph.super_tannakian_matter_receipt.v5"
        and matter_receipt.get("manifest_sha256") == sha256_json(matter_manifest)
        and matter_manifest.get("current_manifest_sha256")
        == sha256_json(current_manifest)
        and matter_manifest.get("current_receipt_sha256")
        == sha256_json(current_receipt),
        "MATTER_CHAIN",
        "the matter packet does not resolve to the current packet",
    )
    require(
        global_receipt.get("issue") == 567
        and global_receipt.get("schema") == "oph.axis_center_descent_receipt.v4"
        and global_receipt.get("manifest_sha256") == sha256_json(global_manifest)
        and global_manifest.get("matter_receipt_sha256")
        == sha256_json(matter_receipt)
        and global_receipt["kernel_on_realized_tensors"]["kernel_order"] == 6,
        "GLOBAL_FORM_CHAIN",
        "the global-form packet does not resolve to the matter packet",
    )

    scalar_body = {
        key: value
        for key, value in scalar_boundary.items()
        if key != "manifest_sha256"
    }
    scalar_verdict = scalar_boundary.get("scalar_response_multiplicity", {}).get(
        "verdict", {}
    )
    require(
        scalar_boundary.get("manifest_sha256")
        == "sha256:" + sha256_json(scalar_body)
        and scalar_boundary.get("schema")
        == "oph.multiplicity_window_certificate.v1"
        and 616 in scalar_boundary.get("issues", [])
        and scalar_verdict.get("scalar_existence") == "not_source_determined"
        and scalar_verdict.get("scalar_multiplicity") == "independence_limited",
        "SCALAR_CHAIN",
        "the scalar nonselection boundary does not resolve",
    )


def verify_manifest(manifest: Mapping[str, Any]) -> dict[str, Any]:
    require(
        manifest.get("schema") == MANIFEST_SCHEMA
        and manifest.get("issue") == ISSUE
        and manifest.get("status") == STATUS
        and manifest.get("claim_class") == "conditional_open_interface",
        "MANIFEST_HEADER",
        "the local-carrier manifest header has drifted",
    )
    body = {key: value for key, value in manifest.items() if key != "manifest_sha256"}
    require(
        manifest.get("manifest_sha256") == "sha256:" + sha256_json(body),
        "MANIFEST_HASH",
        "the local-carrier manifest self-hash does not recompute",
    )

    upstream, pins = load_upstreams()
    verify_upstream_chain(upstream)
    require(
        manifest.get("upstream_pins") == pins,
        "UPSTREAM_PINS",
        "the local-carrier manifest does not contain the exact closed upstream pins",
    )

    carrier = e565.validate_carrier(upstream["carrier_manifest"])
    _, rotations, _ = e565.group_certificate(carrier)
    gram = e565.gram_matrix(carrier)
    expected_refinement = e565.validate_refinement(
        upstream["carrier_manifest"], carrier, rotations, gram
    )
    orbit = sorted({rotation[0] for rotation in rotations})
    require(
        orbit == list(range(12)),
        "SCREEN_TRANSITIVITY",
        "the independent carrier replay is not transitive",
    )
    screen = manifest.get("screen_order_unit_line", {})
    require(
        screen.get("basis") == [1] * 12
        and screen.get("fixed_space_dimension") == 1
        and screen.get("rotation_checks") == 60
        and screen.get("orbit_of_p00") == orbit
        and screen.get("normalized_trace_weights") == ["1/12"] * 12,
        "SCREEN_LINE",
        "the screen invariant line does not match the independent replay",
    )

    rotation_set = set(rotations)
    expected_maps = []
    for row in upstream["carrier_manifest"]["refinement_tower"]["maps"]:
        permutation = e565.parse_port_permutation(row["port_map"], carrier)
        require(
            permutation in rotation_set,
            "SCREEN_REFINEMENT",
            "an upstream refinement map is not a listed rotation",
        )
        expected_maps.append(
            {
                "source": row["source"],
                "target": row["target"],
                "permutation": list(permutation),
                "order_unit_fixed": True,
            }
        )
    refinement = screen.get("refinement", {})
    require(
        refinement.get("maps") == expected_maps
        and refinement.get("order_unit_natural") is True
        and all(
            refinement.get(key) == value
            for key, value in expected_refinement.items()
        ),
        "SCREEN_REFINEMENT",
        "the manifest refinement replay has drifted",
    )

    matter_receipt = upstream["matter_receipt"]
    matter_gate = matter_receipt["physical_source_gate"]
    require(
        matter_gate["declared_scalar_content_source_bound"] is False
        and matter_gate["scalar_economy_source_bound"] is False,
        "SCALAR_PROMOTION",
        "the matter packet may not be relabeled as a scalar source",
    )
    fields = matter_receipt["realized_package"]["fields"]
    weak_copies = fields["Q"]["dimension"] // 2 + fields["L"]["dimension"] // 2
    constraint_rank = center_constraint_rank(weak_copies)
    center_dimension = weak_copies * weak_copies - constraint_rank
    require(
        weak_copies == 4 and constraint_rank == 15 and center_dimension == 1,
        "EW_CENTER",
        "the independent weak-multiplicity center replay failed",
    )
    electroweak = manifest.get("electroweak_order_unit_line", {})
    require(
        electroweak.get("copy_derivation", {}).get(
            "total_weak_doublet_copies"
        )
        == weak_copies
        and electroweak.get("multiplicity_algebra") == "End(C^4)"
        and electroweak.get("central_constraint_rank") == constraint_rank
        and electroweak.get("center_dimension") == center_dimension
        and electroweak.get("scalar_boundary", {}).get(
            "physical_scalar_selected"
        )
        is False,
        "EW_LINE",
        "the electroweak order-unit line does not match the independent replay",
    )

    intertwiner = manifest.get("local_order_unit_intertwiner", {})
    required_intertwiner = {
        "linear",
        "positive",
        "unital",
        "normalized_trace_preserving",
        "bijective",
        "unique_among_unital_linear_maps_between_the_two_lines",
        "refinement_natural_on_the_pinned_finite_towers",
    }
    require(
        intertwiner.get("domain") == "span(1_12)"
        and intertwiner.get("codomain") == "span(I_4)"
        and all(intertwiner.get(key) is True for key in required_intertwiner)
        and "no physical common-load identity" in intertwiner.get("scope", ""),
        "INTERTWINER",
        "the local line isomorphism or its scope boundary has drifted",
    )

    firewall = manifest.get("source_firewall", {})
    require(
        firewall.get("allowed_issue_inputs") == [565, 566, 314, 567, 616]
        and firewall.get("issue_505_consumed") is False
        and firewall.get("cosmic_capacity_consumed") is False
        and firewall.get("measured_cosmology_consumed") is False
        and firewall.get("measured_particle_targets_consumed") is False
        and firewall.get("target_residual_consumed") is False
        and firewall.get("physical_scalar_assumed") is False
        and firewall.get("physical_common_carrier_assumed") is False,
        "SOURCE_FIREWALL",
        "the local precursor consumed a forbidden source or promoted an open premise",
    )

    promotion = manifest.get("promotion", {})
    require(
        set(promotion)
        == {
            "promotion_allowed",
            "physical_common_load_carrier_identified",
            "physical_scalar_selected",
            "physical_higgs_carrier_identified",
            "screen_load_law_intertwined_with_ew_readout",
            "physical_normalization_derived",
            "N_consumed",
            "R_EW_evaluated",
        }
        and all(value is False for value in promotion.values()),
        "PROMOTION_GUARD",
        "every physical, scalar, capacity, and bridge promotion must remain false",
    )
    require(
        manifest.get("open_gates")
        == [
            {
                "gate": "source_selected_scalar_carrier",
                "status": "open",
                "owners": [636],
            },
            {
                "gate": "physical_common_load_semantics",
                "status": "open",
                "owners": [631],
            },
            {
                "gate": (
                    "physical_screen_load_to_electroweak_readout_intertwiner"
                ),
                "status": "open",
                "owners": [631],
            },
            {
                "gate": "N_dependent_bridge_equation",
                "status": "open",
                "owners": [505, 547],
            },
        ],
        "OPEN_GATES",
        "the physical, scalar, or N-dependent boundary was closed without a source",
    )
    require(
        all(
            row.get("expected_failure") is True
            and row.get("failed") is True
            for row in manifest.get("controls", {}).values()
        )
        and len(manifest.get("controls", {})) == 6,
        "CONTROL_GUARD",
        "the fail-closed control battery is incomplete",
    )
    return {
        "screen_fixed_space_dimension": 1,
        "screen_rotation_checks": len(rotations),
        "screen_refinement_maps": len(expected_maps),
        "weak_doublet_copies": weak_copies,
        "weak_multiplicity_center_constraint_rank": constraint_rank,
        "weak_multiplicity_center_dimension": center_dimension,
        "positive_unital_intertwiner_unique": True,
        "physical_common_carrier_identified": False,
        "physical_scalar_selected": False,
        "issue_505_consumed": False,
        "R_EW_evaluated": False,
    }


def build_receipt(manifest: Mapping[str, Any]) -> dict[str, Any]:
    replay = verify_manifest(manifest)
    receipt = {
        "schema": RECEIPT_SCHEMA,
        "issue": ISSUE,
        "manifest_sha256": sha256_json(manifest),
        "manifest_self_hash": manifest["manifest_sha256"],
        "upstream_bundle_digest": sha256_json(manifest["upstream_pins"]),
        "independent_replay": replay,
        "status": STATUS,
        "claim_class": "conditional_open_interface",
        "promotion_allowed": False,
        "open_boundary": (
            "the finite normalized order-unit lines and their unique abstract "
            "isomorphism are certified; physical common-load semantics, scalar "
            "selection, physical normalization, capacity, and the bridge "
            "equation remain open"
        ),
    }
    return receipt


def verify_receipt(
    manifest: Mapping[str, Any], receipt: Mapping[str, Any]
) -> None:
    require(
        receipt == build_receipt(manifest),
        "RECEIPT_MISMATCH",
        "the stored independent receipt is stale, malformed, or tampered",
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--receipt", type=Path, default=RECEIPT_PATH)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)

    manifest = load_json(args.manifest)
    expected = build_receipt(manifest)
    if args.write:
        write_json(args.receipt, expected)
        status = "WROTE"
    else:
        stored = load_json(args.receipt)
        verify_receipt(manifest, stored)
        status = "PASS"
    print(
        json.dumps(
            {
                "status": status,
                "manifest": str(args.manifest),
                "receipt": str(args.receipt),
                "promotion_allowed": False,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
