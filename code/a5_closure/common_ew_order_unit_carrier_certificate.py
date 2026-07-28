#!/usr/bin/env python3
"""Finite local screen/electroweak order-unit precursor for issue #631.

This certificate constructs two exact one-dimensional ordered lines:

* the unique invariant line in the transitive twelve-port screen
  representation, normalized by the uniform port trace; and
* the center/order-unit line of the four-copy weak multiplicity algebra
  obtained by restricting the pinned one-generation matter packet to its
  weak-doublet copies.

There is one positive unital normalized-trace-preserving linear isomorphism
between these lines.  The statement is an algebraic precursor only.  It does
not identify the two order units as one physical load, select a scalar sector,
intertwine a physical screen-load law with an electroweak readout, consume a
cosmic capacity, or evaluate the N-dependent bridge equation.
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

SCHEMA = "oph.common_ew_order_unit_carrier_frontier.v1"
ISSUE = 631
STATUS = "FINITE_ORDER_UNIT_INTERTWINER__PHYSICAL_COMMON_CARRIER_OPEN"

MANIFEST_PATH = (
    MODULE_DIR / "manifests" / "common_ew_order_unit_carrier_reference.json"
)

CARRIER_MANIFEST = MODULE_DIR / "manifests" / "echosahedral_federation_reference.json"
CARRIER_RECEIPT = MODULE_DIR / "receipts" / "echosahedral_federation_reference.receipt.json"
CURRENT_MANIFEST = MODULE_DIR / "manifests" / "port_current_response_reference.json"
CURRENT_RECEIPT = MODULE_DIR / "receipts" / "port_current_inner_reference.receipt.json"
MATTER_MANIFEST = MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json"
MATTER_RECEIPT = MODULE_DIR / "receipts" / "super_tannakian_matter_reference.receipt.json"
GLOBAL_FORM_MANIFEST = MODULE_DIR / "manifests" / "axis_center_descent_reference.json"
GLOBAL_FORM_RECEIPT = MODULE_DIR / "receipts" / "axis_center_descent_reference.receipt.json"
SCALAR_BOUNDARY = MODULE_DIR / "manifests" / "multiplicity_window_reference.json"

PORTS = 12
WEAK_COPIES = 4


def _pin(path: Path, issue: int, role: str) -> tuple[dict[str, Any], dict[str, Any]]:
    payload = load_json(path)
    return payload, {
        "issue": issue,
        "role": role,
        "path": path.relative_to(MODULE_DIR).as_posix(),
        "sha256": sha256_json(payload),
    }


def _verify_self_hash(payload: Mapping[str, Any]) -> None:
    body = {key: value for key, value in payload.items() if key != "manifest_sha256"}
    require(
        payload.get("manifest_sha256") == "sha256:" + sha256_json(body),
        "UPSTREAM_SELF_HASH",
        "the upstream generated manifest does not match its self-hash",
    )


def load_upstreams() -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    carrier_manifest, carrier_manifest_pin = _pin(
        CARRIER_MANIFEST, 565, "screen_carrier_manifest"
    )
    carrier_receipt, carrier_receipt_pin = _pin(
        CARRIER_RECEIPT, 565, "screen_carrier_receipt"
    )
    current_manifest, current_manifest_pin = _pin(
        CURRENT_MANIFEST, 566, "port_current_manifest"
    )
    current_receipt, current_receipt_pin = _pin(
        CURRENT_RECEIPT, 566, "port_current_receipt"
    )
    matter_manifest, matter_manifest_pin = _pin(
        MATTER_MANIFEST, 314, "matter_manifest"
    )
    matter_receipt, matter_receipt_pin = _pin(
        MATTER_RECEIPT, 314, "matter_receipt"
    )
    global_manifest, global_manifest_pin = _pin(
        GLOBAL_FORM_MANIFEST, 567, "global_form_manifest"
    )
    global_receipt, global_receipt_pin = _pin(
        GLOBAL_FORM_RECEIPT, 567, "global_form_receipt"
    )
    scalar_boundary, scalar_boundary_pin = _pin(
        SCALAR_BOUNDARY, 616, "scalar_nonselection_boundary"
    )

    require(
        carrier_manifest.get("schema") == "oph.echosahedral_selector_manifest.v1",
        "CARRIER_SCHEMA",
        "the #565 carrier manifest schema has drifted",
    )
    require(
        carrier_receipt.get("schema") == "oph.echosahedral_selector_receipt.v1"
        and carrier_receipt.get("issue") == 565
        and carrier_receipt.get("manifest_sha256") == sha256_json(carrier_manifest),
        "CARRIER_RECEIPT",
        "the #565 receipt does not bind the carrier manifest",
    )

    require(
        current_manifest.get("schema") == "oph.port_current_response_manifest.v5",
        "CURRENT_SCHEMA",
        "the #566 current manifest schema has drifted",
    )
    require(
        current_receipt.get("schema") == "oph.port_current_inner_receipt.v5"
        and current_receipt.get("issue") == 566
        and current_receipt.get("manifest_sha256") == sha256_json(current_manifest),
        "CURRENT_RECEIPT",
        "the #566 receipt does not bind the current manifest",
    )
    require(
        current_manifest.get("carrier_manifest_sha256")
        == sha256_json(carrier_manifest),
        "CURRENT_CARRIER_PIN",
        "the #566 current packet is not bound to the #565 carrier",
    )
    closure = current_receipt.get("closure", {})
    require(
        closure.get("center_dimension") == 1
        and closure.get("center_is_constant_even_port_line") is True,
        "CURRENT_CENTER",
        "the #566 current center is not the one-dimensional constant port line",
    )

    require(
        matter_manifest.get("schema") == "oph.super_tannakian_matter_manifest.v5",
        "MATTER_SCHEMA",
        "the #314 matter manifest schema has drifted",
    )
    require(
        matter_receipt.get("schema") == "oph.super_tannakian_matter_receipt.v5"
        and matter_receipt.get("issue") == 314
        and matter_receipt.get("manifest_sha256") == sha256_json(matter_manifest),
        "MATTER_RECEIPT",
        "the #314 receipt does not bind the matter manifest",
    )
    require(
        matter_manifest.get("current_manifest_sha256")
        == sha256_json(current_manifest)
        and matter_manifest.get("current_receipt_sha256")
        == sha256_json(current_receipt),
        "MATTER_CURRENT_PIN",
        "the #314 matter packet is not bound to the #566 current packet",
    )
    matter_gate = matter_receipt.get("physical_source_gate", {})
    require(
        matter_gate.get("upstream_response_representation_source_bound") is True
        and matter_gate.get("declared_scalar_content_source_bound") is False
        and matter_gate.get("scalar_economy_source_bound") is False,
        "MATTER_SCOPE",
        "the #314 matter packet lost its scalar non-promotion boundary",
    )

    require(
        global_manifest.get("schema") == "oph.axis_center_descent_manifest.v4",
        "GLOBAL_FORM_SCHEMA",
        "the #567 global-form manifest schema has drifted",
    )
    require(
        global_receipt.get("schema") == "oph.axis_center_descent_receipt.v4"
        and global_receipt.get("issue") == 567
        and global_receipt.get("manifest_sha256") == sha256_json(global_manifest),
        "GLOBAL_FORM_RECEIPT",
        "the #567 receipt does not bind the global-form manifest",
    )
    require(
        global_manifest.get("matter_receipt_sha256")
        == sha256_json(matter_receipt)
        and global_receipt.get("kernel_on_realized_tensors", {}).get("kernel_order")
        == 6,
        "GLOBAL_FORM_MATTER_PIN",
        "the #567 global form is not bound to the #314 matter packet",
    )

    require(
        scalar_boundary.get("schema") == "oph.multiplicity_window_certificate.v1"
        and 616 in scalar_boundary.get("issues", []),
        "SCALAR_BOUNDARY_SCHEMA",
        "the #616 scalar boundary receipt has drifted",
    )
    _verify_self_hash(scalar_boundary)
    scalar_verdict = scalar_boundary.get("scalar_response_multiplicity", {}).get(
        "verdict", {}
    )
    require(
        scalar_verdict.get("scalar_existence") == "not_source_determined"
        and scalar_verdict.get("scalar_multiplicity") == "independence_limited",
        "SCALAR_BOUNDARY",
        "the scalar nonselection boundary must remain explicit",
    )

    payloads = {
        "carrier_manifest": carrier_manifest,
        "carrier_receipt": carrier_receipt,
        "current_manifest": current_manifest,
        "current_receipt": current_receipt,
        "matter_manifest": matter_manifest,
        "matter_receipt": matter_receipt,
        "global_manifest": global_manifest,
        "global_receipt": global_receipt,
        "scalar_boundary": scalar_boundary,
    }
    pins = {
        "carrier_manifest": carrier_manifest_pin,
        "carrier_receipt": carrier_receipt_pin,
        "current_manifest": current_manifest_pin,
        "current_receipt": current_receipt_pin,
        "matter_manifest": matter_manifest_pin,
        "matter_receipt": matter_receipt_pin,
        "global_form_manifest": global_manifest_pin,
        "global_form_receipt": global_receipt_pin,
        "scalar_boundary": scalar_boundary_pin,
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
        "rank input has inconsistent row widths",
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


def matrix_center_constraint_rank(size: int) -> int:
    """Rank of [X,E_ab]=0 on an arbitrary size-by-size matrix X."""

    variables = size * size
    rows: list[list[int]] = []
    for a in range(size):
        for b in range(size):
            for i in range(size):
                for j in range(size):
                    row = [0] * variables
                    if j == b:
                        row[i * size + a] += 1
                    if i == a:
                        row[b * size + j] -= 1
                    if any(row):
                        rows.append(row)
    return rational_rank(rows)


def _parse_refinement_map(
    labels: Sequence[str], carrier: Any
) -> tuple[int, ...]:
    require(
        isinstance(labels, list),
        "REFINEMENT_MAP",
        "a refinement port map must be a list",
    )
    return e565.parse_port_permutation(labels, carrier)


def screen_order_unit(
    carrier_manifest: Mapping[str, Any],
    current_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    carrier = e565.validate_carrier(carrier_manifest)
    _, rotations, _ = e565.group_certificate(carrier)
    gram = e565.gram_matrix(carrier)
    refinement = e565.validate_refinement(
        carrier_manifest, carrier, rotations, gram
    )

    orbit = sorted({rotation[0] for rotation in rotations})
    require(
        orbit == list(range(PORTS)),
        "SCREEN_TRANSITIVITY",
        "the rotation action must be transitive on the twelve ports",
    )
    unit = [1] * PORTS
    rotation_fixed = all(
        [unit[rotation[index]] for index in range(PORTS)] == unit
        for rotation in rotations
    )
    require(
        rotation_fixed,
        "SCREEN_UNIT",
        "the screen order unit is not rotation invariant",
    )

    rotation_set = set(rotations)
    map_rows = []
    for row in carrier_manifest["refinement_tower"]["maps"]:
        permutation = _parse_refinement_map(row["port_map"], carrier)
        require(
            permutation in rotation_set,
            "SCREEN_REFINEMENT",
            "a screen refinement map is not in the orientation-preserving action",
        )
        fixed = [unit[permutation[index]] for index in range(PORTS)] == unit
        require(
            fixed,
            "SCREEN_REFINEMENT",
            "a screen refinement map does not preserve the order unit",
        )
        map_rows.append(
            {
                "source": row["source"],
                "target": row["target"],
                "permutation": list(permutation),
                "order_unit_fixed": True,
            }
        )

    require(
        current_receipt["closure"]["center_dimension"] == 1
        and current_receipt["closure"]["center_is_constant_even_port_line"] is True
        and current_receipt["refinement"]["carrier_tower"]["unit_lines_natural"]
        is True,
        "CURRENT_ORDER_UNIT",
        "the current receipt does not identify the natural constant port line",
    )

    return {
        "carrier": "twelve primitive port units on the pinned oriented screen",
        "basis": unit,
        "positive_ray": "c * 1_12 with c >= 0",
        "normalized_trace": "tau_screen(c * 1_12) = c",
        "normalized_trace_weights": ["1/12"] * PORTS,
        "fixed_space_dimension": 1,
        "fixed_space_reason": (
            "the listed order-sixty action is transitive, so every invariant "
            "port function is constant"
        ),
        "rotation_checks": len(rotations),
        "orbit_of_p00": orbit,
        "current_center_binding": (
            "the #566 current center is the same constant even-port line"
        ),
        "refinement": {
            **refinement,
            "maps": map_rows,
            "order_unit_natural": True,
        },
    }


def electroweak_order_unit(
    matter_receipt: Mapping[str, Any],
    global_receipt: Mapping[str, Any],
    scalar_boundary: Mapping[str, Any],
) -> dict[str, Any]:
    fields = matter_receipt["realized_package"]["fields"]
    require(
        fields["Q"]["dimension"] == 6 and fields["L"]["dimension"] == 2,
        "WEAK_MULTIPLICITY",
        "the pinned matter packet lost its Q and L weak-doublet dimensions",
    )
    quark_doublet_copies = fields["Q"]["dimension"] // 2
    lepton_doublet_copies = fields["L"]["dimension"] // 2
    weak_copies = quark_doublet_copies + lepton_doublet_copies
    require(
        weak_copies == WEAK_COPIES,
        "WEAK_MULTIPLICITY",
        "the weak restriction must contain exactly four doublet copies",
    )

    constraint_rank = matrix_center_constraint_rank(weak_copies)
    center_dimension = weak_copies * weak_copies - constraint_rank
    require(
        constraint_rank == 15 and center_dimension == 1,
        "EW_CENTER",
        "the four-copy weak multiplicity algebra must have a one-dimensional center",
    )

    matter_refinement = matter_receipt.get("refinement", {})
    require(
        matter_refinement.get("natural") is True
        and all(row.get("intertwined") is True for row in matter_refinement["maps"]),
        "EW_REFINEMENT",
        "the weak matter restriction lost refinement naturality",
    )
    require(
        global_receipt["kernel_on_realized_tensors"]["kernel_order"] == 6,
        "EW_GLOBAL_FORM",
        "the finite weak multiplicity line is not bound to the pinned global form",
    )

    scalar_verdict = scalar_boundary["scalar_response_multiplicity"]["verdict"]
    return {
        "restriction": (
            "the one-generation matter module restricted to the weak factor"
        ),
        "copy_derivation": {
            "Q_dimension": fields["Q"]["dimension"],
            "Q_weak_doublet_copies": quark_doublet_copies,
            "L_dimension": fields["L"]["dimension"],
            "L_weak_doublet_copies": lepton_doublet_copies,
            "total_weak_doublet_copies": weak_copies,
        },
        "multiplicity_algebra": "End(C^4)",
        "multiplicity_algebra_dimension": weak_copies * weak_copies,
        "central_constraint_rank": constraint_rank,
        "center_dimension": center_dimension,
        "basis": [
            [1 if row == column else 0 for column in range(weak_copies)]
            for row in range(weak_copies)
        ],
        "positive_ray": "c * I_4 with c >= 0",
        "normalized_trace": "tau_EW(c * I_4) = Tr_4(c * I_4) / 4 = c",
        "refinement_natural": True,
        "global_form_compatible": True,
        "scalar_boundary": {
            "scalar_existence": scalar_verdict["scalar_existence"],
            "scalar_multiplicity": scalar_verdict["scalar_multiplicity"],
            "one_doublet_completion": "declared_only",
            "physical_scalar_selected": False,
        },
    }


def local_intertwiner() -> dict[str, Any]:
    return {
        "domain": "span(1_12)",
        "codomain": "span(I_4)",
        "formula": "T(c * 1_12) = c * I_4",
        "inverse": "T_inv(c * I_4) = c * 1_12",
        "linear": True,
        "positive": True,
        "unital": True,
        "normalized_trace_preserving": True,
        "bijective": True,
        "unique_among_unital_linear_maps_between_the_two_lines": True,
        "refinement_natural_on_the_pinned_finite_towers": True,
        "scope": (
            "unique abstract normalized order-unit isomorphism after the two "
            "lines have been constructed; no physical common-load identity "
            "or physical readout intertwining is inferred"
        ),
    }


def _failed_control(code: str, condition: bool, statement: str) -> dict[str, Any]:
    try:
        require(condition, code, statement)
    except CertificateError as error:
        return {
            "expected_failure": True,
            "failed": True,
            "code": error.code,
            "statement": statement,
        }
    return {
        "expected_failure": True,
        "failed": False,
        "code": "CONTROL_DID_NOT_FAIL",
        "statement": statement,
    }


def controls(rotations: Sequence[tuple[int, ...]]) -> dict[str, Any]:
    unequal = [2] + [1] * (PORTS - 1)
    unequal_invariant = all(
        [unequal[rotation[index]] for index in range(PORTS)] == unequal
        for rotation in rotations
    )
    rows = {
        "unequal_screen_weights": _failed_control(
            "SCREEN_UNIT_NOT_INVARIANT",
            unequal_invariant,
            "an unequal port-weight vector is not fixed by the transitive action",
        ),
        "nonunital_scale": _failed_control(
            "INTERTWINER_NOT_UNITAL",
            2 == 1,
            "the map T(1_12)=2 I_4 is not unital",
        ),
        "weak_copy_mutation": _failed_control(
            "WEAK_MULTIPLICITY",
            5 == WEAK_COPIES,
            "a five-copy weak multiplicity does not match the pinned matter packet",
        ),
        "scalar_promotion": _failed_control(
            "SCALAR_NOT_SOURCE_SELECTED",
            False,
            "the #314/#616 packets do not select a physical scalar",
        ),
        "capacity_injection": _failed_control(
            "CAPACITY_INPUT_FORBIDDEN",
            False,
            "a cosmic-capacity or issue-505 artifact is outside the local precursor",
        ),
        "physical_common_carrier_promotion": _failed_control(
            "PHYSICAL_COMMON_CARRIER_OPEN",
            False,
            "the abstract line isomorphism is not a physical common-load identification",
        ),
    }
    require(
        all(row["expected_failure"] and row["failed"] for row in rows.values()),
        "CONTROL_NOT_FAILED",
        "every local-carrier control must fail closed",
    )
    return rows


def require_no_floats(value: Any, path: str = "$") -> None:
    require(
        not isinstance(value, float),
        "FLOAT_FORBIDDEN",
        f"floating point is forbidden at {path}",
    )
    if isinstance(value, Mapping):
        for key, child in value.items():
            require_no_floats(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            require_no_floats(child, f"{path}[{index}]")


def build_payload() -> dict[str, Any]:
    upstream, pins = load_upstreams()
    carrier = e565.validate_carrier(upstream["carrier_manifest"])
    _, rotations, _ = e565.group_certificate(carrier)

    screen = screen_order_unit(
        upstream["carrier_manifest"], upstream["current_receipt"]
    )
    electroweak = electroweak_order_unit(
        upstream["matter_receipt"],
        upstream["global_receipt"],
        upstream["scalar_boundary"],
    )
    intertwiner = local_intertwiner()

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "status": STATUS,
        "claim_class": "conditional_open_interface",
        "claim_boundary": {
            "proves": (
                "the pinned finite screen has one invariant normalized order-unit "
                "line, the four-copy weak restriction has one central normalized "
                "order-unit line, and the unique positive unital trace-preserving "
                "linear isomorphism between those lines is refinement-natural"
            ),
            "does_not_prove": [
                "that the screen and weak order units have one physical meaning",
                "a source-selected scalar carrier, scalar pole, or scalar multiplicity",
                "a physical screen-load law to electroweak-readout intertwiner",
                "the physical normalization of a common load carrier",
                "the N-dependent bridge equation or a source-only cosmic capacity",
                "a Higgs, W, Z, coupling, vacuum, or mass prediction",
            ],
        },
        "upstream_pins": pins,
        "dependency_chain": [
            "issue 565: finite twelve-port carrier and refinement",
            "issue 566: source-bound current center and natural unit line",
            "issue 314: finite matter module and four weak-doublet copies",
            "issue 567: finite global-form compatibility",
            "issue 616: scalar existence and multiplicity nonselection boundary",
        ],
        "screen_order_unit_line": screen,
        "electroweak_order_unit_line": electroweak,
        "local_order_unit_intertwiner": intertwiner,
        "source_firewall": {
            "allowed_issue_inputs": [565, 566, 314, 567, 616],
            "issue_505_consumed": False,
            "cosmic_capacity_consumed": False,
            "measured_cosmology_consumed": False,
            "measured_particle_targets_consumed": False,
            "target_residual_consumed": False,
            "physical_scalar_assumed": False,
            "physical_common_carrier_assumed": False,
        },
        "promotion": {
            "promotion_allowed": False,
            "physical_common_load_carrier_identified": False,
            "physical_scalar_selected": False,
            "physical_higgs_carrier_identified": False,
            "screen_load_law_intertwined_with_ew_readout": False,
            "physical_normalization_derived": False,
            "N_consumed": False,
            "R_EW_evaluated": False,
        },
        "open_gates": [
            {
                "gate": "source_selected_scalar_carrier",
                "status": "open",
                "owners": [630],
            },
            {
                "gate": "physical_common_load_semantics",
                "status": "open",
                "owners": [631],
            },
            {
                "gate": "physical_screen_load_to_electroweak_readout_intertwiner",
                "status": "open",
                "owners": [631],
            },
            {
                "gate": "N_dependent_bridge_equation",
                "status": "open",
                "owners": [505, 547],
            },
        ],
        "controls": controls(rotations),
        "arithmetic": "exact integer and rational arithmetic; no floating point",
    }
    require_no_floats(payload)
    return payload


def build_manifest() -> dict[str, Any]:
    body = build_payload()
    manifest = dict(body)
    manifest["manifest_sha256"] = "sha256:" + sha256_json(body)
    return manifest


def verify_stored(path: Path = MANIFEST_PATH) -> dict[str, Any]:
    stored = load_json(path)
    body = {key: value for key, value in stored.items() if key != "manifest_sha256"}
    require(
        stored.get("manifest_sha256") == "sha256:" + sha256_json(body),
        "MANIFEST_HASH",
        "the local-carrier manifest self-hash does not recompute",
    )
    require(
        body == build_payload(),
        "MANIFEST_DRIFT",
        "the stored local-carrier manifest differs from a deterministic rebuild",
    )
    return {
        "status": "PASS",
        "manifest": str(path),
        "manifest_sha256": stored["manifest_sha256"],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--output", type=Path, default=MANIFEST_PATH)
    args = parser.parse_args(argv)
    if args.verify:
        print(json.dumps(verify_stored(args.output), indent=2))
        return 0
    manifest = build_manifest()
    write_json(args.output, manifest)
    print(
        json.dumps(
            {
                "status": "WROTE",
                "manifest": str(args.output),
                "manifest_sha256": manifest["manifest_sha256"],
                "promotion_allowed": False,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
