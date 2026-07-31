"""Conditional seam-character menu and matter-action boundary (#627).

The original v1 certificate silently identified an abstract seam coefficient
group with the pure hypercharge character.  That produced the exact fixed
dimensions 7, 3, and 1, but it was not the diagonal kernel action of the
declared color-weak-hypercharge fixture.  That kernel acts trivially on every
state in the declared local matter table.

This corrected certificate keeps three objects separate:

* the pinned order-six seam branch, whose general A1-A3 exhaustiveness remains
  open after the #624 audit;
* the complete menu of characters of Z2, Z3, and Z6 on the declared
  hypercharge spectrum, including the trivial character;
* the conditional common kernel computed by #567 from the declared current
  and matter tables, which fixes all fifteen states.

No character or 2-representation is promoted to the physical seam action.
Neither is the computed kernel promoted to a source-selected physical global
form.  Those selections and their synchronization with the line/flux
attachment remain named conditional interfaces outside this bounded
classification.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from math import gcd
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

SCHEMA = "oph.seam_grammar_matter_classification_certificate.v2"
MANIFEST_PATH = (
    MODULE_DIR / "manifests" / "seam_grammar_matter_classification_reference.json"
)
SEAM_RECEIPT_PATH = MODULE_DIR / "receipts" / "noncentral_seam_reduction_reference.receipt.json"
MATTER_RECEIPT_PATH = MODULE_DIR / "receipts" / "super_tannakian_matter_reference.receipt.json"
RESPONSE_ARTIFACT_NAME = "charged_response_semantic_artifact.json"
AXIS_RECEIPT_PATH = MODULE_DIR / "receipts" / "axis_center_descent_reference.receipt.json"

ISSUE = 627


def pin_file(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    payload = load_json(path)
    return payload, {
        "path": path.relative_to(MODULE_DIR.parent.parent).as_posix(),
        "sha256": sha256_json(payload),
    }


# ---------------------------------------------------------------------------
# Clause one: the finite order-six branch, with its scope boundary
# ---------------------------------------------------------------------------


def grammar_branch_status(seam_receipt: Mapping[str, Any]) -> dict[str, Any]:
    lift = seam_receipt["complete_coefficient_lift"]["classification"]
    require(
        lift["z6"]["classification"]
        == "identified_with_order_six_column_within_pinned_lanes",
        "GRAMMAR_BRANCH",
        "the pinned receipt must carry its order-six named branch",
    )
    require(
        lift["z6"]["faithful_embeddings_into_seam_class_lane"] == [1, 5],
        "GRAMMAR_BRANCH",
        "the order-six branch must retain its two orientation-related embeddings",
    )
    require(
        lift["z7"]["classification"]
        == "excluded_within_two_pinned_lanes_only"
        and lift["z7"]["homs_into_seam_class_lane"] == [0]
        and lift["z7"]["homs_into_register_lane"] == [0],
        "GRAMMAR_BRANCH",
        "the pinned two-lane Z7 calculation has drifted",
    )
    return {
        "named_branch": "Z6",
        "source": "pinned finite order-six central-column receipt",
        "target_free": True,
        "z7_result_scope": (
            "Z7 has no nontrivial homomorphism into either pinned lane. This "
            "does not prove that the two lanes exhaust every complete A1-A3 "
            "coefficient construction"
        ),
        "general_grammar_classification": "open",
    }


# ---------------------------------------------------------------------------
# Clause three: the realized-module classification, exact
# ---------------------------------------------------------------------------


def hypercharge_character_menu(charge_spectrum: Mapping[str, int]) -> dict[str, Any]:
    """Enumerate every pure-hypercharge character, without selecting one."""

    states = []
    total = 0
    for y_str, multiplicity in sorted(charge_spectrum.items()):
        y = Fraction(y_str)
        q6 = int(6 * y) % 6
        require(6 * y == int(6 * y), "CHARGE_GRID", "every hypercharge must sit on the sixth-integer grid")
        states.append({"hypercharge": y_str, "q6": q6, "multiplicity": int(multiplicity)})
        total += int(multiplicity)
    require(total == 15, "MODULE_SIZE", "the realized module must carry fifteen states")

    table = []
    expected_dimensions = {2: [7, 15], 3: [3, 15], 6: [1, 3, 7, 15]}
    for order in (2, 3, 6):
        characters = []
        for exponent in range(order):
            residues = [
                (exponent * state["q6"]) % order for state in states
            ]
            common = order
            for value in residues:
                common = gcd(common, value)
            fixed = sum(
                state["multiplicity"]
                for state, residue in zip(states, residues)
                if residue == 0
            )
            characters.append(
                {
                    "character_exponent": exponent,
                    "faithful_on_module": common == 1,
                    "fixed_subspace_dimension": fixed,
                }
            )
        dimensions = sorted(
            {row["fixed_subspace_dimension"] for row in characters}
        )
        require(
            dimensions == expected_dimensions[order],
            "CHARACTER_MENU",
            "the exact hypercharge-character fixed-space menu has drifted",
        )
        table.append(
            {
                "group_order": order,
                "characters": characters,
                "fixed_dimension_menu": dimensions,
            }
        )
    require(
        any(
            row["faithful_on_module"] and row["fixed_subspace_dimension"] == 1
            for row in table[-1]["characters"]
        ),
        "CHARACTER_MENU",
        "the order-six pure-hypercharge menu must retain its faithful diagnostic",
    )
    return {
        "interpretation": (
            "conditional pure-U(1)_Y characters; no seam-to-hypercharge "
            "homomorphism is selected"
        ),
        "states": states,
        "groups": table,
    }


def diagonal_kernel_action(
    axis_receipt: Mapping[str, Any],
    matter_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    """Verify the conditional diagonal kernel on the declared matter table."""

    generator = tuple(axis_receipt["kernel_on_realized_tensors"]["cyclic_generator"])
    require(
        generator == (1, 1, 1),
        "DIAGONAL_KERNEL",
        "the declared-table diagonal kernel generator must be (1,1,1)",
    )
    weights = axis_receipt["realized_weight_table"]
    fields = matter_receipt["realized_package"]["fields"]
    matter_labels = ("Q", "u_c", "e_c", "d_c", "L")
    rows = []
    total = 0
    for label in matter_labels:
        weight = weights[label]
        phase = (
            2 * generator[0] * int(weight["triality"])
            + 3 * generator[1] * int(weight["duality"])
            + generator[2] * int(weight["q"])
        ) % 6
        require(
            phase == 0,
            "DIAGONAL_KERNEL",
            "the diagonal quotient kernel must fix every realized matter field",
        )
        dimension = int(fields[label]["dimension"])
        total += dimension
        rows.append(
            {
                "field": label,
                "multiplicity": dimension,
                "phase_sixths": phase,
            }
        )
    require(total == 15, "DIAGONAL_KERNEL", "the realized matter dimension must be fifteen")
    return {
        "generator_color_weak_hypercharge": list(generator),
        "fields": rows,
        "module_dimension": total,
        "faithful_on_module": False,
        "fixed_subspace_dimension": total,
        "subgroup_restrictions": [
            {
                "group_order": order,
                "faithful_on_module": False,
                "fixed_subspace_dimension": total,
            }
            for order in (2, 3, 6)
        ],
        "boundary": (
            "this is the common stabilizer action computed from the declared "
            "current and matter tables. It does not select a physical global "
            "form or identify a seam flux or higher-sector 2-representation "
            "with that kernel"
        ),
    }


def conditional_scope_status(
    axis_receipt: Mapping[str, Any],
    matter_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    """Fail closed unless both upstream packets retain their open source gates."""

    matter_conditional = matter_receipt["conditional_algebraic_gate"]
    matter_physical = matter_receipt["physical_source_gate"]
    axis_conditional = axis_receipt["conditional_algebraic_gate"]
    axis_physical = axis_receipt["physical_global_form_gate"]
    require(
        matter_conditional["passed"] is True
        and axis_conditional["passed"] is True,
        "CONDITIONAL_SCOPE",
        "the conditional matter and axis arithmetic gates must both pass",
    )
    require(
        matter_physical["passed"] is False
        and matter_physical["matter_lift_source_bound"] is False
        and matter_physical["upstream_current_representation_source_bound"] is False,
        "PHYSICAL_SCOPE",
        "the matter packet must retain its open physical-source boundary",
    )
    require(
        axis_physical["passed"] is False
        and axis_physical["upstream_response_physically_source_bound"] is False
        and axis_physical["upstream_matter_physically_source_bound"] is False
        and axis_physical["same_source_loop_to_tensor_kernel_identification"] is False,
        "PHYSICAL_SCOPE",
        "the axis packet must retain its open physical-global-form boundary",
    )
    return {
        "conditional_current_matter_kernel_arithmetic": True,
        "physical_matter_lift_source_bound": False,
        "physical_global_form_source_selected": False,
        "same_source_seam_to_tensor_kernel_identification": False,
        "scope": (
            "declared current, matter, and weight-table fixture; no physical "
            "matter action, global form, or seam-to-kernel identification"
        ),
    }


def mechanism_classification(
    seam_receipt: Mapping[str, Any],
    characters: Mapping[str, Any],
    diagonal_action: Mapping[str, Any],
) -> dict[str, Any]:
    """Classify each row without silently supplying its matter character."""

    rows = seam_receipt["two_type_sector_classification"]["mechanism_table"]
    by_order = {
        row["group_order"]: row for row in characters["groups"]
    }
    classified = []
    for row in rows:
        entry: dict[str, Any] = {
            "module": row["module"],
            "classification": row["classification"],
            "contractible": row["contractible"],
        }
        sector = row.get("sector")
        if (
            not row["contractible"]
            and isinstance(sector, Mapping)
            and "coefficient_order" in sector
        ):
            order = int(sector["coefficient_order"])
            if 6 % order == 0:
                entry["realized_module_action"] = {
                    "status": "supplied_character_required",
                    "pure_hypercharge_character_menu": by_order[order][
                        "characters"
                    ],
                    "diagonal_kernel_restriction": next(
                        item
                        for item in diagonal_action["subgroup_restrictions"]
                        if item["group_order"] == order
                    ),
                    "reading": (
                        "the abstract coefficient group does not choose "
                        "between these actions"
                    ),
                }
            else:
                entry["realized_module_action"] = {
                    "through": "none",
                    "reading": (
                        "the coefficient order lies outside this finite "
                        "comparison menu, so its action is not classified here"
                    ),
                }
        elif row["contractible"]:
            entry["realized_module_action"] = {
                "through": "none",
                "reading": "a contractible mechanism carries no sector to act",
            }
        else:
            entry["realized_module_action"] = {
                "through": "none",
                "reading": (
                    "the pinned row carries no central coefficient order "
                    f"({sector.get('reason', 'no sector data')}), so no "
                    "charge character acts"
                ),
            }
        classified.append(entry)
    return {
        "rows": classified,
        "statement": (
            "each applicable explicitly tested coefficient row carries an "
            "exact character menu. No matter character, flux action, or "
            "2-representation is selected"
        ),
    }


# ---------------------------------------------------------------------------
# Clause two: refinement transport
# ---------------------------------------------------------------------------


def validate_refinement_artifact(artifact: Mapping[str, Any]) -> None:
    require(
        artifact.get("schema") == "oph.charged_response_semantic_artifact.v3",
        "REFINEMENT_SCHEMA",
        "the pinned response artifact schema has drifted",
    )
    body = {
        key: value for key, value in artifact.items() if key != "artifact_sha256"
    }
    require(
        artifact.get("artifact_sha256") == "sha256:" + sha256_json(body),
        "REFINEMENT_SELF_HASH",
        "the pinned response artifact self-hash does not match",
    )
    maps = artifact["physical_refinement_maps"]
    require(
        maps["levels_measured"] == 3
        and maps["per_level_defect_port_count"] == 12
        and maps["equivariant_rotation_count"] == 60,
        "REFINEMENT_METADATA",
        "the measured finite refinement metadata has drifted",
    )
    level_pairs = []
    for payload in maps["port_persistence_maps"]:
        map_body = {
            key: value for key, value in payload.items() if key != "map_hash"
        }
        require(
            payload.get("map_hash") == "sha256:" + sha256_json(map_body),
            "REFINEMENT_MAP_HASH",
            "a stored refinement-map hash does not match",
        )
        require(
            payload["port_map"] == list(range(12)),
            "REFINEMENT_IDENTITY",
            "the pinned persistence maps must be the measured identity maps",
        )
        level_pairs.append((payload["source_level"], payload["target_level"]))
    require(
        level_pairs == [(0, 1), (0, 2)],
        "REFINEMENT_LEVEL_PAIRS",
        "the pinned response artifact must carry the two expected level maps",
    )


def refinement_transport() -> dict[str, Any]:
    artifact, pin = pin_file(MODULE_DIR / "manifests" / RESPONSE_ARTIFACT_NAME)
    validate_refinement_artifact(artifact)
    return {
        "pin": pin,
        "classification": (
            "the selected finite response artifact carries hash-verified "
            "identity port persistence from level zero to levels one and "
            "two. This proves port-set persistence for that artifact, not "
            "general seam-character or physical action transport"
        ),
    }


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------


def control_charge_mutation() -> dict[str, Any]:
    """A same-size off-grid charge spectrum must fail the character menu."""

    doctored = {"-1/2": 2, "-2/3": 3, "1": 1, "1/3": 3, "1/7": 6}
    try:
        hypercharge_character_menu(doctored)
    except CertificateError as error:
        return {"expected_failure": True, "failed": True, "code": error.code}
    return {"expected_failure": True, "failed": False}


def control_seam_receipt_mutation() -> dict[str, Any]:
    seam_receipt = load_json(SEAM_RECEIPT_PATH)
    doctored = json.loads(json.dumps(seam_receipt))
    z7 = doctored["complete_coefficient_lift"]["classification"]["z7"]
    z7["classification"] = "admitted"
    z7["homs_into_seam_class_lane"] = [1]
    try:
        grammar_branch_status(doctored)
    except CertificateError as error:
        return {"expected_failure": True, "failed": True, "code": error.code}
    return {"expected_failure": True, "failed": False}


def control_refinement_mutation() -> dict[str, Any]:
    artifact = load_json(MODULE_DIR / "manifests" / RESPONSE_ARTIFACT_NAME)
    doctored = json.loads(json.dumps(artifact))
    row = doctored["physical_refinement_maps"]["port_persistence_maps"][0]
    row["port_map"][0], row["port_map"][1] = row["port_map"][1], row["port_map"][0]
    map_body = {key: value for key, value in row.items() if key != "map_hash"}
    row["map_hash"] = "sha256:" + sha256_json(map_body)
    artifact_body = {
        key: value for key, value in doctored.items() if key != "artifact_sha256"
    }
    doctored["artifact_sha256"] = "sha256:" + sha256_json(artifact_body)
    try:
        validate_refinement_artifact(doctored)
    except CertificateError as error:
        return {"expected_failure": True, "failed": True, "code": error.code}
    return {"expected_failure": True, "failed": False}


def control_action_conflation(
    characters: Mapping[str, Any], diagonal: Mapping[str, Any]
) -> dict[str, Any]:
    order_six = next(
        row for row in characters["groups"] if row["group_order"] == 6
    )
    faithful = next(
        row
        for row in order_six["characters"]
        if row["faithful_on_module"]
        and row["fixed_subspace_dimension"] == 1
    )
    try:
        require(
            faithful["fixed_subspace_dimension"]
            == diagonal["fixed_subspace_dimension"],
            "ACTION_CONFLATION",
            "the pure-hypercharge character is not the diagonal quotient-kernel action",
        )
    except CertificateError:
        return {"expected_failure": True, "failed": True, "code": "ACTION_CONFLATION"}
    return {"expected_failure": True, "failed": False}


def validate_unselected_actions(mechanisms: Mapping[str, Any]) -> None:
    for row in mechanisms["rows"]:
        action = row["realized_module_action"]
        require(
            "selected_character_exponent" not in action,
            "SELECTION_PROMOTION",
            "no matter character may be selected without the open physical interface",
        )


def control_selection_promotion(mechanisms: Mapping[str, Any]) -> dict[str, Any]:
    """A doctored character promotion must fail the open interface."""

    doctored = json.loads(json.dumps(mechanisms))
    target = next(
        row
        for row in doctored["rows"]
        if row["module"] == "Z6 -> 1"
    )
    target["realized_module_action"]["selected_character_exponent"] = 1
    try:
        validate_unselected_actions(doctored)
    except CertificateError as error:
        return {
            "expected_failure": True,
            "failed": True,
            "code": error.code,
        }
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_payload() -> dict[str, Any]:
    seam_receipt, seam_pin = pin_file(SEAM_RECEIPT_PATH)
    matter_receipt, matter_pin = pin_file(MATTER_RECEIPT_PATH)
    axis_receipt, axis_pin = pin_file(AXIS_RECEIPT_PATH)

    scope = conditional_scope_status(axis_receipt, matter_receipt)
    branch = grammar_branch_status(seam_receipt)
    characters = hypercharge_character_menu(
        matter_receipt["realized_package"]["charge_spectrum"]
    )
    diagonal = diagonal_kernel_action(axis_receipt, matter_receipt)
    mechanisms = mechanism_classification(
        seam_receipt, characters, diagonal
    )
    validate_unselected_actions(mechanisms)
    transport = refinement_transport()

    boundary = seam_receipt["two_type_sector_classification"]["selection_boundary"]
    require(
        "627" in str(boundary),
        "SELECTION_BOUNDARY",
        "the pinned selection boundary must name this issue",
    )

    controls = {
        "charge_mutation": control_charge_mutation(),
        "seam_receipt_mutation": control_seam_receipt_mutation(),
        "refinement_mutation": control_refinement_mutation(),
        "action_conflation": control_action_conflation(characters, diagonal),
        "selection_promotion": control_selection_promotion(mechanisms),
    }
    for name, verdict in controls.items():
        require(
            verdict["expected_failure"] is True and verdict["failed"] is True,
            "CONTROL_NOT_FAILED",
            f"control {name} did not record its required failure",
        )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "claim_boundary": (
            "The pinned finite order-six seam branch, the complete "
            "pure-hypercharge character menus, and the conditional common "
            "kernel action of the declared current and matter tables are "
            "separately classified. The "
            "7/3/1 fixed-space values belong only to particular nontrivial "
            "pure-hypercharge characters; the declared-table kernel fixes all "
            "fifteen matter states. No physical global form, seam character, "
            "flux action, or 2-representation is selected, and general "
            "seam-grammar exhaustiveness remains open."
        ),
        "upstream_pins": {
            "seam_classification_receipt": seam_pin,
            "matter_receipt": matter_pin,
            "diagonal_global_form_receipt": axis_pin,
        },
        "conditional_scope": scope,
        "grammar_branch": branch,
        "hypercharge_character_menu": characters,
        "diagonal_kernel_action": diagonal,
        "mechanism_classification": mechanisms,
        "refinement_transport": transport,
        "matter_action_interface": {
            "id": "physical_sector_mechanism_selection",
            "class": "conditional_open_interface",
            "owner_issue": 569,
            "statement": (
                "selection requires a source-derived seam mechanism, a "
                "character or 2-representation, its relation to the diagonal "
                "table kernel and line/flux data, a source-selected physical "
                "global form, refinement transport, and the physical "
                "family-chain receipts"
            ),
        },
        "controls": controls,
        "bounded_exit": "exact_named_character_and_diagonal_action_classification",
    }
    return payload


def build_manifest() -> dict[str, Any]:
    payload = build_payload()
    manifest = dict(payload)
    manifest["manifest_sha256"] = "sha256:" + sha256_json(payload)
    return manifest


def verify_stored() -> dict[str, Any]:
    stored = load_json(MANIFEST_PATH)
    body = {key: value for key, value in stored.items() if key != "manifest_sha256"}
    require(
        stored.get("manifest_sha256") == "sha256:" + sha256_json(body),
        "MANIFEST_HASH",
        "stored manifest hash does not match its body",
    )
    require(
        body == build_payload(),
        "MANIFEST_DRIFT",
        "stored manifest does not match a deterministic rebuild",
    )
    return {"status": "PASS", "manifest": str(MANIFEST_PATH)}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Seam grammar and matter classification certificate for issue #627")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--output", type=Path, default=MANIFEST_PATH)
    args = parser.parse_args(argv)
    if args.verify:
        print(json.dumps(verify_stored(), indent=2))
        return 0
    manifest = build_manifest()
    write_json(args.output, manifest)
    print(json.dumps({"status": "WROTE", "manifest": str(args.output), "manifest_sha256": manifest["manifest_sha256"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
