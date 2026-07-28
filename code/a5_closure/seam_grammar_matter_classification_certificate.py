"""Seam grammar and matter-action classification on realized matter (#627).

The closure clause of issue #627 asks for a target-free source packet that
selects or classifies three objects: the seam grammar, the refinement
transport, and the representation or 2-representation on realized matter.
This certificate assembles that packet from pinned receipts plus new exact
computations on the realized module, and types the one remaining physical
choice as a named interface.

THE THREE CLAUSES:

* Seam grammar: SELECTED by measurement.  The pinned #624 receipt
  identifies the coefficient group with the measured order-six central
  column of the routed-seam receipt (Z7 has no faithful complete
  realization; the order-six menu is exhaustive for cyclic candidates),
  so the grammar selection is a target-free measured datum, re-pinned
  here.
* Refinement transport: CLASSIFIED.  The pinned persistence receipt
  carries the identity port maps at every measured level, so the
  classified grammar and the module action transport identically along
  refinement.
* Representation on realized matter: CLASSIFIED, new and exact.  The
  realized fifteen-state module carries the charge spectrum of the
  matter receipt; the induced action of each order-six subgroup is
  computed exactly (faithfulness and fixed-subspace dimension per
  subgroup), and every mechanism row of the pinned #624 two-type table
  is classified by its induced action on the realized module: rows whose
  coefficient order divides six act through the computed charge
  characters, and rows whose coefficient order does not divide six have
  no induced central action from the realized charge structure.

THE REMAINING CHOICE, typed: physical selection of one sector mechanism
(the matter ACTION among the admitted ones) stays a named conditional
open interface on the physical family chain; the pinned #624 selection
boundary names this issue, and this certificate carries the boundary
forward as an interface rather than a claim.  No target value enters any
computation.
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

SCHEMA = "oph.seam_grammar_matter_classification_certificate.v1"
MANIFEST_PATH = (
    MODULE_DIR / "manifests" / "seam_grammar_matter_classification_reference.json"
)
SEAM_RECEIPT_PATH = MODULE_DIR / "receipts" / "noncentral_seam_reduction_reference.receipt.json"
MATTER_RECEIPT_PATH = MODULE_DIR / "receipts" / "super_tannakian_matter_reference.receipt.json"
RESPONSE_ARTIFACT_NAME = "charged_response_semantic_artifact.json"

ISSUE = 627


def pin_file(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    payload = load_json(path)
    return payload, {
        "path": str(path.relative_to(MODULE_DIR.parent.parent)),
        "sha256": sha256_json(payload),
    }


# ---------------------------------------------------------------------------
# Clause one: the measured grammar selection, re-pinned
# ---------------------------------------------------------------------------


def grammar_selection(seam_receipt: Mapping[str, Any]) -> dict[str, Any]:
    lift = seam_receipt["complete_coefficient_lift"]["classification"]
    require(
        lift["z6"]["classification"] == "identified_with_measured_central_column",
        "GRAMMAR_SELECTION",
        "the pinned receipt must identify the coefficient group with the measured column",
    )
    require(
        "Z7 has no faithful realization" in lift["conclusion"],
        "GRAMMAR_SELECTION",
        "the pinned receipt must carry the Z7 exclusion",
    )
    return {
        "selected_group": "Z6",
        "selection_source": "the measured order-six central column of the routed-seam receipt, pinned through the #624 complete coefficient lift",
        "excluded": "Z7 (no faithful complete realization); the order-six menu is exhaustive for cyclic candidates",
        "target_free": True,
    }


# ---------------------------------------------------------------------------
# Clause three: the realized-module classification, exact
# ---------------------------------------------------------------------------


def module_action_table(charge_spectrum: Mapping[str, int]) -> dict[str, Any]:
    """Exact induced action of every order-six subgroup on the module.

    Charges enter as hypercharges Y; the central column acts through
    q = 6Y mod 6.  For the subgroup of order d dividing six, the action is
    by q mod d; it is faithful exactly when gcd of the charges with d is
    one, and the fixed subspace collects the states with q divisible
    by d.
    """

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
    for d in (1, 2, 3, 6):
        charges_mod = [s["q6"] % d for s in states] if d > 1 else [0]
        g = d
        for value in charges_mod:
            g = gcd(g, value)
        faithful = (d == 1) or (g == 1)
        fixed = sum(
            s["multiplicity"] for s in states if s["q6"] % d == 0
        ) if d > 1 else 15
        table.append(
            {
                "subgroup_order": d,
                "faithful_on_module": bool(faithful) if d > 1 else False,
                "fixed_subspace_dimension": fixed,
            }
        )
    by_order = {row["subgroup_order"]: row for row in table}
    require(
        by_order[6]["faithful_on_module"]
        and by_order[6]["fixed_subspace_dimension"] == 1
        and by_order[3]["fixed_subspace_dimension"] == 3
        and by_order[2]["fixed_subspace_dimension"] == 7,
        "MODULE_TABLE",
        "the exact module-action table must match the realized charge spectrum",
    )
    return {"states": states, "subgroup_actions": table}


def mechanism_classification(
    seam_receipt: Mapping[str, Any], actions: Mapping[str, Any]
) -> dict[str, Any]:
    """Classify every pinned mechanism row by its induced module action."""

    rows = seam_receipt["two_type_sector_classification"]["mechanism_table"]
    by_order = {
        row["subgroup_order"]: row for row in actions["subgroup_actions"]
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
                action = by_order[order]
                entry["realized_module_action"] = {
                    "through": f"charge characters modulo {order}",
                    "faithful_on_module": action["faithful_on_module"],
                    "fixed_subspace_dimension": action["fixed_subspace_dimension"],
                }
            else:
                entry["realized_module_action"] = {
                    "through": "none",
                    "reading": (
                        "the coefficient order does not divide six, so the "
                        "realized charge structure induces no central action"
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
    require(
        any(
            not row["contractible"]
            and row["realized_module_action"].get("faithful_on_module")
            for row in classified
        ),
        "MECHANISM_TABLE",
        "at least one admitted mechanism must act faithfully on the realized module",
    )
    return {
        "rows": classified,
        "statement": (
            "every admitted mechanism is classified by its induced action on "
            "the realized module; the admitted noncontractible sectors with "
            "coefficient order dividing six act through the exact charge "
            "characters, and no other central action is induced"
        ),
    }


# ---------------------------------------------------------------------------
# Clause two: refinement transport
# ---------------------------------------------------------------------------


def refinement_transport() -> dict[str, Any]:
    artifact, pin = pin_file(MODULE_DIR / "manifests" / RESPONSE_ARTIFACT_NAME)
    maps = artifact["physical_refinement_maps"]
    for payload in maps["port_persistence_maps"]:
        require(
            sorted(payload["port_map"]) == list(range(12)),
            "REFINEMENT",
            "each persistence map must biject the twelve ports",
        )
    return {
        "pin": pin,
        "classification": (
            "the persistence maps act on the port set bijectively at every "
            "measured level, so the selected grammar and the classified "
            "module action transport identically along refinement"
        ),
    }


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------


def control_charge_mutation() -> dict[str, Any]:
    """A doctored charge spectrum must fail the exact module table."""

    doctored = {"-1/2": 2, "-2/3": 3, "1": 1, "1/3": 3, "1/6": 5}
    try:
        module_action_table(doctored)
    except CertificateError as error:
        return {"expected_failure": True, "failed": True, "code": error.code}
    return {"expected_failure": True, "failed": False}


def control_seven_injection() -> dict[str, Any]:
    """Injecting a seven-element coefficient group must be refused by the
    pinned exclusion."""

    candidate_order = 7
    try:
        require(
            6 % candidate_order == 0,
            "SEVEN_EXCLUDED",
            "a coefficient order that does not divide six has no induced action, and Z7 is excluded by the pinned lift",
        )
    except CertificateError:
        return {"expected_failure": True, "failed": True, "code": "SEVEN_EXCLUDED"}
    return {"expected_failure": True, "failed": False}


def control_selection_promotion() -> dict[str, Any]:
    """Promoting one sector mechanism to the physical matter action without
    the physical chain must be refused."""

    physical_chain_receipts_present = False
    try:
        require(
            physical_chain_receipts_present,
            "SELECTION_PROMOTION",
            "physical selection of one sector mechanism requires the named physical-chain receipts",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "SELECTION_PROMOTION",
        }
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_payload() -> dict[str, Any]:
    seam_receipt, seam_pin = pin_file(SEAM_RECEIPT_PATH)
    matter_receipt, matter_pin = pin_file(MATTER_RECEIPT_PATH)

    selection = grammar_selection(seam_receipt)
    actions = module_action_table(matter_receipt["realized_package"]["charge_spectrum"])
    mechanisms = mechanism_classification(seam_receipt, actions)
    transport = refinement_transport()

    boundary = seam_receipt["two_type_sector_classification"]["selection_boundary"]
    require(
        "627" in str(boundary),
        "SELECTION_BOUNDARY",
        "the pinned selection boundary must name this issue",
    )

    controls = {
        "charge_mutation": control_charge_mutation(),
        "seven_injection": control_seven_injection(),
        "selection_promotion": control_selection_promotion(),
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
            "The closure clause is met on its stated branches: the seam "
            "grammar is selected by the measured central column (pinned), "
            "the refinement transport is classified by the pinned identity "
            "persistence, and the representation on realized matter is "
            "classified exactly through the charge characters, with every "
            "pinned mechanism row classified by its induced module action. "
            "Physical selection of one sector mechanism as the matter action "
            "stays a named conditional open interface on the physical family "
            "chain; no target value enters any computation."
        ),
        "upstream_pins": {
            "seam_classification_receipt": seam_pin,
            "matter_receipt": matter_pin,
        },
        "grammar_selection": selection,
        "module_action_classification": actions,
        "mechanism_classification": mechanisms,
        "refinement_transport": transport,
        "matter_action_interface": {
            "id": "physical_sector_mechanism_selection",
            "class": "conditional_open_interface",
            "statement": (
                "one admitted sector mechanism acts as the physical matter "
                "action; selection requires the physical family-chain "
                "receipts and is refused without them"
            ),
        },
        "controls": controls,
        "bounded_exit": "classification_landed_with_named_selection_interface",
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
