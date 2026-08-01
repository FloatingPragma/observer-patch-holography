#!/usr/bin/env python3
"""Independent verifier for the bounded generation-register packet.

This verifier imports no producer module. It recomputes the finite capacity
rows and control classification, verifies the receipt and certificate
self-hashes, checks every upstream byte pin, and enforces the scientific scope.
In particular, it fails if sampled controls are promoted to a universal
all-rung A1--A3 membership theorem or an executable-to-Lean bridge.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "complete_packet_capacity_lift_receipt.json"
CERTIFICATE_PATH = RUNTIME / "complete_packet_capacity_lift_certificate.json"
FIXED_PACKET_PATH = RUNTIME / "source_derived_public_checkpoint_packet.json"
FAMILY_PRODUCER_PATH = HERE / "capacity_indexed_source_family.py"
SPEC_PATH = HERE / "F_READBACK_SPEC.md"

EXPECTED_SCHEMA = "oph.complete_packet_capacity_lift.v2"
EXPECTED_CERTIFICATE_SCHEMA = "oph.complete_packet_capacity_lift_certificate.v2"
EXPECTED_VERDICT = "BOUNDED_GENERATION_REGISTER_COUNTERMODEL__UNIVERSAL_MEMBERSHIP_OPEN"
EXPECTED_DIRECT_N_STATUS = "NOT_EVALUABLE_INCOMPLETE_CAPACITY_SOURCE_ANTECEDENT"
EXPECTED_SOURCE_RULE = "oph.public-record-capacity.generation-register-lift.v1"
EXPECTED_SAMPLE_RUNGS = [1, 2, 3, 4, 5, 6]
EXPECTED_BRANCHES = [
    "reversible_identity",
    "copy_collapse_erasure",
    "capped_two_class",
    "hidden_spectator",
    "parity_oscillation",
]
EXPECTED_PASSERS = [
    "reversible_identity",
    "copy_collapse_erasure",
    "capped_two_class",
]
EXPECTED_EXCLUSIONS = {
    "hidden_spectator": ["a3_state_determining"],
    "parity_oscillation": ["a2_natural", "extension_no_new_confusability"],
}
EXPECTED_CONTROL_KEYS = {
    "terminal_fiber",
    "sections",
    "histories",
    "publicness_frozen",
    "manifest_closed",
    "source_ancestry_complete",
    "kernels_marginal_consistent",
    "a2_natural",
    "a3_state_determining",
    "extension_no_new_confusability",
    "refinement_stable",
    "sewing_exact",
}
EXPECTED_MUTATION_KEYS = {
    "self_read_injection_detected",
    "kernel_tamper_changes_capacity",
    "ancestry_mislabel_detected",
    "oscillation_extension_square_fails",
    "spectator_determinacy_fails",
    "history_completeness_checked",
    "sewing_tamper_detected",
    "refinement_tamper_detected",
    "section_connectivity_required",
    "all_mutations_detected",
}
EXPECTED_FORMULAS = {
    "reversible_identity": "every positive rung",
    "copy_collapse_erasure": "rung 1 only",
    "capped_two_class": "rungs 1 and 2 only",
}
BASE_PUBLIC_ATOMS = 24


class VerificationError(SystemExit):
    def __init__(self, code: str, detail: str) -> None:
        super().__init__(f"{code}: {detail}")


def _canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def _tagged_sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise VerificationError("RECEIPT_MISSING", str(path))
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise VerificationError("JSON_INVALID", f"{path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise VerificationError("JSON_SHAPE", str(path))
    return payload


def _verify_self_hash(payload: Mapping[str, Any], key: str, label: str) -> None:
    committed = payload.get(key)
    unhashed = dict(payload)
    unhashed.pop(key, None)
    expected = _tagged_sha256(_canonical_json_bytes(unhashed))
    if committed != expected:
        raise VerificationError(f"{label}_SELF_HASH", f"{committed} vs {expected}")


def _oriented_slots() -> list[str]:
    packet = _load(FIXED_PACKET_PATH)
    registry = packet.get("observer_registry")
    if not isinstance(registry, list) or len(registry) != 12:
        raise VerificationError("FIXED_PACKET_SHAPE", "observer registry drifted")
    if packet.get("capacity_dimension") != BASE_PUBLIC_ATOMS:
        raise VerificationError(
            "FIXED_PACKET_DIMENSION", str(packet.get("capacity_dimension"))
        )
    supports = packet.get("projection_supports")
    if isinstance(supports, Mapping) and len(supports) == BASE_PUBLIC_ATOMS:
        return sorted(str(slot) for slot in supports)
    basis = packet.get("carrier_projection_manifest", {}).get("basis_order")
    if isinstance(basis, list) and len(basis) == BASE_PUBLIC_ATOMS:
        return [str(slot) for slot in basis]
    raise VerificationError("FIXED_PACKET_SLOTS", "cannot recover slot basis")


def _records(slots: list[str], k: int) -> list[str]:
    return [f"{slot}|copy={copy}" for slot in slots for copy in range(k)]


def _kernel_image(branch_id: str, record: str, k: int) -> str:
    slot, copy_part = record.rsplit("|copy=", 1)
    copy = int(copy_part)
    if branch_id in {"reversible_identity", "hidden_spectator"}:
        return record
    if branch_id == "copy_collapse_erasure":
        return f"{slot}|copy=0"
    if branch_id == "capped_two_class":
        return f"{slot}|copy={min(copy, 1)}"
    if branch_id == "parity_oscillation":
        return record if k % 2 == 1 else f"{slot}|copy=0"
    raise VerificationError("UNKNOWN_BRANCH", branch_id)


def _capacity(branch_id: str, slots: list[str], k: int) -> int:
    return len({_kernel_image(branch_id, record, k) for record in _records(slots, k)})


def _expected_control(branch: str, k: int, control: str) -> bool:
    if control == "source_ancestry_complete":
        if branch in {"reversible_identity", "hidden_spectator"}:
            return True
        if branch == "copy_collapse_erasure":
            return k == 1
        if branch == "capped_two_class":
            return k <= 2
        if branch == "parity_oscillation":
            return k % 2 == 1
    if control == "a2_natural":
        return branch != "parity_oscillation" or k == 1
    if control == "a3_state_determining":
        return branch != "hidden_spectator"
    if control == "extension_no_new_confusability":
        return not (branch == "parity_oscillation" and k >= 3 and k % 2 == 1)
    return True


def _verify_control_table(receipt: Mapping[str, Any]) -> None:
    assessment = receipt.get("sampled_control_assessment", {})
    table = assessment.get("control_table")
    if (
        not isinstance(table, list)
        or [row.get("branch_id") for row in table] != EXPECTED_BRANCHES
    ):
        raise VerificationError("CONTROL_BRANCHES", str(table))

    derived_passers: list[str] = []
    derived_exclusions: dict[str, list[str]] = {}
    for branch_row in table:
        branch = str(branch_row["branch_id"])
        rows = branch_row.get("per_rung_controls")
        if not isinstance(rows, Mapping) or list(rows) != [
            str(k) for k in EXPECTED_SAMPLE_RUNGS
        ]:
            raise VerificationError("CONTROL_RUNGS", branch)
        transported_pass = True
        ancestry_pass = True
        failed_controls: set[str] = set()
        for k in EXPECTED_SAMPLE_RUNGS:
            row = rows[str(k)]
            if not isinstance(row, Mapping) or set(row) != EXPECTED_CONTROL_KEYS:
                raise VerificationError("CONTROL_KEYS", f"{branch} rung {k}")
            for control in EXPECTED_CONTROL_KEYS:
                expected = _expected_control(branch, k, control)
                if row.get(control) is not expected:
                    raise VerificationError(
                        "CONTROL_VALUE", f"{branch} rung {k} {control}"
                    )
                if control == "source_ancestry_complete":
                    ancestry_pass = ancestry_pass and expected
                else:
                    transported_pass = transported_pass and expected
                    if not expected:
                        failed_controls.add(control)
        if (
            branch_row.get("passes_all_sampled_transported_controls")
            is not transported_pass
        ):
            raise VerificationError("CONTROL_PASS_FLAG", branch)
        if branch_row.get("passes_source_ancestry_on_all_sampled_rungs") is not (
            transported_pass and ancestry_pass
        ):
            raise VerificationError("SOURCE_CLOSED_FLAG", branch)
        if transported_pass:
            derived_passers.append(branch)
        else:
            derived_exclusions[branch] = sorted(failed_controls)

    if derived_passers != EXPECTED_PASSERS:
        raise VerificationError("CONTROL_PASSERS", str(derived_passers))
    if (
        assessment.get("branches_passing_all_sampled_transported_controls")
        != derived_passers
    ):
        raise VerificationError("PASSER_DECLARATION", str(assessment))
    if (
        assessment.get("branches_failing_sampled_transported_controls")
        != derived_exclusions
    ):
        raise VerificationError("EXCLUSION_DECLARATION", str(derived_exclusions))
    if derived_exclusions != EXPECTED_EXCLUSIONS:
        raise VerificationError("EXCLUSION_DRIFT", str(derived_exclusions))
    if assessment.get("does_not_certify_universal_membership") is not True:
        raise VerificationError("SAMPLED_SCOPE_PROMOTED", str(assessment))


def _verify_capacity_arithmetic(receipt: Mapping[str, Any], slots: list[str]) -> None:
    arithmetic = receipt.get("all_rung_capacity_arithmetic", {})
    rows = arithmetic.get("branch_zero_sets")
    if (
        not isinstance(rows, list)
        or [row.get("branch_id") for row in rows] != EXPECTED_PASSERS
    ):
        raise VerificationError("ARITHMETIC_BRANCHES", str(rows))
    for row in rows:
        branch = str(row["branch_id"])
        recomputed = [
            k
            for k in EXPECTED_SAMPLE_RUNGS
            if _capacity(branch, slots, k) == BASE_PUBLIC_ATOMS * k
        ]
        if row.get("sampled_zero_rungs") != recomputed:
            raise VerificationError("ZERO_SET_DRIFT", f"{branch}: {recomputed}")
        if row.get("all_rung_formula") != EXPECTED_FORMULAS[branch]:
            raise VerificationError("FORMULA_DRIFT", branch)

    identity_rows = arithmetic.get("sampled_identity_rows")
    if (
        not isinstance(identity_rows, list)
        or [row.get("rung") for row in identity_rows] != EXPECTED_SAMPLE_RUNGS
    ):
        raise VerificationError("IDENTITY_ROWS", str(identity_rows))
    for row in identity_rows:
        k = int(row["rung"])
        expected = _capacity("reversible_identity", slots, k)
        if row != {
            "rung": k,
            "branch_id": "reversible_identity",
            "public_dimension": BASE_PUBLIC_ATOMS * k,
            "zero_error_capacity": expected,
            "slack_zero": True,
        }:
            raise VerificationError("IDENTITY_ROW_DRIFT", str(row))
    required = {
        "identity_capacity_formula": "M0(k)=24*k for every positive k",
        "identity_slack_zero_set": "every positive rung",
        "identity_has_unique_positive_slack_zero": False,
        "inequivalent_formula_zero_sets": True,
        "all_positive_rung_arithmetic_proved_in_lean": True,
    }
    for key, expected in required.items():
        if arithmetic.get(key) != expected:
            raise VerificationError("ARITHMETIC_SCOPE", key)


def _verify_open_membership(receipt: Mapping[str, Any]) -> None:
    sample_scope = receipt.get("sample_scope", {})
    if sample_scope != {
        "kind": "finite_audit_only",
        "rungs": EXPECTED_SAMPLE_RUNGS,
        "universal_all_rung_membership_inferred": False,
    }:
        raise VerificationError("SAMPLE_SCOPE", str(sample_scope))
    status = receipt.get("source_contract_status", {})
    required_false = (
        "universal_all_rung_membership_proved",
        "executable_lean_membership_bridge_proved",
        "complete_a1_a3_source_class_nonidentifiability_proved",
        "positive_unique_zero_proved",
    )
    if any(status.get(key) is not False for key in required_false):
        raise VerificationError("UNIVERSAL_MEMBERSHIP_PROMOTED", str(status))
    if status.get("direct_n_status") != EXPECTED_DIRECT_N_STATUS:
        raise VerificationError("DIRECT_N_STATUS", str(status.get("direct_n_status")))
    missing = status.get("missing_obligations")
    if not isinstance(missing, list) or len(missing) != 3 or not all(missing):
        raise VerificationError("MISSING_OBLIGATIONS", str(missing))
    conclusion = receipt.get("bounded_conclusion", {})
    if conclusion != {
        "status": EXPECTED_VERDICT,
        "unique_direct_n_selected": False,
        "cosmic_value_emitted": False,
        "complete_source_exit_attained": False,
    }:
        raise VerificationError("BOUNDED_CONCLUSION", str(conclusion))


def _verify_supporting_evidence(receipt: Mapping[str, Any], slots: list[str]) -> None:
    target_flags = receipt.get("target_cleanliness")
    if (
        not isinstance(target_flags, Mapping)
        or set(target_flags)
        != {
            "measured_cosmological_constant_read",
            "observed_horizon_radius_read",
            "electroweak_target_read",
            "desired_capacity_read",
            "external_fit_read",
            "self_read_predicate_injected",
        }
        or any(target_flags.values())
    ):
        raise VerificationError("TARGET_CLEANLINESS", str(target_flags))

    mutations = receipt.get("mutation_controls")
    if not isinstance(mutations, Mapping) or set(mutations) != EXPECTED_MUTATION_KEYS:
        raise VerificationError("MUTATION_KEYS", str(mutations))
    if any(value is not True for value in mutations.values()):
        raise VerificationError("MUTATION_CONTROLS", str(mutations))

    cross = receipt.get("bounded_family_cross_check", {})
    if cross.get("consistent") is not True or cross.get("mismatch_count") != 0:
        raise VerificationError("BOUNDED_CROSS_CHECK", str(cross))
    expected_pairs = [(branch, k) for branch in EXPECTED_PASSERS for k in (1, 2, 3, 4)]
    rows = cross.get("rows")
    if (
        not isinstance(rows, list)
        or [(row.get("branch_id"), row.get("rung")) for row in rows] != expected_pairs
    ):
        raise VerificationError("CROSS_CHECK_ROWS", str(rows))
    for row in rows:
        expected = _capacity(str(row["branch_id"]), slots, int(row["rung"]))
        if not (
            row.get("bounded_capacity") == row.get("complete_lift_capacity") == expected
            and row.get("agree") is True
        ):
            raise VerificationError("CROSS_CHECK_ROW", str(row))

    witnesses = receipt.get("exclusion_witnesses", {})
    spectator = witnesses.get("hidden_spectator", {}).get("witness", {})
    parity = witnesses.get("parity_oscillation", {})
    if spectator.get("raw_family_multiplicity") != 4096:
        raise VerificationError("SPECTATOR_WITNESS", str(spectator))
    if parity.get("a2_extension_square", {}).get("failure_counts_by_rung") != {
        "2": 24,
        "3": 48,
    }:
        raise VerificationError("OSCILLATION_A2_WITNESS", str(parity))
    if parity.get("no_new_confusability", {}).get("new_edge_count") != 72:
        raise VerificationError("OSCILLATION_EXTENSION_WITNESS", str(parity))


def _verify_upstream_pins(receipt: Mapping[str, Any]) -> None:
    pins = receipt.get("upstream_pins")
    expected = {
        "fixed_packet_sha256": _tagged_sha256(FIXED_PACKET_PATH.read_bytes()),
        "family_producer_sha256": _tagged_sha256(FAMILY_PRODUCER_PATH.read_bytes()),
        "readback_spec_sha256": _tagged_sha256(SPEC_PATH.read_bytes()),
    }
    if pins != expected:
        raise VerificationError("UPSTREAM_PIN", f"{pins} vs {expected}")


def verify() -> None:
    receipt = _load(RECEIPT_PATH)
    if receipt.get("schema") != EXPECTED_SCHEMA or receipt.get("issue") != 551:
        raise VerificationError("RECEIPT_IDENTITY", str(receipt.get("schema")))
    if receipt.get("source_rule_id") != EXPECTED_SOURCE_RULE:
        raise VerificationError("SOURCE_RULE", str(receipt.get("source_rule_id")))
    if receipt.get("scientific_verdict") != EXPECTED_VERDICT:
        raise VerificationError("VERDICT_DRIFT", str(receipt.get("scientific_verdict")))
    _verify_self_hash(receipt, "receipt_sha256", "RECEIPT")

    slots = _oriented_slots()
    _verify_control_table(receipt)
    _verify_capacity_arithmetic(receipt, slots)
    _verify_open_membership(receipt)
    _verify_supporting_evidence(receipt, slots)
    _verify_upstream_pins(receipt)

    certificate = _load(CERTIFICATE_PATH)
    if certificate.get("schema") != EXPECTED_CERTIFICATE_SCHEMA:
        raise VerificationError("CERTIFICATE_SCHEMA", str(certificate.get("schema")))
    if certificate.get("issue") != 551:
        raise VerificationError("CERTIFICATE_ISSUE", str(certificate.get("issue")))
    if certificate.get("scientific_verdict") != EXPECTED_VERDICT:
        raise VerificationError(
            "CERTIFICATE_VERDICT", str(certificate.get("scientific_verdict"))
        )
    if certificate.get("receipt_sha256") != receipt.get("receipt_sha256"):
        raise VerificationError("CERTIFICATE_PIN", "receipt pin mismatch")
    if certificate.get("sampled_control_passers") != EXPECTED_PASSERS:
        raise VerificationError("CERTIFICATE_PASSERS", str(certificate))
    status = receipt["source_contract_status"]
    if certificate.get("open_bridges") != status["missing_obligations"]:
        raise VerificationError("CERTIFICATE_OPEN_BRIDGES", str(certificate))
    if certificate.get("direct_n_status") != EXPECTED_DIRECT_N_STATUS:
        raise VerificationError("CERTIFICATE_DIRECT_N", str(certificate))
    _verify_self_hash(certificate, "certificate_sha256", "CERTIFICATE")

    print("BOUNDED_PACKET_LIFT_INDEPENDENT_VALID")


if __name__ == "__main__":
    verify()
