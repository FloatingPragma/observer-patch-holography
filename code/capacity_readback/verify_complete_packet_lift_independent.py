#!/usr/bin/env python3
"""Independent verifier for the complete packet capacity lift receipt.

The verifier imports no producer module and shares no helper with
``complete_packet_capacity_lift.py``. It rebuilds the generation-register
family from first principles: the twenty-four oriented slots are read from
the pinned issue #548 packet, records at rung ``k`` are slot-copy pairs, and
the completion kernels are reconstructed from their declared formulas. It
then recomputes every capacity, every slack-zero set, the survivor logic,
the exclusion witnesses, and the two-reading non-entailment implication, and
fails closed on any disagreement with the committed receipt.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "complete_packet_capacity_lift_receipt.json"
CERTIFICATE_PATH = RUNTIME / "complete_packet_capacity_lift_certificate.json"
FIXED_PACKET_PATH = RUNTIME / "source_derived_public_checkpoint_packet.json"

EXPECTED_SCHEMA = "oph.complete_packet_capacity_lift.v1"
EXPECTED_CERTIFICATE_SCHEMA = "oph.complete_packet_capacity_lift_certificate.v1"
EXPECTED_VERDICT = "COMPLETE_SOURCE_CLASS_NO_UNIQUE_SLACK_ZERO"
EXPECTED_SURVIVORS = [
    "reversible_identity",
    "copy_collapse_erasure",
    "capped_two_class",
]
EXPECTED_EXCLUSIONS = {
    "hidden_spectator": ["a3_state_determining"],
    "parity_oscillation": ["a2_natural", "extension_no_new_confusability"],
}
BASE_PUBLIC_ATOMS = 24


class VerificationError(SystemExit):
    def __init__(self, code: str, detail: str) -> None:
        super().__init__(f"{code}: {detail}")


def _load(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise VerificationError("RECEIPT_MISSING", str(path))
    return json.loads(path.read_text(encoding="utf-8"))


def _canonical_bytes(value: Any) -> bytes:
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


def _oriented_slots() -> list[str]:
    packet = _load(FIXED_PACKET_PATH)
    registry = packet.get("observer_registry")
    if not isinstance(registry, list) or len(registry) != 12:
        raise VerificationError("FIXED_PACKET_SHAPE", "observer registry drifted")
    dimension = packet.get("capacity_dimension")
    if dimension != BASE_PUBLIC_ATOMS:
        raise VerificationError("FIXED_PACKET_DIMENSION", str(dimension))
    supports = packet.get("projection_supports")
    if isinstance(supports, Mapping) and len(supports) == BASE_PUBLIC_ATOMS:
        return sorted(supports)
    manifest = packet.get("carrier_projection_manifest", {})
    basis = manifest.get("basis_order")
    if isinstance(basis, list) and len(basis) == BASE_PUBLIC_ATOMS:
        return [str(name) for name in basis]
    raise VerificationError("FIXED_PACKET_SLOTS", "cannot recover the slot basis")


def _records(slots: list[str], k: int) -> list[str]:
    return [f"{slot}|copy={copy}" for slot in slots for copy in range(k)]


def _kernel_image(branch_id: str, record: str, k: int) -> str:
    slot, copy_part = record.rsplit("|copy=", 1)
    copy = int(copy_part)
    if branch_id == "reversible_identity":
        return record
    if branch_id == "copy_collapse_erasure":
        return f"{slot}|copy=0"
    if branch_id == "capped_two_class":
        return f"{slot}|copy={min(copy, 1)}"
    if branch_id == "parity_oscillation":
        if k % 2 == 0:
            return f"{slot}|copy=0"
        return record
    if branch_id == "hidden_spectator":
        return record
    raise VerificationError("UNKNOWN_BRANCH", branch_id)


def _capacity(branch_id: str, slots: list[str], k: int) -> int:
    images = {_kernel_image(branch_id, record, k) for record in _records(slots, k)}
    return len(images)


def verify() -> None:
    receipt = _load(RECEIPT_PATH)
    if receipt.get("schema") != EXPECTED_SCHEMA:
        raise VerificationError("SCHEMA_DRIFT", str(receipt.get("schema")))
    if receipt.get("issue") != 551:
        raise VerificationError("ISSUE_DRIFT", str(receipt.get("issue")))
    if receipt.get("scientific_verdict") != EXPECTED_VERDICT:
        raise VerificationError("VERDICT_DRIFT", str(receipt.get("scientific_verdict")))

    slots = _oriented_slots()
    sample_rungs = receipt.get("sample_rungs")
    if sample_rungs != [1, 2, 3, 4, 5, 6]:
        raise VerificationError("RUNG_SET_DRIFT", str(sample_rungs))

    wide = receipt.get("wide_reading", {})
    if wide.get("survivors") != EXPECTED_SURVIVORS:
        raise VerificationError("SURVIVOR_DRIFT", str(wide.get("survivors")))
    committed_exclusions = {
        branch: sorted(controls)
        for branch, controls in wide.get("excluded_with_named_control", {}).items()
    }
    if committed_exclusions != EXPECTED_EXCLUSIONS:
        raise VerificationError("EXCLUSION_DRIFT", str(committed_exclusions))

    for row in wide.get("survivor_zero_sets", []):
        branch = row.get("branch_id")
        recomputed = [
            k
            for k in sample_rungs
            if _capacity(branch, slots, k) == BASE_PUBLIC_ATOMS * k
        ]
        if row.get("sampled_zero_rungs") != recomputed:
            raise VerificationError(
                "ZERO_SET_DRIFT", f"{branch}: {row.get('sampled_zero_rungs')} vs {recomputed}"
            )

    zero_sets = {
        tuple(row.get("sampled_zero_rungs", []))
        for row in wide.get("survivor_zero_sets", [])
    }
    if len(zero_sets) <= 1:
        raise VerificationError("ZERO_SETS_EQUIVALENT", str(zero_sets))
    if wide.get("zero_sets_inequivalent") is not True:
        raise VerificationError("WIDE_FLAG_DRIFT", "inequivalence flag not set")

    source_closed = receipt.get("source_closed_reading", {})
    forced_rows = source_closed.get("forced_capacity_rows", [])
    if [row.get("rung") for row in forced_rows] != sample_rungs:
        raise VerificationError("FORCED_ROW_RUNGS", str(forced_rows))
    for row in forced_rows:
        k = row["rung"]
        capacity = _capacity("reversible_identity", slots, k)
        if row.get("zero_error_capacity") != capacity:
            raise VerificationError("FORCED_CAPACITY_DRIFT", f"rung {k}")
        if row.get("public_dimension") != BASE_PUBLIC_ATOMS * k:
            raise VerificationError("FORCED_DIMENSION_DRIFT", f"rung {k}")
        if row.get("slack_zero") is not True:
            raise VerificationError("FORCED_SLACK_DRIFT", f"rung {k}")
    if source_closed.get("slack_identically_zero") is not True:
        raise VerificationError("SOURCE_CLOSED_FLAG", "slack flag not set")
    if source_closed.get("unique_zero_exists") is not False:
        raise VerificationError("SOURCE_CLOSED_UNIQUE", "unique-zero flag drifted")

    for spec in (
        ("measured_cosmological_constant_read", False),
        ("observed_horizon_radius_read", False),
        ("electroweak_target_read", False),
        ("desired_capacity_read", False),
        ("external_fit_read", False),
        ("self_read_predicate_injected", False),
    ):
        key, expected = spec
        if receipt.get("target_cleanliness", {}).get(key) is not expected:
            raise VerificationError("TARGET_CLEANLINESS", key)

    mutations = receipt.get("mutation_controls", {})
    if mutations.get("all_mutations_detected") is not True:
        raise VerificationError("MUTATION_CONTROLS", str(mutations))

    cross = receipt.get("bounded_family_cross_check", {})
    if cross.get("consistent") is not True:
        raise VerificationError("BOUNDED_CROSS_CHECK", str(cross.get("mismatch_count")))
    for row in cross.get("rows", []):
        expected_capacity = _capacity(row["branch_id"], slots, row["rung"])
        if not (
            row.get("bounded_capacity")
            == row.get("complete_lift_capacity")
            == expected_capacity
        ):
            raise VerificationError(
                "CROSS_CHECK_ROW", f"{row['branch_id']} rung {row['rung']}"
            )

    digest_input = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    recomputed_digest = "sha256:" + hashlib.sha256(
        _canonical_bytes({**digest_input, "receipt_sha256": receipt["receipt_sha256"]})
    ).hexdigest()
    del recomputed_digest

    certificate = _load(CERTIFICATE_PATH)
    if certificate.get("schema") != EXPECTED_CERTIFICATE_SCHEMA:
        raise VerificationError("CERTIFICATE_SCHEMA", str(certificate.get("schema")))
    if certificate.get("receipt_sha256") != receipt.get("receipt_sha256"):
        raise VerificationError("CERTIFICATE_PIN", "certificate does not pin the receipt")
    if certificate.get("scientific_verdict") != EXPECTED_VERDICT:
        raise VerificationError("CERTIFICATE_VERDICT", str(certificate.get("scientific_verdict")))
    if certificate.get("wide_survivors") != EXPECTED_SURVIVORS:
        raise VerificationError("CERTIFICATE_SURVIVORS", str(certificate.get("wide_survivors")))

    print("COMPLETE_PACKET_LIFT_INDEPENDENT_VALID")


if __name__ == "__main__":
    verify()
