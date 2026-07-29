#!/usr/bin/env python3
"""Consume the issue 551 verdict in the direct issue 505 N equation."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

from capacity_indexed_source_family import (
    CERTIFICATE_PATH,
    PROJECTION_PATH,
    canonical_json_bytes,
    tagged_sha256,
)


HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
FIXED_CERTIFICATE_PATH = (
    RUNTIME / "source_derived_public_checkpoint_certificate.json"
)
INDEPENDENT_RECEIPT_PATH = (
    RUNTIME / "capacity_indexed_source_family_independent_receipt.json"
)
INDEPENDENT_CUSTODY_PATH = (
    RUNTIME / "capacity_indexed_source_family_independent_custody.json"
)
OUTPUT_PATH = RUNTIME / "direct_n_closure_verdict.json"
SCHEMA = "oph.direct_n_closure_verdict.v1"
STATUS = "NOT_EVALUABLE_INCOMPLETE_CAPACITY_SOURCE_ANTECEDENT"
FAMILY_SCHEMA = "oph.capacity_indexed_source_family_certificate.v1"
PROJECTION_SCHEMA = "oph.capacity_indexed_source_family_projection.v1"
INDEPENDENT_SCHEMA = (
    "oph.capacity_indexed_source_family_independent_receipt.v1"
)
CUSTODY_SCHEMA = (
    "oph.capacity_indexed_source_family_independent_custody.v1"
)
EXPECTED_BRANCH_IDS = {
    "capped_two_class",
    "copy_collapse_erasure",
    "hidden_spectator",
    "reversible_identity",
}
EXPECTED_CUSTODY_PATHS = {
    "data/capacity_readback/capacity_indexed_source_family_projection.json",
    "data/capacity_readback/capacity_indexed_source_family_independent_receipt.json",
    "oph_fpe/cosmology/capacity_indexed_family_verifier.py",
    "schemas/cosmology/capacity_indexed_source_family_projection.schema.json",
    "tests/test_capacity_indexed_family_verifier.py",
}
HEX40 = re.compile(r"[0-9a-f]{40}")
TAGGED_SHA256 = re.compile(r"sha256:[0-9a-f]{64}")
RAW_SHA256 = re.compile(r"[0-9a-f]{64}")


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _pin(path: Path) -> dict[str, Any]:
    data = path.read_bytes()
    return {
        "path": path.relative_to(HERE.parent.parent).as_posix(),
        "sha256": "sha256:" + hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
    }


def _raw_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def _require_digest(value: object, label: str) -> str:
    if not isinstance(value, str) or TAGGED_SHA256.fullmatch(value) is None:
        raise ValueError(f"{label} is not a tagged SHA-256 digest")
    return value


def _require_raw_digest(value: object, label: str) -> str:
    if not isinstance(value, str) or RAW_SHA256.fullmatch(value) is None:
        raise ValueError(f"{label} is not a raw SHA-256 digest")
    return value


def _validate_family(
    projection: dict[str, Any],
    family: dict[str, Any],
) -> None:
    if projection.get("schema") != PROJECTION_SCHEMA:
        raise ValueError("capacity-indexed projection schema drift")
    if projection.get("scientific_verdict") != (
        "BOUNDED_COMPLETION_CLASS_NONIDENTIFIABLE"
    ):
        raise ValueError("capacity-indexed projection verdict drift")
    if projection.get("physical_n_closure_promoted") is not False:
        raise ValueError("projection promoted physical N closure")
    target_flags = projection.get("target_cleanliness")
    if not isinstance(target_flags, dict) or any(target_flags.values()):
        raise ValueError("capacity-indexed projection is not target-clean")

    branches = projection.get("branches")
    if not isinstance(branches, list) or not branches:
        raise ValueError("capacity-indexed projection has no branches")
    branch_ids = [branch.get("branch_id") for branch in branches]
    if set(branch_ids) != EXPECTED_BRANCH_IDS or len(branch_ids) != len(
        EXPECTED_BRANCH_IDS
    ):
        raise ValueError("capacity-indexed projection branch grammar drift")

    if family.get("schema") != FAMILY_SCHEMA or family.get("issue") != 551:
        raise ValueError("capacity-indexed family identity drift")
    if family.get("status") != "BOUNDED_COMPLETION_CLASS_NONIDENTIFIABLE":
        raise ValueError("bounded capacity-completion verdict is not attained")
    required_flags = {
        "target_clean": True,
        "source_family_all_positive_rungs": True,
        "branch_grammar_complete_for_declared_counterfamily": True,
        "zero_sets_differ": True,
        "physical_n_closure_promoted": False,
        "full_a1_a3_source_class_nonidentifiability_proved": False,
    }
    for key, expected in required_flags.items():
        if family.get(key) is not expected:
            raise ValueError(f"capacity-indexed family {key} drift")
    if family.get("direct_n_status") != STATUS:
        raise ValueError("capacity-indexed direct-N status drift")
    if family.get("projection_sha256") != tagged_sha256(
        canonical_json_bytes(projection)
    ):
        raise ValueError("capacity-indexed projection pin mismatch")

    certificate_digest = family.get("certificate_sha256")
    _require_digest(certificate_digest, "capacity-indexed certificate pin")
    unhashed = dict(family)
    unhashed.pop("certificate_sha256", None)
    if certificate_digest != tagged_sha256(canonical_json_bytes(unhashed)):
        raise ValueError("capacity-indexed certificate self-pin mismatch")

    receipts = family.get("branch_receipts")
    if not isinstance(receipts, dict) or set(receipts) != EXPECTED_BRANCH_IDS:
        raise ValueError("capacity-indexed branch receipts are incomplete")
    for branch_id, branch_receipts in receipts.items():
        if not isinstance(branch_receipts, dict):
            raise ValueError(f"{branch_id} receipt packet is malformed")
        for family_name in (
            "continuation_composition",
            "sewing",
            "extension",
        ):
            rows = branch_receipts.get(family_name)
            if not isinstance(rows, list) or not rows:
                raise ValueError(
                    f"{branch_id} {family_name} receipts are empty"
                )
            if any(
                not isinstance(row, dict) or row.get("status") != "PASS"
                for row in rows
            ):
                raise ValueError(
                    f"{branch_id} {family_name} receipt did not pass"
                )

    projected_zero_sets = {
        branch["branch_id"]: branch.get("claimed_bounded_zero_set")
        for branch in branches
    }
    if family.get("bounded_zero_sets") != projected_zero_sets:
        raise ValueError("capacity-indexed bounded zero sets drift")
    if len(
        {
            canonical_json_bytes(rows)
            for rows in projected_zero_sets.values()
        }
    ) <= 1:
        raise ValueError("capacity-indexed zero sets do not differ")


def _validate_independent_replay(
    projection: dict[str, Any],
    independent: dict[str, Any],
    custody: dict[str, Any],
) -> None:
    if (
        independent.get("schema") != INDEPENDENT_SCHEMA
        or independent.get("issue") != 551
        or independent.get("status") != "PASS"
    ):
        raise ValueError("independent finite replay did not pass")
    if independent.get("target_clean") is not True:
        raise ValueError("independent finite replay is not target-clean")
    if independent.get("scientific_verdict_replayed") != (
        "BOUNDED_COMPLETION_CLASS_NONIDENTIFIABLE"
    ):
        raise ValueError("independent replay verdict drift")
    if set(independent.get("branch_ids_replayed", [])) != EXPECTED_BRANCH_IDS:
        raise ValueError("independent replay branch coverage is incomplete")
    scope = independent.get("scope")
    if scope != {
        "all_positive_integer_rungs_proved": False,
        "finite_sample_replay": True,
        "full_a1_a3_packet_lift_replayed": False,
        "physical_n_closure_promoted": False,
        "producer_implementation_independent": True,
    }:
        raise ValueError("independent replay scope drift")
    if independent.get("projection_sha256") != _raw_sha256(PROJECTION_PATH):
        raise ValueError("independent replay used different projection bytes")
    if independent.get("shared_source_signature_sha256") != projection.get(
        "shared_source_signature_sha256"
    ):
        raise ValueError("independent replay source signature drift")
    if independent.get("upstream_pins") != projection.get("upstream_pins"):
        raise ValueError("independent replay upstream pins drift")
    _require_digest(
        independent.get("full_replay_report_sha256"),
        "independent replay report pin",
    )

    if (
        custody.get("schema") != CUSTODY_SCHEMA
        or custody.get("issue") != 551
        or custody.get("status") != "PASS"
    ):
        raise ValueError("independent replay custody is incomplete")
    commit = custody.get("commit")
    if (
        not isinstance(commit, str)
        or HEX40.fullmatch(commit) is None
        or commit == "0" * 40
    ):
        raise ValueError("independent replay custody commit is invalid")
    if custody.get("repository") != (
        "https://github.com/muellerberndt/oph-physics-sim"
    ):
        raise ValueError("independent replay custody repository drift")

    artifacts = custody.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise ValueError("independent replay custody artifacts are empty")
    artifact_map: dict[str, str] = {}
    for row in artifacts:
        if (
            not isinstance(row, dict)
            or set(row) != {"path", "sha256"}
            or not isinstance(row.get("path"), str)
            or row["path"] in artifact_map
        ):
            raise ValueError("independent replay custody artifact is malformed")
        artifact_map[row["path"]] = _require_raw_digest(
            row.get("sha256"),
            f"custody artifact {row['path']}",
        )
    if set(artifact_map) != EXPECTED_CUSTODY_PATHS:
        raise ValueError("independent replay custody artifact set drift")
    if artifact_map[
        "data/capacity_readback/capacity_indexed_source_family_projection.json"
    ] != _raw_sha256(PROJECTION_PATH).removeprefix("sha256:"):
        raise ValueError("custody projection hash mismatch")
    if artifact_map[
        "data/capacity_readback/capacity_indexed_source_family_independent_receipt.json"
    ] != _raw_sha256(INDEPENDENT_RECEIPT_PATH).removeprefix("sha256:"):
        raise ValueError("custody receipt hash mismatch")
    if artifact_map[
        "oph_fpe/cosmology/capacity_indexed_family_verifier.py"
    ] != str(independent.get("independent_verifier_sha256")).removeprefix(
        "sha256:"
    ):
        raise ValueError("custody verifier hash mismatch")
    if artifact_map[
        "schemas/cosmology/capacity_indexed_source_family_projection.schema.json"
    ] != str(independent.get("schema_sha256")).removeprefix("sha256:"):
        raise ValueError("custody schema hash mismatch")

    verification = custody.get("verification")
    if not isinstance(verification, dict):
        raise ValueError("independent replay custody verification is missing")
    if (
        verification.get("finite_projection_replayed_independently") is not True
        or verification.get("all_positive_integer_rungs_proved_by_replay")
        is not False
        or verification.get("full_a1_a3_packet_lift_replayed") is not False
    ):
        raise ValueError("independent replay custody scope drift")


def _classify_fiber_branch(
    name: str,
    row: dict[str, Any],
    dimension: int,
) -> dict[str, Any]:
    capacity_set = row.get("terminal_fiber_capacity_set")
    existentially_closed = (
        isinstance(capacity_set, list) and dimension in capacity_set
    )
    robustly_closed = capacity_set == [dimension]
    return {
        "branch": name,
        "fiber_kind": row["fiber_kind"],
        "parent_status": row["status"],
        "capacity_set": capacity_set,
        "existentially_closed": existentially_closed,
        "robustly_closed": robustly_closed,
    }


def build_verdict() -> dict[str, Any]:
    fixed = _load(FIXED_CERTIFICATE_PATH)
    projection = _load(PROJECTION_PATH)
    family = _load(CERTIFICATE_PATH)
    independent = _load(INDEPENDENT_RECEIPT_PATH)
    custody = _load(INDEPENDENT_CUSTODY_PATH)

    if fixed.get("status") != "PASS" or fixed.get("issue") != 548:
        raise ValueError("fixed-cutoff parent is not attained")
    _validate_family(projection, family)
    _validate_independent_replay(projection, independent, custody)

    fiber_controls = fixed["controls"]["terminal_fibers"]
    required_fibers = {"empty", "incomplete", "ambiguous", "singleton"}
    if set(fiber_controls) - {"status"} != required_fibers:
        raise ValueError("fixed-cutoff fiber classifier is incomplete")
    expected_fibers = {
        "empty": ("EMPTY", "NO_CAPACITY_READBACK", None),
        "incomplete": ("INCOMPLETE", "INCOMPLETE_TERMINAL_FIBER", None),
        "ambiguous": ("AMBIGUOUS", "AMBIGUOUS_CAPACITY_READBACK", [1, 24]),
        "singleton": ("SINGLETON", "PASS", [24]),
    }
    for name, (kind, status, capacity_set) in expected_fibers.items():
        row = fiber_controls[name]
        if (
            row.get("fiber_kind") != kind
            or row.get("status") != status
            or row.get("terminal_fiber_capacity_set") != capacity_set
        ):
            raise ValueError(f"fixed-cutoff {name} fiber classifier drift")
    if fixed["capacity"]["exact_zero_error_capacity"] != 24:
        raise ValueError("fixed-cutoff capacity drift")

    pins = [
        _pin(FIXED_CERTIFICATE_PATH),
        _pin(PROJECTION_PATH),
        _pin(CERTIFICATE_PATH),
        _pin(INDEPENDENT_RECEIPT_PATH),
        _pin(INDEPENDENT_CUSTODY_PATH),
    ]
    verdict = {
        "schema": SCHEMA,
        "issues": [505, 551],
        "status": STATUS,
        "official_equation": "N_out=log(M0(U_N_in)); closure requires N_out=N_in",
        "typed_coordinates": {
            "trial_carrier_dimension": "D",
            "input_coordinate": "N_in=log(D)",
            "correctable_record_count": "M0",
            "output_coordinate": "N_out=log(M0)",
        },
        "fixed_cutoff_result": {
            "D": 24,
            "M0": 24,
            "status": "SOURCE_DERIVED_FIXED_CUTOFF_PHYSICAL_PACKET",
            "cosmic_value_selected": False,
        },
        "capacity_indexed_result": {
            "source_class": projection["shared_source"][
                "branch_completion_scope"
            ],
            "branch_ids": [
                branch["branch_id"] for branch in projection["branches"]
            ],
            "zero_sets_differ": family["zero_sets_differ"],
            "unique_source_zero_entailed": False,
            "strange_loop_identity_rejected": False,
        },
        "fiber_classifier": {
            name: fiber_controls[name]["fiber_kind"]
            for name in sorted(required_fibers)
        },
        "fiber_closure_branches": {
            name: _classify_fiber_branch(name, fiber_controls[name], 24)
            for name in sorted(required_fibers)
        },
        "source_controls": {
            "target_clean": family["target_clean"],
            "constructor_reads_desired_capacity": False,
            "identity_control": "reversible_identity"
            in family["branch_receipts"],
            "erasure_control": "copy_collapse_erasure"
            in family["branch_receipts"],
            "multiple_zero_control": "capped_two_class"
            in family["branch_receipts"],
            "hidden_spectator_control": "hidden_spectator"
            in family["branch_receipts"],
            "cyclic_permutation_control": (
                fixed["controls"]["cyclic_permutation"]["status"] == "PASS"
            ),
            "marginal_coupling_control": (
                fixed["controls"]["alternative_joint_coupling"]["status"]
                == "PASS"
            ),
            "full_support_noise_control": (
                fixed["controls"]["full_support_noise"]["status"] == "PASS"
            ),
            "ambiguous_fiber_control": (
                fiber_controls["ambiguous"]["fiber_kind"] == "AMBIGUOUS"
            ),
            "capacity_extension_control": all(
                receipt["status"] == "PASS"
                for branch in family["branch_receipts"].values()
                for receipt in branch["extension"]
            ),
            "sewing_control": all(
                receipt["status"] == "PASS"
                for branch in family["branch_receipts"].values()
                for receipt in branch["sewing"]
            ),
            "continuation_composition_control": all(
                receipt["status"] == "PASS"
                for branch in family["branch_receipts"].values()
                for receipt in branch["continuation_composition"]
            ),
            "independent_finite_replay": independent["scope"][
                "finite_sample_replay"
            ],
            "all_rung_lean_theorem": (
                "OPH.CapacityNonidentifiability."
                "boundedCompletionClass_doesNotForceUniqueZero"
            ),
        },
        "comparison_boundary": {
            "direct_numeric_N_emitted": False,
            "cosmological_comparison_permitted": False,
            "horizon_record_attachment_evaluable": False,
            "electroweak_bridge_may_repair_result": False,
            "next_forecast_action": (
                "keep direct N outside scoring while the complete A1-A3 "
                "packet lift remains absent"
            ),
        },
        "remaining_positive_route": family["remaining_positive_route"],
        "parent_pins": pins,
        "claim_boundary": (
            "The declared base-agreement, positivity, and carrier-bound "
            "completion class does not select a unique direct N. The complete "
            "A1-A3 capacity source antecedent has not been lifted across the "
            "regulator family, so issues 551 and 505 remain open."
        ),
    }
    verdict["verdict_sha256"] = tagged_sha256(canonical_json_bytes(verdict))
    return verdict


def write_runtime() -> None:
    OUTPUT_PATH.write_bytes(canonical_json_bytes(build_verdict()))


def verify_runtime() -> None:
    expected = canonical_json_bytes(build_verdict())
    if OUTPUT_PATH.read_bytes() != expected:
        raise SystemExit("direct N closure verdict is stale")
    if _load(OUTPUT_PATH)["status"] != STATUS:
        raise SystemExit("direct N closure verdict did not attain its bounded exit")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        write_runtime()
    if args.verify:
        verify_runtime()
    if not args.write and not args.verify:
        print(canonical_json_bytes(build_verdict()).decode("ascii"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
