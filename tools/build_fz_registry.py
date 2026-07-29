"""Build and validate the frozen-prediction ladder surface (issue #607).

The machine-readable register is ``claims/frozen_prediction_register.json``.
This tool validates it fail-closed and renders
``docs/FROZEN_PREDICTION_LADDER.md``; ``--check`` fails when the committed
page differs from the render, and the mandatory suite runs that check.

Fail-closed rules: frozen rows carry a parseable, non-future UTC freeze time,
custody, a typed attestation state, content hash, kill band, and comparison
protocol; pending rows carry an owning issue that is open in the committed
snapshot and a milestone. Retrospective results occupy a separate collection;
their former reservations cannot also occur as ladder rows. The issue-506
record is checked against a fresh replay of its canonical producer as well as
its recomputed payload digest. Committed custody contracts bind the source-side
FZ-02 receipt and Lean module even in an isolated clone. When the sibling
oph-meta custody checkout is present, the tool additionally verifies every
manifest artifact, detached OpenTimestamps digest, attestation class, and the
append-only FZ-02 erratum. When it is absent, the tool reports
``external_custody_not_present`` explicitly rather than claiming that the
external artifact set was verified.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
REGISTER_PATH = ROOT / "claims" / "frozen_prediction_register.json"
SURFACE_PATH = ROOT / "docs" / "FROZEN_PREDICTION_LADDER.md"
SNAPSHOT_PATH = ROOT / "tracking" / "open_issues" / "open_problem_ledger.json"
FZ02_RECEIPT_PATH = (
    ROOT
    / "code"
    / "a5_closure"
    / "receipts"
    / "a5_angular_multiplet_reference.receipt.json"
)
FZ04_VERDICT_REL = (
    "code/particles/alpha_hvp_audit/outputs/alpha_hvp_class_verdict.json"
)
FZ04_VERDICT_PATH = ROOT / FZ04_VERDICT_REL
FZ04_BUILDER_PATH = (
    ROOT
    / "code"
    / "particles"
    / "alpha_hvp_audit"
    / "build_alpha_hvp_verdict.py"
)
DEFAULT_CUSTODY_ROOT = ROOT.parent

SCHEMA = "oph.frozen_prediction_register.v3"
STATUSES = {
    "frozen_attested",
    "frozen_stamped_upgrade_pending",
    "standing_frozen",
    "registered_pending_freeze",
    "resource_deferred",
}
FROZEN_STATUSES = {
    "frozen_attested",
    "frozen_stamped_upgrade_pending",
    "standing_frozen",
}

ROW_KEYS = {
    "id",
    "content",
    "status",
    "frozen_utc",
    "custody",
    "attestation",
    "content_sha256",
    "kill_band",
    "comparison_protocol",
    "owning_issue",
    "milestone",
}
REGISTER_KEYS = {
    "schema",
    "issue",
    "generated_surface",
    "policy",
    "external_custody_contracts",
    "retrospective_results",
    "rows",
}
RETROSPECTIVE_RESULT_KEYS = {
    "id",
    "former_ladder_reservation",
    "content",
    "status",
    "payload_path",
    "payload_sha256",
    "comparison_protocol",
    "evidential_boundary",
    "owning_issue",
    "milestone",
}
COMMON_CONTRACT_KEYS = {
    "rows",
    "custody_path",
    "registration_manifest",
    "registration_manifest_sha256",
    "attestation_state",
    "artifact_sha256",
    "in_repo_artifact_sha256",
}
FZ02_CONTRACT_EXTRA_KEYS = {
    "custody_commit",
    "custody_commit_utc",
    "source_commit",
    "custody_erratum",
    "custody_erratum_sha256",
    "target_file",
    "target_block_sha256",
    "target_payload_sha256",
}
ATTESTATION_STATES = {"calendar_pending", "bitcoin_attested"}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")

FZ04_SCOPE = {
    "comparison_timing": "retrospective",
    "prospective_freeze": False,
    "independent_hvp_implementation_supplied": False,
    "empirical_input_promoted_to_source_output": False,
    "physical_alpha_prediction_emitted": False,
}
FZ04_CLAIM_BOUNDARY = (
    "The multi-class independent alpha/HVP test is not evaluable. One "
    "byte-pinned KNT19 accounting row is compatible under a secondary "
    "arithmetic replay. Raw-dispersive, independent-code, and lattice-HVP "
    "classes lack frozen repository ingests. The result is retrospective and "
    "supplies neither a prospective freeze nor a physical OPH alpha prediction."
)

# DetachedTimestampFile header followed by the SHA-256 operation tag and the
# 32-byte digest of the paired file. Reading this prefix does not require the
# optional opentimestamps Python package or network access.
OTS_DETACHED_HEADER = bytes.fromhex(
    "004f70656e54696d657374616d707300" "0050726f6f6600bf89e2e884e8929401"
)
OTS_SHA256_TAG = b"\x08"
OTS_PENDING_ATTESTATION_TAG = bytes.fromhex("83dfe30d2ef90c8e")
OTS_BITCOIN_ATTESTATION_TAG = bytes.fromhex("0588960d73d71901")


def fail(message: str) -> None:
    raise SystemExit(f"frozen-prediction register: {message}")


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        fail(f"missing input {path.relative_to(ROOT)}")
    except json.JSONDecodeError as error:
        fail(f"invalid JSON in {path.relative_to(ROOT)}: {error}")
    raise AssertionError("unreachable")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode("utf-8")


def rebuild_issue506_verdict() -> dict[str, Any]:
    """Replay the canonical producer instead of trusting its stored digest."""

    spec = importlib.util.spec_from_file_location(
        "_oph_alpha_hvp_verdict_builder", FZ04_BUILDER_PATH
    )
    if spec is None or spec.loader is None:
        fail("cannot import the issue-506 canonical verdict producer")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
        rebuilt = module.build_verdict()
    except Exception as error:
        fail(f"issue-506 canonical verdict replay failed: {error}")
    if not isinstance(rebuilt, dict):
        fail("issue-506 canonical verdict producer did not return an object")
    return rebuilt


def require_sha256(value: Any, where: str) -> str:
    if not isinstance(value, str) or not SHA256_RE.fullmatch(value):
        fail(f"{where} must be a lowercase SHA-256 hex digest")
    return value


def parse_utc(value: Any, where: str, *, reject_future: bool = True) -> datetime:
    if not isinstance(value, str) or not value.endswith("Z"):
        fail(f"{where} must be an ISO-8601 UTC timestamp ending in Z")
    try:
        parsed = datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError:
        fail(f"{where} is not a valid ISO-8601 UTC timestamp")
    if parsed.utcoffset() != timezone.utc.utcoffset(parsed):
        fail(f"{where} must be UTC")
    if reject_future and parsed > datetime.now(timezone.utc):
        fail(f"{where} cannot be in the future")
    return parsed


def validate_hash_mapping(mapping: Any, where: str) -> dict[str, str]:
    if not isinstance(mapping, dict):
        fail(f"{where} must be an object")
    for path, digest in mapping.items():
        if not isinstance(path, str) or not path or Path(path).is_absolute():
            fail(f"{where} keys must be nonempty relative paths")
        require_sha256(digest, f"{where}[{path!r}]")
    return mapping


def validate_custody_contracts(
    register: dict, rows_by_id: dict[str, dict]
) -> dict[str, dict]:
    contracts = register.get("external_custody_contracts")
    if not isinstance(contracts, dict) or set(contracts) != {"FZ-01", "FZ-02", "FZ-10"}:
        fail("external_custody_contracts must contain exactly FZ-01, FZ-02, and FZ-10")

    claimed_rows: set[str] = set()
    for contract_id, contract in contracts.items():
        expected_keys = (
            COMMON_CONTRACT_KEYS | FZ02_CONTRACT_EXTRA_KEYS
            if contract_id == "FZ-02"
            else COMMON_CONTRACT_KEYS
        )
        if not isinstance(contract, dict) or set(contract) != expected_keys:
            fail(f"external_custody_contracts[{contract_id}]: keys mismatch")
        row_ids = contract["rows"]
        if (
            not isinstance(row_ids, list)
            or not row_ids
            or any(not isinstance(row_id, str) for row_id in row_ids)
        ):
            fail(f"external_custody_contracts[{contract_id}].rows must be nonempty")
        for row_id in row_ids:
            if row_id not in rows_by_id:
                fail(
                    f"external_custody_contracts[{contract_id}] names unknown row {row_id}"
                )
            if row_id in claimed_rows:
                fail(f"custody row {row_id} is claimed by more than one contract")
            claimed_rows.add(row_id)

        custody_path = contract["custody_path"]
        manifest_name = contract["registration_manifest"]
        if (
            not isinstance(custody_path, str)
            or Path(custody_path).is_absolute()
            or not isinstance(manifest_name, str)
            or Path(manifest_name).name != manifest_name
        ):
            fail(f"external_custody_contracts[{contract_id}] paths must be relative")
        require_sha256(
            contract["registration_manifest_sha256"],
            f"external_custody_contracts[{contract_id}].registration_manifest_sha256",
        )
        artifacts = validate_hash_mapping(
            contract["artifact_sha256"],
            f"external_custody_contracts[{contract_id}].artifact_sha256",
        )
        in_repo = validate_hash_mapping(
            contract["in_repo_artifact_sha256"],
            f"external_custody_contracts[{contract_id}].in_repo_artifact_sha256",
        )
        state = contract["attestation_state"]
        if state not in ATTESTATION_STATES:
            fail(
                f"external_custody_contracts[{contract_id}] has unknown attestation state"
            )

        for relative_path, expected_hash in in_repo.items():
            source_path = ROOT / relative_path
            if not source_path.is_file():
                fail(f"{contract_id} missing in-repo custody artifact {relative_path}")
            if sha256_file(source_path) != expected_hash:
                fail(
                    f"{contract_id} in-repo custody artifact hash mismatch: {relative_path}"
                )

        for row_id in row_ids:
            row = rows_by_id[row_id]
            attestation = str(row["attestation"] or "").lower()
            if state == "calendar_pending":
                if row["status"] != "frozen_stamped_upgrade_pending":
                    fail(f"{row_id} calendar-pending custody requires pending status")
                if "pending" not in attestation or "bitcoin" not in attestation:
                    fail(f"{row_id} must state that its Bitcoin upgrade is pending")
                if (
                    "complete bitcoin" in attestation
                    or "bitcoin-attested" in attestation
                ):
                    fail(f"{row_id} must not claim a completed Bitcoin attestation")
            elif row["status"] not in {"frozen_attested", "standing_frozen"}:
                fail(f"{row_id} Bitcoin custody requires an attested frozen status")
            elif "bitcoin" not in attestation:
                fail(f"{row_id} must name its Bitcoin attestation")

        if contract_id == "FZ-02":
            for key in (
                "registration_manifest_sha256",
                "custody_erratum_sha256",
                "target_block_sha256",
                "target_payload_sha256",
            ):
                require_sha256(
                    contract[key], f"external_custody_contracts[FZ-02].{key}"
                )
            for key in ("custody_commit", "source_commit"):
                if not isinstance(contract[key], str) or not GIT_COMMIT_RE.fullmatch(
                    contract[key]
                ):
                    fail(
                        f"external_custody_contracts[FZ-02].{key} must be a full commit"
                    )
            custody_time = parse_utc(
                contract["custody_commit_utc"],
                "external_custody_contracts[FZ-02].custody_commit_utc",
            )
            row_time = parse_utc(rows_by_id["FZ-02"]["frozen_utc"], "FZ-02.frozen_utc")
            if row_time != custody_time:
                fail("FZ-02 frozen_utc must equal the corrected custody commit time")
            if contract["custody_commit"] != "1e7d7c73dadeef9aa10ec60061a85cee8426c5b1":
                fail(
                    "FZ-02 custody commit must equal the append-only correction record"
                )
            if contract["source_commit"] != "091658ce585c107a260e7b980352be904d2419b2":
                fail(
                    "FZ-02 source commit must contain the frozen receipt and Lean module"
                )
            if contract["target_file"] not in artifacts:
                fail("FZ-02 target_file must be present in its artifact hash contract")
            for key in ("custody_erratum", "target_file"):
                value = contract[key]
                if not isinstance(value, str) or Path(value).name != value:
                    fail(f"external_custody_contracts[FZ-02].{key} must be a file name")

    return contracts


def validate_retrospective_results(register: dict) -> set[str]:
    results = register.get("retrospective_results")
    if not isinstance(results, list) or not results:
        fail("retrospective_results must be a nonempty list")

    seen_ids: set[str] = set()
    former_reservations: set[str] = set()
    for index, result in enumerate(results):
        where = f"retrospective_results[{index}] ({result.get('id')})"
        if not isinstance(result, dict) or set(result) != RETROSPECTIVE_RESULT_KEYS:
            fail(f"{where}: keys mismatch")
        result_id = result["id"]
        if (
            not isinstance(result_id, str)
            or not result_id
            or result_id in seen_ids
        ):
            fail(f"{where}: id must be a unique nonempty string")
        seen_ids.add(result_id)
        former = result["former_ladder_reservation"]
        if (
            not isinstance(former, str)
            or re.fullmatch(r"FZ-\d{2}", former) is None
            or former in former_reservations
        ):
            fail(f"{where}: former_ladder_reservation must be a unique FZ id")
        former_reservations.add(former)
        if result["status"] != "retrospective_not_evaluable":
            fail(f"{where}: unsupported retrospective status")
        for key in (
            "content",
            "payload_path",
            "comparison_protocol",
            "evidential_boundary",
            "milestone",
        ):
            if not isinstance(result[key], str) or not result[key].strip():
                fail(f"{where}: {key} must be nonempty")
        payload_path = Path(result["payload_path"])
        if payload_path.is_absolute() or ".." in payload_path.parts:
            fail(f"{where}: payload_path must be a repository-relative path")
        require_sha256(result["payload_sha256"], f"{where}.payload_sha256")
        if not isinstance(result["owning_issue"], int):
            fail(f"{where}: owning_issue must identify the closed source issue")

    fz04 = results[0]
    if (
        len(results) != 1
        or fz04["id"] != "RR-506-ALPHA-HVP"
        or fz04["former_ladder_reservation"] != "FZ-04"
        or fz04["owning_issue"] != 506
        or fz04["payload_path"] != FZ04_VERDICT_REL
    ):
        fail("the issue-506 retrospective result binding is malformed")

    verdict = load_json(FZ04_VERDICT_PATH)
    if (
        verdict.get("schema") != "oph.alpha_hvp_class_verdict.v2"
        or verdict.get("issue") != 506
        or verdict.get("row_class")
        != "retrospective_empirical_same_scheme_accounting_audit"
        or verdict.get("verdict")
        != "MULTI_CLASS_NOT_EVALUABLE__ONE_RECORDED_ACCOUNTING_REPLAY_COMPATIBLE"
    ):
        fail("the issue-506 payload has the wrong retrospective verdict identity")
    if verdict.get("scope") != FZ04_SCOPE:
        fail("the issue-506 payload scope differs from the bounded retrospective scope")
    if verdict.get("claim_boundary") != FZ04_CLAIM_BOUNDARY:
        fail("the issue-506 payload claim boundary differs from the bounded statement")

    payload_without_digest = {
        key: value for key, value in verdict.items() if key != "verdict_sha256"
    }
    computed_hash = sha256_bytes(canonical_json_bytes(payload_without_digest))
    reported_hash = verdict.get("verdict_sha256")
    if reported_hash != f"sha256:{computed_hash}":
        fail(
            "the issue-506 payload self-digest does not equal its canonical "
            f"content hash {computed_hash}"
        )
    if fz04["payload_sha256"] != computed_hash:
        fail(
            "the issue-506 retrospective payload hash does not equal the "
            f"canonical content hash {computed_hash}"
        )
    rebuilt_verdict = rebuild_issue506_verdict()
    if verdict != rebuilt_verdict:
        fail(
            "the issue-506 stored payload does not equal the canonical producer replay"
        )
    return former_reservations


def validate(register: dict) -> list[dict]:
    if set(register) != REGISTER_KEYS:
        fail("top-level keys mismatch")
    if register.get("schema") != SCHEMA:
        fail(f"schema must equal {SCHEMA}")
    if register.get("issue") != 607:
        fail("issue must equal 607")
    rows = register.get("rows")
    if not isinstance(rows, list) or not rows:
        fail("rows must be a nonempty list")

    snapshot = load_json(SNAPSHOT_PATH)
    open_issues = {row["number"] for row in snapshot["rows"]}

    seen_ids = [row.get("id") for row in rows]
    if (
        any(not isinstance(row_id, str) for row_id in seen_ids)
        or len(set(seen_ids)) != len(seen_ids)
        or any(re.fullmatch(r"FZ-\d{2}", row_id) is None for row_id in seen_ids)
        or seen_ids != sorted(seen_ids)
    ):
        fail("rows must carry unique ascending FZ identifiers")

    former_reservations = validate_retrospective_results(register)
    overlap = set(seen_ids) & former_reservations
    if overlap:
        fail(
            "retrospective reservations must not appear as prospective ladder "
            f"rows: {sorted(overlap)}"
        )
    allocated = sorted(set(seen_ids) | former_reservations)
    expected_allocated = [f"FZ-{index:02d}" for index in range(1, 11)]
    if allocated != expected_allocated:
        fail(
            "ladder rows and explicitly retired reservations must account for "
            f"{expected_allocated}, got {allocated}"
        )

    rows_by_id: dict[str, dict] = {}
    for index, row in enumerate(rows):
        where = f"rows[{index}] ({row.get('id')})"
        if set(row) != ROW_KEYS:
            fail(f"{where}: keys mismatch")
        if row["status"] not in STATUSES:
            fail(f"{where}: unknown status {row['status']}")
        for key in ("content", "kill_band", "comparison_protocol"):
            if not isinstance(row[key], str) or not row[key].strip():
                fail(f"{where}: {key} must be nonempty")
        if row["status"] in FROZEN_STATUSES:
            for key in ("frozen_utc", "custody", "attestation", "content_sha256"):
                if not isinstance(row[key], str) or not row[key].strip():
                    fail(f"{where}: a frozen row requires {key}")
            parse_utc(row["frozen_utc"], f"{where}.frozen_utc")
        elif row["status"] == "registered_pending_freeze":
            owning = row["owning_issue"]
            if not isinstance(owning, int):
                fail(f"{where}: a pending row requires an owning issue")
            if owning not in open_issues:
                fail(f"{where}: owning issue #{owning} is not open in the snapshot")
            if not isinstance(row["milestone"], str) or not row["milestone"].strip():
                fail(f"{where}: a pending row requires a milestone")
        else:
            if row["owning_issue"] is not None:
                fail(f"{where}: a resource-deferred row cannot retain an open owner")
            if not isinstance(row["milestone"], str) or not row["milestone"].strip():
                fail(f"{where}: a resource-deferred row requires a disposition")
        rows_by_id[row["id"]] = row

    validate_custody_contracts(register, rows_by_id)

    fz02 = rows[1]
    receipt = load_json(FZ02_RECEIPT_PATH)
    live_hash = receipt.get("receipt_sha256")
    if fz02["content_sha256"] != live_hash:
        fail(
            "the FZ-02 content hash does not equal the live angular-multiplet "
            f"receipt hash {live_hash}"
        )
    return rows


def ots_binding(path: Path) -> dict[str, Any]:
    data = path.read_bytes()
    prefix = OTS_DETACHED_HEADER + OTS_SHA256_TAG
    if not data.startswith(prefix) or len(data) < len(prefix) + 32:
        fail(f"{path}: unsupported or malformed detached OpenTimestamps proof")
    digest_start = len(prefix)
    return {
        "file_sha256": data[digest_start : digest_start + 32].hex(),
        "has_pending_calendar": OTS_PENDING_ATTESTATION_TAG in data,
        "has_bitcoin_attestation": OTS_BITCOIN_ATTESTATION_TAG in data,
    }


def verify_ots_binding(path: Path, expected_hash: str, expected_state: str) -> None:
    if not path.is_file():
        fail(f"missing detached OpenTimestamps proof {path}")
    proof = ots_binding(path)
    if proof["file_sha256"] != expected_hash:
        fail(f"OpenTimestamps digest mismatch for {path}")
    if expected_state == "calendar_pending":
        if not proof["has_pending_calendar"] or proof["has_bitcoin_attestation"]:
            fail(f"{path} must be calendar-pending and not yet Bitcoin-attested")
    elif expected_state == "bitcoin_attested":
        if not proof["has_bitcoin_attestation"]:
            fail(f"{path} does not contain a Bitcoin block attestation")
    else:
        fail(f"unknown expected OpenTimestamps state {expected_state}")


def fenced_target_hashes(path: Path) -> tuple[str, str]:
    payload = path.read_bytes()
    start_marker = b"<!-- FZ02-TARGET-BEGIN -->"
    end_marker = b"<!-- FZ02-TARGET-END -->"
    try:
        start = payload.index(start_marker)
        start_line_end = payload.index(b"\n", start) + 1
        end = payload.index(end_marker, start_line_end)
    except ValueError:
        fail(f"{path}: missing FZ-02 target fence")
    end_with_marker = end + len(end_marker)
    if payload[end_with_marker : end_with_marker + 2] == b"\r\n":
        end_with_marker += 2
    elif payload[end_with_marker : end_with_marker + 1] == b"\n":
        end_with_marker += 1
    return (
        sha256_bytes(payload[start:end_with_marker]),
        sha256_bytes(payload[start_line_end:end]),
    )


def verify_external_custody(
    register: dict,
    custody_root: Path = DEFAULT_CUSTODY_ROOT,
) -> dict[str, Any]:
    contracts = register["external_custody_contracts"]
    custody_base = custody_root / "falsification" / "frozen_targets"
    if not custody_base.is_dir():
        return {
            "state": "external_custody_not_present",
            "custody_root": str(custody_root),
            "contracts": {
                contract_id: {
                    "verification": "external_custody_not_present",
                    "attestation_state": contract["attestation_state"],
                }
                for contract_id, contract in contracts.items()
            },
        }

    results: dict[str, dict[str, str]] = {}
    for contract_id, contract in contracts.items():
        directory = custody_root / contract["custody_path"]
        if not directory.is_dir():
            fail(f"{contract_id} external custody directory is missing: {directory}")
        manifest_path = directory / contract["registration_manifest"]
        if not manifest_path.is_file():
            fail(f"{contract_id} registration manifest is missing")
        if sha256_file(manifest_path) != contract["registration_manifest_sha256"]:
            fail(f"{contract_id} registration manifest hash mismatch")
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as error:
            fail(f"{contract_id} registration manifest is invalid JSON: {error}")
        if manifest.get("artifacts") != contract["artifact_sha256"]:
            fail(
                f"{contract_id} registration manifest differs from the in-repo contract"
            )

        state = contract["attestation_state"]
        for filename, expected_hash in contract["artifact_sha256"].items():
            artifact_path = directory / filename
            if not artifact_path.is_file():
                fail(f"{contract_id} custody artifact is missing: {filename}")
            if sha256_file(artifact_path) != expected_hash:
                fail(f"{contract_id} custody artifact hash mismatch: {filename}")
            verify_ots_binding(Path(str(artifact_path) + ".ots"), expected_hash, state)
        verify_ots_binding(
            Path(str(manifest_path) + ".ots"),
            contract["registration_manifest_sha256"],
            state,
        )

        if contract_id == "FZ-02":
            erratum_path = directory / contract["custody_erratum"]
            if not erratum_path.is_file():
                fail("FZ-02 append-only custody erratum is missing")
            if sha256_file(erratum_path) != contract["custody_erratum_sha256"]:
                fail("FZ-02 custody erratum hash mismatch")
            target_path = directory / contract["target_file"]
            block_hash, payload_hash = fenced_target_hashes(target_path)
            if block_hash != contract["target_block_sha256"]:
                fail("FZ-02 fenced target block hash mismatch")
            if payload_hash != contract["target_payload_sha256"]:
                fail("FZ-02 fenced target payload hash mismatch")

        results[contract_id] = {
            "verification": "verified",
            "attestation_state": state,
        }
    return {
        "state": "verified",
        "custody_root": str(custody_root),
        "contracts": results,
    }


def render(register: dict, rows: list[dict]) -> str:
    lines: list[str] = []
    lines.append("# The frozen-prediction ladder")
    lines.append("")
    lines.append(
        "Generated by `tools/build_fz_registry.py` from"
        " `claims/frozen_prediction_register.json`; edit the JSON, then"
        " regenerate. The standing register was established under issue #607."
    )
    lines.append("")
    lines.append(register["policy"])
    lines.append("")
    lines.append("| Freeze | Content | Status | Frozen (UTC) | Owner | Kill band |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for row in rows:
        owner = (
            f"[#{row['owning_issue']}](https://github.com/FloatingPragma/observer-patch-holography/issues/{row['owning_issue']})"
            f" ({row['milestone']})"
            if row["owning_issue"] is not None
            else row["milestone"]
        )
        if row["frozen_utc"]:
            frozen = row["frozen_utc"]
        elif row["status"] == "resource_deferred":
            frozen = "not registered"
        else:
            frozen = "to freeze"
        lines.append(
            f"| {row['id']} | {row['content']} | {row['status']} | {frozen} |"
            f" {owner} | {row['kill_band']} |"
        )
    lines.append("")
    lines.append("## Retrospective results outside the ladder")
    lines.append("")
    lines.append(
        "These records were evaluated after their comparison inputs were known."
        " They are not freezes, ladder rungs, predictions, or evidence from a"
        " prospective test. A former reservation remains visible only to make"
        " the bookkeeping transition traceable."
    )
    lines.append("")
    lines.append(
        "| Record | Former reservation | Result | Status | Source | Payload hash |"
    )
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for result in register["retrospective_results"]:
        owner = (
            f"[#{result['owning_issue']}](https://github.com/FloatingPragma/"
            "observer-patch-holography/issues/"
            f"{result['owning_issue']}) ({result['milestone']})"
        )
        lines.append(
            f"| {result['id']} | {result['former_ladder_reservation']} |"
            f" {result['content']} | {result['status']} | {owner} |"
            f" `{result['payload_sha256']}` |"
        )
    lines.append("")
    for result in register["retrospective_results"]:
        lines.append(
            f"- **{result['id']} protocol**: {result['comparison_protocol']}"
        )
        lines.append(
            f"  Evidential boundary: {result['evidential_boundary']}"
        )
        lines.append(f"  Payload: `{result['payload_path']}`.")
    lines.append("")
    lines.append("## Custody and attestation")
    lines.append("")
    for row in rows:
        if row["status"] in FROZEN_STATUSES:
            lines.append(f"- **{row['id']}**: {row['custody']} {row['attestation']}.")
            lines.append(
                f"  Content hash: `{row['content_sha256']}`."
                f" Comparison protocol: {row['comparison_protocol']}"
            )
    lines.append("")
    lines.append("## Custody verification contracts")
    lines.append("")
    lines.append("| Contract | Rows | Required attestation state | External custody |")
    lines.append("| --- | --- | --- | --- |")
    state_labels = {
        "bitcoin_attested": "Bitcoin block attestation present",
        "calendar_pending": "calendar commitments present; Bitcoin upgrade pending",
    }
    for contract_id, contract in register["external_custody_contracts"].items():
        lines.append(
            f"| {contract_id} | {', '.join(contract['rows'])} |"
            f" {state_labels[contract['attestation_state']]} |"
            f" `{contract['custody_path']}` |"
        )
    lines.append("")
    fz02_contract = register["external_custody_contracts"]["FZ-02"]
    lines.append(
        "FZ-02 is bound to oph-meta custody commit"
        f" `{fz02_contract['custody_commit']}` at"
        f" `{fz02_contract['custody_commit_utc']}` and source commit"
        f" `{fz02_contract['source_commit']}`. Its append-only erratum corrects"
        " the original timestamp, source-commit, and whole-file-versus-fenced-"
        "block hash metadata without modifying any stamped artifact."
    )
    lines.append("")
    lines.append(
        "The validator always checks the committed in-repo hash contracts. In"
        " the coordinated oph-meta workspace it also resolves the sibling"
        " custody directories, recomputes every manifest and artifact hash,"
        " checks each detached `.ots` digest, and distinguishes pending calendar"
        " commitments from Bitcoin block attestations. In an isolated source"
        " clone it reports `external_custody_not_present`; that classification"
        " is clean-clone-safe but is not an external-artifact verification."
        " The local structural check does not contact a Bitcoin node; independent"
        " chain verification remains the job of `ots verify` after an upgrade."
    )
    lines.append("")
    lines.append(
        "Pending rows freeze at their milestones, before their comparison data"
    )
    lines.append(
        "is examined; the register validation requires each pending row to name"
    )
    lines.append("a live owning issue and fails closed otherwise.")
    lines.append(
        "Retrospective results are validated and rendered in their separate"
        " section. Their former reservations do not occur in the ladder table."
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    register = load_json(REGISTER_PATH)
    rows = validate(register)
    custody = verify_external_custody(register)
    surface = render(register, rows)
    if args.check:
        committed = (
            SURFACE_PATH.read_text(encoding="utf-8") if SURFACE_PATH.is_file() else ""
        )
        if committed != surface:
            print(
                "frozen-prediction register: docs/FROZEN_PREDICTION_LADDER.md is "
                "stale; run python tools/build_fz_registry.py",
                file=sys.stderr,
            )
            return 1
        print(
            "frozen-prediction register: external custody "
            f"{custody['state']} at {custody['custody_root']}"
        )
        print("frozen-prediction register: surface is current")
        return 0
    SURFACE_PATH.write_text(surface, encoding="utf-8", newline="\n")
    print(
        "frozen-prediction register: external custody "
        f"{custody['state']} at {custody['custody_root']}"
    )
    print(f"frozen-prediction register: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
