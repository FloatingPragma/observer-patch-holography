#!/usr/bin/env python3
"""Verify the archived OPH-FPE local-domain receipt family without simulator imports."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "manifest.json"

PIN_FIELDS = {
    "arrays": "arrays_sha256",
    "receipt": "receipt_sha256",
    "stage2_receipt": "stage2_receipt_sha256",
    "stage3_receipt": "stage3_receipt_sha256",
    "stage4_receipt": "stage4_receipt_sha256",
    "defect_sector_receipt": "defect_sector_receipt_sha256",
    "clock_unit_verdict": "clock_unit_verdict_sha256",
    "classical_realization_receipt": "classical_realization_receipt_sha256",
    "matter_attachment_receipt": "matter_attachment_receipt_sha256",
    "source_gap_receipt": "source_gap_receipt_sha256",
}

EXPECTED_SCHEMAS = {
    "stage1_receipt.json": "oph.local-domain-stage1.v1",
    "stage2_receipt.json": "oph.local-domain-stage2.v1",
    "stage3_receipt.json": "oph.local-domain-stage3.v1",
    "stage4_receipt.json": "oph.local-domain-stage4.v1",
    "defect_sector_receipt.json": "oph.local-domain-defect-sector-spectra.v1",
    "clock_unit_verdict.json": "oph.local-domain-clock-unit-verdict.v1",
    "classical_realization_receipt.json":
        "oph.local-domain-classical-realization.v1",
    "matter_attachment_receipt.json": "oph.local-domain-matter-attachment.v1",
    "source_gap_receipt.json": "oph.source-clock-gap.v1",
}


class VerificationError(RuntimeError):
    """Raised when archive custody or a theorem-level boundary changes."""


def _reject_constant(value: str) -> None:
    raise VerificationError(f"non-finite JSON constant: {value}")


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise VerificationError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_no_duplicates,
        parse_constant=_reject_constant,
    )
    if not isinstance(value, dict):
        raise VerificationError(f"{path.name} is not a JSON object")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


def _assert_finite(value: Any, where: str = "$") -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise VerificationError(f"non-finite number at {where}")
    if isinstance(value, dict):
        for key, child in value.items():
            _assert_finite(child, f"{where}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _assert_finite(child, f"{where}[{index}]")


def verify() -> dict[str, Any]:
    manifest = _load_json(MANIFEST)
    if manifest.get("schema") != "oph.local-domain-stage1.manifest.v1":
        raise VerificationError("manifest schema changed")

    verified_files: list[str] = []
    documents: dict[str, dict[str, Any]] = {}
    for name_field, hash_field in PIN_FIELDS.items():
        name = manifest.get(name_field)
        expected_hash = manifest.get(hash_field)
        if not isinstance(name, str) or Path(name).name != name:
            raise VerificationError(f"unsafe or missing manifest path: {name_field}")
        if not isinstance(expected_hash, str) or not expected_hash.startswith("sha256:"):
            raise VerificationError(f"missing manifest digest: {hash_field}")
        path = ROOT / name
        if not path.is_file():
            raise VerificationError(f"missing archived file: {name}")
        if _sha256(path) != expected_hash:
            raise VerificationError(f"archive digest mismatch: {name}")
        verified_files.append(name)
        if path.suffix == ".json":
            document = _load_json(path)
            _assert_finite(document)
            documents[name] = document

    for name, schema in EXPECTED_SCHEMAS.items():
        if documents[name].get("schema") != schema:
            raise VerificationError(f"schema changed: {name}")

    stage1 = documents["stage1_receipt.json"]
    stage4 = documents["stage4_receipt.json"]
    if stage1.get("event_count") != 2304:
        raise VerificationError("stage1 event count changed")
    if stage1.get("verdict") != "ATTAINED":
        raise VerificationError("stage1 finite verdict changed")
    if stage1.get("physical_promotion_allowed") is not False:
        raise VerificationError("stage1 physical-promotion boundary changed")
    if stage4.get("bundle_resolution", {}).get("passed") is not True:
        raise VerificationError("stage4 bundle resolution failed")
    replay = stage4.get("producer_semantic_replay", {})
    if replay.get("all_semantic_exact") is not True:
        raise VerificationError("stage4 semantic replay changed")
    if replay.get("producer_independence") is not False:
        raise VerificationError("stage4 producer-independence boundary changed")

    for name, document in documents.items():
        if "physical_promotion_allowed" in document:
            if document["physical_promotion_allowed"] is not False:
                raise VerificationError(
                    f"physical-promotion boundary changed: {name}"
                )

    return {
        "status": "VERIFIED_RER_LOCAL_DOMAIN_ARCHIVE",
        "manifest_sha256": _sha256(MANIFEST),
        "archived_files_verified": len(verified_files),
        "event_count": stage1["event_count"],
        "stage1_verdict": stage1["verdict"],
        "stage4_bundle_passed": True,
        "producer_independence_claimed": False,
        "physical_promotion_allowed": False,
    }


def main() -> int:
    try:
        result = verify()
    except (OSError, UnicodeError, json.JSONDecodeError, VerificationError) as exc:
        print(json.dumps({"status": "FAIL", "reason": str(exc)}, sort_keys=True))
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
