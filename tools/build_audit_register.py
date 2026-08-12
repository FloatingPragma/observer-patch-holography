#!/usr/bin/env python3
"""Validate and render the append-only V3 audit register."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import strict_json

ROOT = Path(__file__).resolve().parents[1]
REGISTER_PATH = ROOT / "tracking" / "audit_register.json"
LEDGER_PATH = ROOT / "tracking" / "observation_ledger.json"
SURFACE_PATH = ROOT / "docs" / "AUDIT_REGISTER.md"

SCHEMA = "oph.audit_register.v1"
ISSUE = 738
AUDIT_ID_RE = re.compile(r"^AUD-[A-Z0-9-]+$")
COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
DATE_RE = re.compile(r"^20[0-9]{2}-[01][0-9]-[0-3][0-9]$")
SEVERITIES = {"critical", "high", "medium", "low", "boundary"}
DISPOSITIONS = {"fixed", "tracked_open", "verified_boundary"}

TOP_KEYS = {"schema", "issue", "policy", "records"}
RECORD_KEYS = {
    "id",
    "started_on",
    "completed_on",
    "baseline_commit",
    "audited_head",
    "audit_class",
    "reviewers",
    "scope",
    "reviewed_observation_rows",
    "attained_rows_reviewed",
    "qualifies_for",
    "does_not_qualify_for",
    "findings",
    "evidence",
    "limitations",
}
REVIEWER_KEYS = {"name", "task", "role"}
FINDING_KEYS = {"id", "severity", "disposition", "summary", "owner_issues"}


def fail(message: str) -> None:
    raise SystemExit(f"audit register: {message}")


def _display_path(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def load_json(path: Path) -> dict:
    try:
        return strict_json.load(path)
    except FileNotFoundError:
        fail(f"missing input {_display_path(path)}")
    except (json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"invalid JSON in {_display_path(path)}: {error}")
    raise AssertionError("unreachable")


def _strings(where: str, value: object, *, allow_empty: bool = False) -> list[str]:
    if not isinstance(value, list) or (not value and not allow_empty):
        fail(
            f"{where} must be a {'possibly empty ' if allow_empty else 'nonempty '}list"
        )
    if any(not isinstance(item, str) or not item.strip() for item in value):
        fail(f"{where} entries must be nonempty strings")
    if len(value) != len(set(value)):
        fail(f"{where} must be duplicate-free")
    return value


def validate(data: dict) -> list[dict]:
    if not isinstance(data, dict) or set(data) != TOP_KEYS:
        fail(f"top-level keys must equal {sorted(TOP_KEYS)}")
    if data["schema"] != SCHEMA or data["issue"] != ISSUE:
        fail(f"schema and issue must equal {SCHEMA} and {ISSUE}")
    if not isinstance(data["policy"], str) or not data["policy"].strip():
        fail("policy must be a nonempty string")
    ledger = load_json(LEDGER_PATH)
    ledger_ids = {
        row.get("id") for row in ledger.get("rows", []) if isinstance(row, dict)
    }
    records = data["records"]
    if not isinstance(records, list) or not records:
        fail("records must be a nonempty list")
    seen: set[str] = set()
    for index, record in enumerate(records):
        where = f"records[{index}]"
        if not isinstance(record, dict) or set(record) != RECORD_KEYS:
            fail(f"{where}: keys must equal {sorted(RECORD_KEYS)}")
        audit_id = record["id"]
        if not isinstance(audit_id, str) or not AUDIT_ID_RE.fullmatch(audit_id):
            fail(f"{where}: malformed audit id")
        if audit_id in seen:
            fail(f"duplicate audit id {audit_id}")
        seen.add(audit_id)
        where = audit_id
        for key in ("started_on", "completed_on"):
            if not isinstance(record[key], str) or not DATE_RE.fullmatch(record[key]):
                fail(f"{where}: malformed {key}")
        for key in ("baseline_commit", "audited_head"):
            if not isinstance(record[key], str) or not COMMIT_RE.fullmatch(record[key]):
                fail(f"{where}: malformed {key}")
        for key in ("audit_class", "scope"):
            if not isinstance(record[key], str) or not record[key].strip():
                fail(f"{where}: {key} must be nonempty")
        reviewers = record["reviewers"]
        if not isinstance(reviewers, list) or not reviewers:
            fail(f"{where}: reviewers must be nonempty")
        tasks: set[str] = set()
        for reviewer in reviewers:
            if not isinstance(reviewer, dict) or set(reviewer) != REVIEWER_KEYS:
                fail(f"{where}: malformed reviewer")
            if any(
                not isinstance(reviewer[key], str) or not reviewer[key].strip()
                for key in REVIEWER_KEYS
            ):
                fail(f"{where}: reviewer fields must be nonempty")
            if reviewer["task"] in tasks:
                fail(f"{where}: reviewer tasks must be distinct")
            tasks.add(reviewer["task"])
        reviewed = _strings(
            f"{where}.reviewed_observation_rows", record["reviewed_observation_rows"]
        )
        attained = _strings(
            f"{where}.attained_rows_reviewed", record["attained_rows_reviewed"]
        )
        if set(reviewed) - ledger_ids or set(attained) - set(reviewed):
            fail(f"{where}: reviewed or attained row is outside the ledger/scope")
        _strings(f"{where}.qualifies_for", record["qualifies_for"])
        _strings(f"{where}.does_not_qualify_for", record["does_not_qualify_for"])
        _strings(f"{where}.limitations", record["limitations"])
        evidence = _strings(f"{where}.evidence", record["evidence"])
        for path in evidence:
            if (
                path.startswith("/")
                or ".." in path.split("/")
                or not (ROOT / path).is_file()
            ):
                fail(f"{where}: evidence path missing or non-relative: {path}")
        findings = record["findings"]
        if not isinstance(findings, list) or not findings:
            fail(f"{where}: findings must be nonempty")
        finding_ids: set[str] = set()
        for finding in findings:
            if not isinstance(finding, dict) or set(finding) != FINDING_KEYS:
                fail(f"{where}: malformed finding")
            finding_id = finding["id"]
            if not isinstance(finding_id, str) or finding_id in finding_ids:
                fail(f"{where}: finding ids must be unique strings")
            finding_ids.add(finding_id)
            if finding["severity"] not in SEVERITIES:
                fail(f"{where}.{finding_id}: invalid severity")
            if finding["disposition"] not in DISPOSITIONS:
                fail(f"{where}.{finding_id}: invalid disposition")
            if (
                not isinstance(finding["summary"], str)
                or not finding["summary"].strip()
            ):
                fail(f"{where}.{finding_id}: summary must be nonempty")
            owners = finding["owner_issues"]
            if (
                not isinstance(owners, list)
                or any(
                    not isinstance(issue, int) or isinstance(issue, bool)
                    for issue in owners
                )
                or len(owners) != len(set(owners))
            ):
                fail(f"{where}.{finding_id}: owner_issues must be unique integers")
    return records


def render(data: dict, records: list[dict]) -> str:
    lines = [
        "# Audit register",
        "",
        "Generated by `tools/build_audit_register.py` from `tracking/audit_register.json`; edit the JSON and regenerate. Standing custody belongs to [issue #738](https://github.com/FloatingPragma/observer-patch-holography/issues/738).",
        "",
        data["policy"],
    ]
    for record in records:
        lines.extend(
            [
                "",
                f"## {record['id']}",
                "",
                f"Audit window: `{record['started_on']}` through `{record['completed_on']}`. Baseline `{record['baseline_commit']}`; audited head `{record['audited_head']}`.",
                "",
                f"Class: {record['audit_class']}",
                "",
                record["scope"],
                "",
                "Reviewers:",
                "",
                *[
                    f"- **{item['name']}** (`{item['task']}`): {item['role']}"
                    for item in record["reviewers"]
                ],
                "",
                "| Finding | Severity | Disposition | Owners | Summary |",
                "| --- | --- | --- | --- | --- |",
                *[
                    f"| {item['id']} | {item['severity']} | {item['disposition']} | "
                    + (
                        ", ".join(f"#{issue}" for issue in item["owner_issues"])
                        or "none"
                    )
                    + f" | {item['summary']} |"
                    for item in record["findings"]
                ],
                "",
                "Qualifies for: " + ", ".join(record["qualifies_for"]) + ".",
                "",
                "Does not qualify for: "
                + ", ".join(record["does_not_qualify_for"])
                + ".",
                "",
                "Limitations:",
                "",
                *[f"- {item}" for item in record["limitations"]],
                "",
                "Evidence:",
                "",
                *[f"- `{path}`" for path in record["evidence"]],
            ]
        )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = load_json(REGISTER_PATH)
    records = validate(data)
    output = render(data, records).encode("utf-8")
    if args.check:
        if not SURFACE_PATH.is_file() or SURFACE_PATH.read_bytes() != output:
            fail("generated surface is stale")
        print(f"audit register: {len(records)} record(s), surface current")
        return 0
    SURFACE_PATH.write_bytes(output)
    print(f"audit register: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
