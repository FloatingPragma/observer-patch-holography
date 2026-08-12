"""Build and validate the V3 emergent-instrument register surface (issue #738).

The machine-readable register is ``claims/emergent_instrument_register.json``.
This tool validates it fail-closed and renders
``docs/INSTRUMENT_REGISTER_V3.md``; ``--check`` fails when the committed page
differs byte for byte from the render.

The register holds preregistered simulation instruments. Lane #737 owns the
instruments; the standing custody lane #738 owns this register. It is a
separate surface from the frozen-prediction ladder: the ladder is reserved for
physical predictions against external data, and simulation instruments never
enter it.

Fail-closed rules: row ids match INS-<two digits>, ascend strictly, and start
at INS-01; every status comes from the seven-value enum; every ledger row id
names a committed observation-ledger row; every owning issue is a V3 lane
issue in 728..738. A SPECIFIED row carries no freeze artifacts, no freeze
time, and no verdict receipts. A FROZEN or RUNNING row carries a parseable,
non-future UTC freeze time, at least one freeze artifact, a decision rule
naming all three verdict labels, and no verdict receipts. A REPLICATED,
FAILED, or INCONCLUSIVE row additionally carries at least one verdict
receipt. A VOID row issues no verdict and carries freeze artifacts exactly
when it carries a freeze time. Artifact entries pin a relative path and a
SHA-256 digest; an entry whose path resolves inside this repository is hashed
against its pinned digest, and an entry outside it is a recorded pointer
whose bytes are verified in the owning repository's custody checkout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTER_PATH = ROOT / "claims" / "emergent_instrument_register.json"
LEDGER_PATH = ROOT / "tracking" / "observation_ledger.json"
SURFACE_PATH = ROOT / "docs" / "INSTRUMENT_REGISTER_V3.md"

SCHEMA = "oph.emergent_instrument_register.v1"
ISSUE = 738
GENERATED_SURFACE = "docs/INSTRUMENT_REGISTER_V3.md"
REPO_URL = "https://github.com/muellerberndt/reverse-engineering-reality"

STATUSES = (
    "SPECIFIED",
    "FROZEN",
    "RUNNING",
    "REPLICATED",
    "FAILED",
    "INCONCLUSIVE",
    "VOID",
)
FROZEN_PRE_RUN_STATUSES = {"FROZEN", "RUNNING"}
VERDICT_STATUSES = {"REPLICATED", "FAILED", "INCONCLUSIVE"}
VERDICT_LABELS = ("REPLICATED", "FAILED", "INCONCLUSIVE")
LANE_MIN = 728
LANE_MAX = 738

ID_PATTERN = re.compile(r"^INS-\d{2}$")
LEDGER_ID_PATTERN = re.compile(r"^OL-[A-I][1-9]$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
BANNED_CHARACTERS = ("—", "–")

REGISTER_KEYS = {"schema", "issue", "generated_surface", "policy", "rows"}
ROW_KEYS = {
    "id",
    "title",
    "owning_issue",
    "ledger_row",
    "status",
    "spec_pointer",
    "decision_rule",
    "seeds_policy",
    "controls",
    "freeze_artifacts",
    "frozen_utc",
    "verdict_receipts",
}
ARTIFACT_KEYS = {"path", "sha256"}

# Paths under these prefixes are recorded pointers into external
# repositories; every other relative path must exist in this repository
# and match its pinned digest.
EXTERNAL_PATH_PREFIXES = ("oph-physics-sim/",)

STATUS_MEANING = {
    "SPECIFIED": (
        "A design specification is recorded and no freeze has happened. The"
        " row carries no freeze artifacts, no freeze time, and no verdict"
        " receipts; the decision rule may be a pointer into the"
        " specification."
    ),
    "FROZEN": (
        "The decision rule, seeds policy, and controls are pinned before any"
        " run. The row carries the verbatim decision rule, a non-future UTC"
        " freeze time, and hashed freeze artifacts; no verdict receipt"
        " exists."
    ),
    "RUNNING": (
        "The frozen instrument executes under its pinned configuration. The"
        " row carries its freeze fields and no verdict receipts."
    ),
    "REPLICATED": (
        "The frozen decision rule returned its positive verdict. The row"
        " carries freeze artifacts and verdict receipts."
    ),
    "FAILED": (
        "The frozen decision rule returned its negative verdict. The row"
        " carries freeze artifacts and verdict receipts, and the bound"
        " ledger row's emergent rung is demoted with equal prominence."
    ),
    "INCONCLUSIVE": (
        "The frozen decision rule returned no verdict. The row carries"
        " freeze artifacts and verdict receipts, and the bound ledger row is"
        " neither promoted nor demoted by this instrument."
    ),
    "VOID": (
        "The registration is withdrawn or superseded. A VOID row issues no"
        " verdict, and it carries freeze artifacts exactly when it carries a"
        " freeze time."
    ),
}


def fail(message: str) -> None:
    raise SystemExit(f"instrument register: {message}")


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        fail(f"missing input {path.relative_to(ROOT)}")
    except json.JSONDecodeError as error:
        fail(f"invalid JSON in {path.relative_to(ROOT)}: {error}")
    raise AssertionError("unreachable")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_utc(value: object, where: str) -> datetime:
    if not isinstance(value, str) or not value.endswith("Z"):
        fail(f"{where} must be an ISO-8601 UTC timestamp ending in Z")
    try:
        parsed = datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError:
        fail(f"{where} is not a valid ISO-8601 UTC timestamp")
    if parsed > datetime.now(timezone.utc):
        fail(f"{where} cannot be in the future")
    return parsed


def _clean_prose(where: str, field: str, value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        fail(f"{where}: {field} must be a nonempty string")
    for character in BANNED_CHARACTERS:
        if character in value:
            fail(f"{where}: {field} carries a banned dash character")
    return value


def validate_artifact_entries(entries: object, where: str) -> list[dict]:
    if not isinstance(entries, list):
        fail(f"{where} must be a list")
    seen_paths: set[str] = set()
    for index, entry in enumerate(entries):
        entry_where = f"{where}[{index}]"
        if not isinstance(entry, dict) or set(entry) != ARTIFACT_KEYS:
            fail(f"{entry_where} must be an object with exactly path, sha256")
        path = entry["path"]
        if not isinstance(path, str) or not path:
            fail(f"{entry_where}: path must be a nonempty string")
        if path.startswith("/") or ".." in path.split("/"):
            fail(f"{entry_where}: path {path} must be relative")
        if path in seen_paths:
            fail(f"{where}: duplicate artifact path {path}")
        seen_paths.add(path)
        digest = entry["sha256"]
        if not isinstance(digest, str) or not SHA256_RE.fullmatch(digest):
            fail(f"{entry_where}: sha256 must be a lowercase SHA-256 hex digest")
        if any(path.startswith(prefix) for prefix in EXTERNAL_PATH_PREFIXES):
            # A recorded pointer into the simulator repository: the digest
            # is custody data and is checked there, not here.
            continue
        candidate = ROOT / path
        if not candidate.is_file():
            fail(f"{entry_where}: in-repo artifact path {path} does not exist")
        if sha256_file(candidate) != digest:
            fail(f"{entry_where}: in-repo artifact hash mismatch for {path}")
    return entries


def load_ledger_row_ids() -> set[str]:
    ledger = load_json(LEDGER_PATH)
    rows = ledger.get("rows") if isinstance(ledger, dict) else None
    if not isinstance(rows, list) or not rows:
        fail("observation ledger rows are unavailable for cross-checking")
    ids: set[str] = set()
    for row in rows:
        if isinstance(row, dict) and isinstance(row.get("id"), str):
            ids.add(row["id"])
    if not ids:
        fail("observation ledger carries no row ids")
    return ids


def validate(register: dict, ledger_row_ids: set[str] | None = None) -> list[dict]:
    if not isinstance(register, dict):
        fail("register must be an object")
    if set(register) != REGISTER_KEYS:
        fail(
            "top-level keys must be exactly schema, issue, generated_surface,"
            " policy, rows"
        )
    if register["schema"] != SCHEMA:
        fail(f"schema must equal {SCHEMA}")
    if register["issue"] != ISSUE:
        fail(f"issue must equal {ISSUE}")
    if register["generated_surface"] != GENERATED_SURFACE:
        fail(f"generated_surface must equal {GENERATED_SURFACE}")
    _clean_prose("register", "policy", register["policy"])
    rows = register["rows"]
    if not isinstance(rows, list) or not rows:
        fail("rows must be a nonempty list")
    if ledger_row_ids is None:
        ledger_row_ids = load_ledger_row_ids()

    previous_number = 0
    for index, row in enumerate(rows):
        where = f"rows[{index}]"
        if not isinstance(row, dict):
            fail(f"{where}: row must be an object")
        if set(row) != ROW_KEYS:
            missing = ROW_KEYS - set(row)
            extra = set(row) - ROW_KEYS
            fail(
                f"{where}: keys mismatch "
                f"(missing {sorted(missing)}, extra {sorted(extra)})"
            )
        row_id = row["id"]
        if not isinstance(row_id, str) or not ID_PATTERN.match(row_id):
            fail(f"{where}: id must match INS-<two digits>")
        number = int(row_id[4:])
        if index == 0 and number != 1:
            fail("the first row must be INS-01")
        if number <= previous_number:
            fail(f"{where}: row ids must ascend strictly, found {row_id}")
        previous_number = number
        where = row_id

        _clean_prose(where, "title", row["title"])
        _clean_prose(where, "spec_pointer", row["spec_pointer"])
        _clean_prose(where, "decision_rule", row["decision_rule"])
        _clean_prose(where, "seeds_policy", row["seeds_policy"])

        owning_issue = row["owning_issue"]
        if not isinstance(owning_issue, int) or isinstance(owning_issue, bool):
            fail(f"{where}: owning_issue must be an integer")
        if not LANE_MIN <= owning_issue <= LANE_MAX:
            fail(f"{where}: owning_issue must lie in {LANE_MIN}..{LANE_MAX}")

        ledger_row = row["ledger_row"]
        if not isinstance(ledger_row, str) or not LEDGER_ID_PATTERN.match(
            ledger_row
        ):
            fail(f"{where}: ledger_row must match OL-<A..I><digit>")
        if ledger_row not in ledger_row_ids:
            fail(f"{where}: ledger_row {ledger_row} is not on the observation ledger")

        status = row["status"]
        if status not in STATUSES:
            fail(f"{where}: status must be one of {STATUSES}")

        controls = row["controls"]
        if not isinstance(controls, list) or not controls:
            fail(f"{where}: controls must be a nonempty list")
        if len(controls) != len(set(controls)):
            fail(f"{where}: controls must be duplicate-free")
        for position, control in enumerate(controls):
            _clean_prose(where, f"controls[{position}]", control)

        freeze_artifacts = validate_artifact_entries(
            row["freeze_artifacts"], f"{where}.freeze_artifacts"
        )
        verdict_receipts = validate_artifact_entries(
            row["verdict_receipts"], f"{where}.verdict_receipts"
        )

        frozen_utc = row["frozen_utc"]
        if frozen_utc is not None:
            parse_utc(frozen_utc, f"{where}.frozen_utc")
        if frozen_utc is None and freeze_artifacts:
            fail(f"{where}: an unfrozen row carries no freeze artifacts")

        if status == "SPECIFIED":
            if frozen_utc is not None:
                fail(f"{where}: a SPECIFIED row carries no freeze time")
            if freeze_artifacts:
                fail(f"{where}: a SPECIFIED row carries no freeze artifacts")
            if verdict_receipts:
                fail(f"{where}: a SPECIFIED row carries no verdict receipts")
        elif status in FROZEN_PRE_RUN_STATUSES or status in VERDICT_STATUSES:
            if frozen_utc is None:
                fail(f"{where}: a {status} row requires a freeze time")
            if not freeze_artifacts:
                fail(f"{where}: a {status} row requires freeze artifacts")
            for label in VERDICT_LABELS:
                if label not in row["decision_rule"]:
                    fail(
                        f"{where}: a frozen decision rule must name the"
                        f" verdict label {label}"
                    )
            if status in FROZEN_PRE_RUN_STATUSES and verdict_receipts:
                fail(f"{where}: a {status} row carries no verdict receipts")
            if status in VERDICT_STATUSES and not verdict_receipts:
                fail(f"{where}: a {status} row requires verdict receipts")
        else:  # VOID
            if verdict_receipts:
                fail(f"{where}: a VOID row carries no verdict receipts")
            if frozen_utc is not None and not freeze_artifacts:
                fail(
                    f"{where}: a VOID row with a freeze time keeps its freeze"
                    " artifacts"
                )
    return rows


def _issue_link(number: int) -> str:
    return f"[#{number}]({REPO_URL}/issues/{number})"


def _artifact_lines(entries: list[dict], label: str, empty_note: str) -> list[str]:
    if not entries:
        return [f"- {label}: {empty_note}."]
    lines = [f"- {label}:"]
    for entry in entries:
        lines.append(f"  - `{entry['path']}` sha256 `{entry['sha256']}`")
    return lines


def render(register: dict, rows: list[dict]) -> str:
    lines: list[str] = []
    lines.append("# V3 Emergent-Instrument Register")
    lines.append("")
    lines.append(
        "Generated by `tools/build_instrument_register.py` from"
        " `claims/emergent_instrument_register.json`; edit the JSON, then"
        " regenerate."
    )
    lines.append("")
    lines.append(
        f"One row per preregistered simulation instrument of the emergent"
        f" adequacy program. Lane {_issue_link(737)} owns the instruments;"
        f" the standing custody lane {_issue_link(738)} owns this register."
        f" Each instrument binds to exactly one row of the observation ledger"
        f" (`docs/OBSERVATION_LEDGER_V3.md`), and its verdict promotes or"
        f" blocks that row's emergent rung."
    )
    lines.append("")
    lines.append("## Separation from the frozen-prediction ladder")
    lines.append("")
    lines.append(
        "The frozen-prediction ladder (`docs/FROZEN_PREDICTION_LADDER.md`) is"
        " reserved for physical predictions against external data. Simulation"
        " instruments live on this register and never enter that ladder: an"
        " instrument verdict concerns what simulated observers' records"
        " exhibit inside the architecture and makes no physical prediction."
    )
    lines.append("")
    lines.append(
        "Freeze artifacts pin `{path, sha256}` pairs. Paths under"
        " `oph-physics-sim/` are recorded pointers into the simulator"
        " repository, checked by its own custody; every other relative path"
        " must exist in this repository and match its pinned digest, and the"
        " validator fails closed otherwise."
    )
    lines.append("")
    lines.append("## Policy")
    lines.append("")
    lines.append(register["policy"])
    lines.append("")
    lines.append("## Status enum")
    lines.append("")
    for status in STATUSES:
        lines.append(f"- **{status}**: {STATUS_MEANING[status]}")
    lines.append("")
    lines.append("## Coherence rules (validator, fail-closed)")
    lines.append("")
    lines.append(
        "A SPECIFIED row carries no freeze time, no freeze artifacts, and no"
        " verdict receipts. A FROZEN or RUNNING row carries a non-future UTC"
        " freeze time, at least one hashed freeze artifact, a decision rule"
        " naming REPLICATED, FAILED, and INCONCLUSIVE, and no verdict"
        " receipts. A REPLICATED, FAILED, or INCONCLUSIVE row additionally"
        " carries at least one verdict receipt. A VOID row issues no verdict"
        " and carries freeze artifacts exactly when it carries a freeze time."
        " An artifact entry whose path resolves inside this repository is"
        " hashed against its pinned digest; an entry outside it is a recorded"
        " pointer whose bytes are verified in the owning repository's custody"
        " checkout."
    )
    lines.append("")
    lines.append("## Instruments")
    lines.append("")
    lines.append("| Row | Instrument | Ledger row | Lane | Status | Frozen (UTC) |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for row in rows:
        frozen = row["frozen_utc"] if row["frozen_utc"] is not None else "unfrozen"
        lines.append(
            f"| {row['id']} | {row['title']} | {row['ledger_row']} |"
            f" {_issue_link(row['owning_issue'])} | `{row['status']}` |"
            f" {frozen} |"
        )

    for row in rows:
        lines.append("")
        lines.append(f"### {row['id']} {row['title']}")
        lines.append("")
        lines.append(
            f"- Ledger row {row['ledger_row']}; owning lane"
            f" {_issue_link(row['owning_issue'])}; status `{row['status']}`."
        )
        lines.append(f"- Specification: `{row['spec_pointer']}`.")
        lines.append(f"- Decision rule: {row['decision_rule']}")
        lines.append(f"- Seeds policy: {row['seeds_policy']}")
        lines.append("- Controls:")
        for control in row["controls"]:
            lines.append(f"  - {control}")
        lines.extend(
            _artifact_lines(
                row["freeze_artifacts"], "Freeze artifacts", "none (unfrozen row)"
            )
        )
        frozen = row["frozen_utc"] if row["frozen_utc"] is not None else "none"
        lines.append(f"- Frozen (UTC): {frozen}.")
        lines.extend(
            _artifact_lines(
                row["verdict_receipts"], "Verdict receipts", "none (no run verdict)"
            )
        )

    counts = {status: 0 for status in STATUSES}
    for row in rows:
        counts[row["status"]] += 1
    populated = ", ".join(
        f"{counts[status]} {status}" for status in STATUSES if counts[status]
    )
    noun = "instrument" if len(rows) == 1 else "instruments"
    lines.append("")
    lines.append(f"Totals: {len(rows)} {noun}. Status: {populated}.")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when the committed surface differs from the render",
    )
    args = parser.parse_args(argv)

    register = load_json(REGISTER_PATH)
    rows = validate(register)
    surface = render(register, rows).encode("utf-8")
    if args.check:
        committed = SURFACE_PATH.read_bytes() if SURFACE_PATH.is_file() else b""
        if committed != surface:
            print(
                "instrument register: docs/INSTRUMENT_REGISTER_V3.md is stale;"
                " run python tools/build_instrument_register.py",
                file=sys.stderr,
            )
            return 1
        print("instrument register: surface is current")
        return 0
    SURFACE_PATH.write_bytes(surface)
    print(f"instrument register: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
