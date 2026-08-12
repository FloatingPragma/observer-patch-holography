"""Build and validate the Standard Model Lagrangian correspondence table.

Issue #735's headline deliverable: one row per term of the textbook Standard
Model Lagrangian, each classified by what OPH supplies for it. The
machine-readable table is ``tracking/sm_lagrangian_correspondence.json``.
This tool validates it fail-closed and renders
``docs/SM_LAGRANGIAN_CORRESPONDENCE.md``; ``--check`` fails when the
committed page differs from the render.

Fail-closed rules: row ids are contiguous in the SML-01 style; every
classification is one of the four enum values; a derived_conditional or
registered_premise row names at least one premise id and an absent row names
none; every premise id (in the premises list or cited inline in a boundary)
is a row of the fixed program-wide premise register, which the artifact
embeds verbatim; every evidence path resolves to a committed file; a
derived_conditional or partial row cites at least one evidence path; and the
rendered prose carries no banned wording.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLE_PATH = ROOT / "tracking" / "sm_lagrangian_correspondence.json"
SURFACE_PATH = ROOT / "docs" / "SM_LAGRANGIAN_CORRESPONDENCE.md"

SCHEMA = "oph.sm_lagrangian_correspondence.v1"
ISSUE = 735

CLASSIFICATIONS = (
    "derived_conditional",
    "partial",
    "registered_premise",
    "absent",
)
CLASS_LABELS = {
    "derived_conditional": "derived conditional",
    "partial": "partial",
    "registered_premise": "registered premise",
    "absent": "absent",
}

TOP_KEYS = {"schema", "issue", "policy", "premise_register", "rows"}
ROW_KEYS = {"id", "term", "classification", "premises", "evidence", "boundary"}

# The program-wide premise register rows (issue #727). Ids, names, types, and
# dispositions are fixed program-wide; the artifact embeds them verbatim and
# validation rejects any divergence.
PREMISE_REGISTER: tuple[dict, ...] = (
    {
        "id": "PR-01",
        "name": "confluent terminating repair contract",
        "type": "structural_rule",
        "disposition": "axiomatize",
    },
    {
        "id": "PR-02",
        "name": "algebra-state representation of public records",
        "type": "representation_choice",
        "disposition": "axiomatize",
    },
    {
        "id": "PR-03",
        "name": "operational effect additivity",
        "type": "selection_rule",
        "disposition": "remove",
    },
    {
        "id": "PR-04",
        "name": "phase-operation instrument",
        "type": "structural_rule",
        "disposition": "remove",
    },
    {
        "id": "PR-05",
        "name": "counting path reference",
        "type": "selection_rule",
        "disposition": "remove",
    },
    {
        "id": "PR-06",
        "name": "real Legendre enrichment",
        "type": "representation_choice",
        "disposition": "remove",
    },
    {
        "id": "PR-07",
        "name": (
            "declared repair law: common faithful reference and "
            "repaired-visible fibre"
        ),
        "type": "structural_rule",
        "disposition": "axiomatize",
    },
    {
        "id": "PR-08",
        "name": "refinement-uniform low-temperature control",
        "type": "structural_rule",
        "disposition": "remove",
    },
    {
        "id": "PR-09",
        "name": "equal-face-weight and barycentric one-third port-dual rules",
        "type": "selection_rule",
        "disposition": "axiomatize",
    },
    {
        "id": "PR-10",
        "name": "measure-to-metric and nearest-point repair rules",
        "type": "selection_rule",
        "disposition": "axiomatize",
    },
    {
        "id": "PR-11",
        "name": "compact-Lie classification inputs",
        "type": "external_mathematics",
        "disposition": "import",
    },
    {
        "id": "PR-12",
        "name": "matter candidate grammar",
        "type": "selection_rule",
        "disposition": "remove",
    },
    {
        "id": "PR-13",
        "name": "Koide balance premise",
        "type": "selection_rule",
        "disposition": "remove",
    },
    {
        "id": "PR-14",
        "name": "hadronic transport packet",
        "type": "empirical_import",
        "disposition": "import",
    },
    {
        "id": "PR-15",
        "name": "clock and energy calibration anchors",
        "type": "empirical_import",
        "disposition": "import",
    },
    {
        "id": "PR-16",
        "name": "stable causality and open-image conditions",
        "type": "structural_rule",
        "disposition": "remove",
    },
    {
        "id": "PR-17",
        "name": "numerical inputs P and N",
        "type": "numerical_input",
        "disposition": "import",
    },
)

PR_TOKEN = re.compile(r"PR-\d{2}")
_H_WORD = "hon" + "est"
BANNED_WORDS = re.compile(
    r"\b("
    + "|".join(
        (
            _H_WORD,
            _H_WORD + "ly",
            _H_WORD + "y",
            "now",
            "currently",
            "previously",
            "already",
            "recently",
        )
    )
    + r")\b",
    re.IGNORECASE,
)
BANNED_FRAGMENTS = ("—",)


def fail(message: str) -> None:
    raise SystemExit(f"sm correspondence: {message}")


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        fail(f"missing input {path.relative_to(ROOT)}")
    except json.JSONDecodeError as error:
        fail(f"invalid JSON in {path.relative_to(ROOT)}: {error}")
    raise AssertionError("unreachable")


def check_prose(where: str, text: str) -> None:
    for fragment in BANNED_FRAGMENTS:
        if fragment in text:
            fail(f"{where}: banned fragment {fragment!r}")
    match = BANNED_WORDS.search(text)
    if match:
        fail(f"{where}: banned word {match.group(0)!r}")


def validate(data: dict) -> list[dict]:
    if not isinstance(data, dict) or set(data) != TOP_KEYS:
        fail(f"top-level keys must equal {sorted(TOP_KEYS)}")
    if data["schema"] != SCHEMA:
        fail(f"schema must equal {SCHEMA}")
    if data["issue"] != ISSUE:
        fail(f"issue must equal {ISSUE}")
    if not isinstance(data["policy"], str) or not data["policy"].strip():
        fail("policy must be a nonempty string")
    check_prose("policy", data["policy"])

    if data["premise_register"] != list(PREMISE_REGISTER):
        fail("premise_register must match the fixed program-wide rows verbatim")
    register_ids = {entry["id"] for entry in PREMISE_REGISTER}

    rows = data["rows"]
    if not isinstance(rows, list) or not rows:
        fail("rows must be a nonempty list")
    for index, row in enumerate(rows):
        expected_id = f"SML-{index + 1:02d}"
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
        if row["id"] != expected_id:
            fail(f"{where}: id must equal {expected_id}")
        where = row["id"]
        term = row["term"]
        if not isinstance(term, str) or ": " not in term:
            fail(
                f"{where}: term must be a string of the form "
                "'textbook name: schematic expression'"
            )
        check_prose(f"{where}.term", term)
        classification = row["classification"]
        if classification not in CLASSIFICATIONS:
            fail(
                f"{where}: classification {classification!r} is not one of "
                f"{CLASSIFICATIONS}"
            )
        premises = row["premises"]
        if (
            not isinstance(premises, list)
            or any(not isinstance(p, str) for p in premises)
            or len(premises) != len(set(premises))
        ):
            fail(f"{where}: premises must be a duplicate-free list of strings")
        unknown = set(premises) - register_ids
        if unknown:
            fail(f"{where}: unknown premise ids {sorted(unknown)}")
        if classification in ("derived_conditional", "registered_premise"):
            if not premises:
                fail(
                    f"{where}: a {classification} row must name at least one "
                    "premise register row"
                )
        if classification == "absent" and premises:
            fail(f"{where}: an absent row names no premises")
        evidence = row["evidence"]
        if (
            not isinstance(evidence, list)
            or any(not isinstance(p, str) for p in evidence)
            or len(evidence) != len(set(evidence))
        ):
            fail(f"{where}: evidence must be a duplicate-free list of paths")
        for path in evidence:
            if not (ROOT / path).is_file():
                fail(f"{where}: evidence path missing: {path}")
        if classification in ("derived_conditional", "partial") and not evidence:
            fail(
                f"{where}: a {classification} row must cite at least one "
                "evidence path"
            )
        boundary = row["boundary"]
        if (
            not isinstance(boundary, str)
            or not boundary.strip()
            or not boundary.endswith(".")
        ):
            fail(f"{where}: boundary must be a nonempty sentence ending in '.'")
        check_prose(f"{where}.boundary", boundary)
        inline_unknown = set(PR_TOKEN.findall(boundary)) - register_ids
        if inline_unknown:
            fail(
                f"{where}: boundary cites unknown premise ids "
                f"{sorted(inline_unknown)}"
            )
    return rows


def referenced_register_rows(rows: list[dict]) -> list[dict]:
    ids: set[str] = set()
    for row in rows:
        ids.update(row["premises"])
        ids.update(PR_TOKEN.findall(row["boundary"]))
    return [entry for entry in PREMISE_REGISTER if entry["id"] in ids]


def render(data: dict, rows: list[dict]) -> str:
    counts = {cls: 0 for cls in CLASSIFICATIONS}
    for row in rows:
        counts[row["classification"]] += 1

    lines: list[str] = []
    lines.append("# The Standard Model Lagrangian correspondence table")
    lines.append("")
    lines.append(
        "Generated by `tools/build_sm_correspondence.py` from"
        " `tracking/sm_lagrangian_correspondence.json`; edit the JSON, then"
        " regenerate. Issue #735 owns this surface."
    )
    lines.append("")
    lines.append(data["policy"])
    lines.append("")
    lines.append(
        f"Summary: {len(rows)} terms;"
        f" {counts['derived_conditional']} derived conditional,"
        f" {counts['partial']} partial,"
        f" {counts['registered_premise']} registered premise,"
        f" {counts['absent']} absent."
    )
    lines.append("")
    lines.append("| Row | Term | Classification | Premises |")
    lines.append("| --- | --- | --- | --- |")
    for row in rows:
        name, expression = row["term"].split(": ", 1)
        premises = ", ".join(row["premises"]) if row["premises"] else "none"
        lines.append(
            f"| {row['id']} | {name}: `{expression}` |"
            f" `{row['classification']}` | {premises} |"
        )
    lines.append("")
    lines.append("## Term boundaries")
    lines.append("")
    lines.append(
        "One paragraph per row: exactly what is supplied, exactly what is"
        " not, and the committed evidence."
    )
    for row in rows:
        name = row["term"].split(": ", 1)[0]
        evidence = (
            ", ".join(f"`{path}`" for path in row["evidence"])
            if row["evidence"]
            else "none"
        )
        lines.append("")
        lines.append(
            f"**{row['id']}, {name}"
            f" ({CLASS_LABELS[row['classification']]}).**"
            f" {row['boundary']} Evidence: {evidence}."
        )
    lines.append("")
    lines.append("## Premise register key")
    lines.append("")
    lines.append(
        "Premise ids name rows of the program-wide premise register"
        " (issue #727). Rows referenced by this table:"
    )
    lines.append("")
    for entry in referenced_register_rows(rows):
        lines.append(
            f"- `{entry['id']}` {entry['name']}"
            f" (type {entry['type']}, disposition {entry['disposition']})"
        )
    lines.append("")
    lines.append(
        "The classification enum is closed: `derived_conditional` is an"
        " exact theorem shape under the named register rows,"
        " `registered_premise` is consumed as a declared register row,"
        " `partial` carries the attained structure stated in its boundary,"
        " and `absent` has no OPH structure. A row upgrade requires"
        " receipts, a premise list, and a regeneration of this surface."
    )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when the committed surface differs from the render",
    )
    args = parser.parse_args(argv)

    data = load_json(TABLE_PATH)
    rows = validate(data)
    surface = render(data, rows)
    if args.check:
        committed = SURFACE_PATH.read_bytes() if SURFACE_PATH.is_file() else b""
        if committed != surface.encode("utf-8"):
            print(
                "sm correspondence: docs/SM_LAGRANGIAN_CORRESPONDENCE.md is"
                " stale; run python tools/build_sm_correspondence.py",
                file=sys.stderr,
            )
            return 1
        print("sm correspondence: surface is current")
        return 0
    SURFACE_PATH.write_bytes(surface.encode("utf-8"))
    print(f"sm correspondence: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
