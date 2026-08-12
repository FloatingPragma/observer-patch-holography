"""Build and validate the V3 premise register surface (issue #727).

The machine-readable register is ``tracking/premise_register.json``. This tool
validates it fail-closed and renders ``docs/PREMISE_REGISTER_V3.md``;
``--check`` fails when the committed page differs byte-for-byte from the
render.

Fail-closed rules: the row ids, names, types, and dispositions are fixed
program-wide and must equal the expected inventory exactly, in order; every
type comes from the six-value enum and every disposition from the three-value
enum; every consuming lane is a V3 lane issue in 728..738, duplicate-free and
ascending; every evidence path exists in the repository; statements and notes
are nonempty. The rendered page is a generated surface; edit the JSON, then
regenerate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTER_PATH = ROOT / "tracking" / "premise_register.json"
SURFACE_PATH = ROOT / "docs" / "PREMISE_REGISTER_V3.md"

SCHEMA = "oph.premise_register.v1"
ISSUE = 727
ISSUE_URL = "https://github.com/muellerberndt/reverse-engineering-reality/issues"

TYPES = (
    "structural_rule",
    "representation_choice",
    "selection_rule",
    "external_mathematics",
    "empirical_import",
    "numerical_input",
)
DISPOSITIONS = ("remove", "axiomatize", "import")
LANE_MIN = 728
LANE_MAX = 738

ROW_KEYS = {
    "id",
    "name",
    "type",
    "disposition",
    "statement",
    "consuming_lanes",
    "evidence",
    "notes",
}

# Fixed program-wide inventory: (id, name, type, disposition), in order.
EXPECTED_ROWS = (
    ("PR-01", "confluent terminating repair contract", "structural_rule", "axiomatize"),
    (
        "PR-02",
        "algebra-state representation of public records",
        "representation_choice",
        "axiomatize",
    ),
    ("PR-03", "operational effect additivity", "selection_rule", "remove"),
    ("PR-04", "phase-operation instrument", "structural_rule", "remove"),
    ("PR-05", "counting path reference", "selection_rule", "remove"),
    ("PR-06", "real Legendre enrichment", "representation_choice", "remove"),
    (
        "PR-07",
        "declared repair law: common faithful reference and repaired-visible fibre",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-08",
        "refinement-uniform low-temperature control",
        "structural_rule",
        "remove",
    ),
    (
        "PR-09",
        "equal-face-weight and barycentric one-third port-dual rules",
        "selection_rule",
        "axiomatize",
    ),
    (
        "PR-10",
        "measure-to-metric and nearest-point repair rules",
        "selection_rule",
        "axiomatize",
    ),
    ("PR-11", "compact-Lie classification inputs", "external_mathematics", "import"),
    ("PR-12", "matter candidate grammar", "selection_rule", "remove"),
    ("PR-13", "Koide balance premise", "selection_rule", "remove"),
    ("PR-14", "hadronic transport packet", "empirical_import", "import"),
    ("PR-15", "clock and energy calibration anchors", "empirical_import", "import"),
    (
        "PR-16",
        "stable causality and open-image conditions",
        "structural_rule",
        "remove",
    ),
    ("PR-17", "numerical inputs P and N", "numerical_input", "import"),
    (
        "PR-18",
        "spherically symmetric radial readout",
        "structural_rule",
        "remove",
    ),
    ("PR-19", "shell-flux normalization", "structural_rule", "remove"),
    (
        "PR-20",
        "equal source-counting measure on the sixty directed source seams",
        "selection_rule",
        "axiomatize",
    ),
    (
        "PR-21",
        "auxiliary oscillator lift of the seam-current generator",
        "representation_choice",
        "axiomatize",
    ),
    (
        "PR-22",
        "physical-frequency identification for the seam-current symbol",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-23",
        "typed modular flow supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-24",
        "half-sided inclusion and null translation supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-25",
        "null-stress tomography supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-26",
        "generalized-entropy stationarity supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-27",
        "small-ball geometry supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-28",
        "physical tower and calibration supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-29",
        "layered shell-cardinality law",
        "structural_rule",
        "remove",
    ),
    (
        "PR-30",
        "layered shell equidistribution",
        "structural_rule",
        "remove",
    ),
    (
        "PR-31",
        "layered steady sourcing",
        "structural_rule",
        "remove",
    ),
    (
        "PR-32",
        "Cartesian joint carrier of an observer pair",
        "representation_choice",
        "axiomatize",
    ),
    (
        "PR-33",
        "two-slot diamond region map",
        "selection_rule",
        "axiomatize",
    ),
    (
        "PR-34",
        "two-observer slot split",
        "representation_choice",
        "axiomatize",
    ),
    (
        "PR-35",
        "same-source loop-to-kernel identity",
        "structural_rule",
        "remove",
    ),
    (
        "PR-36",
        "family band-selection rules: single complete faithful band and operational cost order",
        "selection_rule",
        "remove",
    ),
    (
        "PR-37",
        "fine-structure closure map selection",
        "selection_rule",
        "remove",
    ),
    (
        "PR-38",
        "same-quantity bridge for the P closure",
        "structural_rule",
        "remove",
    ),
    (
        "PR-39",
        "physical Thomson readback contract",
        "structural_rule",
        "remove",
    ),
    (
        "PR-40",
        "payload-coherent anchor-gap premise",
        "selection_rule",
        "remove",
    ),
    (
        "PR-41",
        "electroweak sector selection premises",
        "selection_rule",
        "remove",
    ),
    (
        "PR-42",
        "horizon record central-charge identification",
        "structural_rule",
        "axiomatize",
    ),
)

DISPOSITION_MEANING = {
    "remove": (
        "A remove row is consumed as a declared premise while its derivation "
        "from the three axioms is open work in the discharge queue (issue "
        "#739). The row leaves the register only through a proved source "
        "derivation; every composition citing it is conditional on it, and "
        "any bounded no-go attached to the row records exactly what the "
        "current axioms fail to supply."
    ),
    "axiomatize": (
        "An axiomatize row is a candidate for promotion into the "
        "architecture. Promotion is a recorded register decision with "
        "rationale and a basis-wide dependency and countermodel audit, never "
        "a silent edit; the row keeps its source-selection no-gos as the "
        "evidence of what the current axioms do not supply, and a promoted "
        "row remains on the register with its decision record."
    ),
    "import": (
        "An import row is a flagged external input that stays on the "
        "register: classical external mathematics, pinned empirical data "
        "packets, and the two numerical inputs P and N. An import is never "
        "presented as a derivation; replacing an import with a source-side "
        "result is a recorded disposition change."
    ),
}


def fail(message: str) -> None:
    raise SystemExit(f"premise register: {message}")


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        fail(f"missing input {path.relative_to(ROOT)}")
    except json.JSONDecodeError as error:
        fail(f"invalid JSON in {path.relative_to(ROOT)}: {error}")
    raise AssertionError("unreachable")


def validate(register: dict) -> list[dict]:
    if not isinstance(register, dict):
        fail("register must be an object")
    if register.get("schema") != SCHEMA:
        fail(f"schema must equal {SCHEMA}")
    if register.get("issue") != ISSUE:
        fail(f"issue must equal {ISSUE}")
    if set(register) != {"schema", "issue", "rows"}:
        fail("top-level keys must be exactly schema, issue, rows")
    rows = register.get("rows")
    if not isinstance(rows, list):
        fail("rows must be a list")
    if len(rows) != len(EXPECTED_ROWS):
        fail(f"rows must contain exactly {len(EXPECTED_ROWS)} entries")

    for index, row in enumerate(rows):
        expected_id, expected_name, expected_type, expected_disposition = (
            EXPECTED_ROWS[index]
        )
        where = f"row {expected_id}"
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
            fail(f"{where}: id must equal {expected_id}, found {row['id']!r}")
        if row["name"] != expected_name:
            fail(f"{where}: name must equal {expected_name!r}")
        if row["type"] not in TYPES:
            fail(f"{where}: type {row['type']!r} is not in the enum {TYPES}")
        if row["type"] != expected_type:
            fail(f"{where}: type must equal {expected_type!r}")
        if row["disposition"] not in DISPOSITIONS:
            fail(
                f"{where}: disposition {row['disposition']!r} is not in the "
                f"enum {DISPOSITIONS}"
            )
        if row["disposition"] != expected_disposition:
            fail(f"{where}: disposition must equal {expected_disposition!r}")

        statement = row["statement"]
        if not isinstance(statement, str) or not statement.strip():
            fail(f"{where}: statement must be a nonempty string")
        notes = row["notes"]
        if not isinstance(notes, str) or not notes.strip():
            fail(f"{where}: notes must be a nonempty string")

        lanes = row["consuming_lanes"]
        if not isinstance(lanes, list) or not lanes:
            fail(f"{where}: consuming_lanes must be a nonempty list")
        for lane in lanes:
            if not isinstance(lane, int) or isinstance(lane, bool):
                fail(f"{where}: consuming lane {lane!r} must be an integer")
            if not LANE_MIN <= lane <= LANE_MAX:
                fail(
                    f"{where}: consuming lane {lane} is outside "
                    f"{LANE_MIN}..{LANE_MAX}"
                )
        if lanes != sorted(set(lanes)):
            fail(f"{where}: consuming_lanes must be strictly ascending")

        evidence = row["evidence"]
        if not isinstance(evidence, list) or not evidence:
            fail(f"{where}: evidence must be a nonempty list")
        if len(evidence) != len(set(evidence)):
            fail(f"{where}: evidence paths must be duplicate-free")
        for path in evidence:
            if not isinstance(path, str) or not path.strip():
                fail(f"{where}: evidence entries must be nonempty strings")
            if not (ROOT / path).exists():
                fail(f"{where}: evidence path missing: {path}")
    return rows


def lane_links(lanes: list[int]) -> str:
    return ", ".join(f"[#{n}]({ISSUE_URL}/{n})" for n in lanes)


def render(rows: list[dict]) -> str:
    counts = {
        disposition: sum(1 for row in rows if row["disposition"] == disposition)
        for disposition in DISPOSITIONS
    }

    lines: list[str] = []
    lines.append("# OPH V3 Premise Register")
    lines.append("")
    lines.append(
        "Generated by `tools/build_premise_register.py` from"
        " `tracking/premise_register.json`; edit the JSON, then regenerate."
        f" [Issue #{ISSUE}]({ISSUE_URL}/{ISSUE}) owns this register."
    )
    lines.append("")
    lines.append(
        "Every named premise the V3 composition lanes consume appears here"
        " exactly once, with its type, its declared disposition, the lanes"
        " that consume it, its exact statement, and its evidence paths. A"
        " premise used by a lane theorem and absent from this register is a"
        " hard audit failure. P and N are the only numerical inputs; every"
        " other number entering a composition is derived, external"
        " mathematics, or a flagged empirical import on this register."
        " Registering a premise legitimizes its use inside conditional"
        " compositions; it does not make the premise true, derived, or"
        " physical. Disposition changes, including promotion to axiom, are"
        " recorded events with rationale, never silent edits."
    )
    lines.append("")
    lines.append("| Row | Premise | Type | Disposition | Consuming lanes |")
    lines.append("| --- | --- | --- | --- | --- |")
    for row in rows:
        lines.append(
            f"| {row['id']} | {row['name']} | `{row['type']}` |"
            f" `{row['disposition']}` | {lane_links(row['consuming_lanes'])} |"
        )
    lines.append("")
    lines.append(
        f"Totals: {len(rows)} premises. {counts['remove']} remove,"
        f" {counts['axiomatize']} axiomatize, {counts['import']} import."
    )
    lines.append("")

    lines.append("## Row statements and evidence")
    lines.append("")
    for row in rows:
        lines.append(f"### {row['id']} {row['name']}")
        lines.append("")
        lines.append(row["statement"])
        lines.append("")
        lines.append(
            f"- Type `{row['type']}`; disposition `{row['disposition']}`;"
            f" consumed by {lane_links(row['consuming_lanes'])}."
        )
        evidence = ", ".join(f"`{path}`" for path in row["evidence"])
        lines.append(f"- Evidence: {evidence}.")
        lines.append(f"- Disposition note: {row['notes']}")
        lines.append("")

    lines.append("## What the dispositions mean")
    lines.append("")
    for disposition in DISPOSITIONS:
        lines.append(f"### {disposition}")
        lines.append("")
        lines.append(DISPOSITION_MEANING[disposition])
        lines.append("")

    lines.append(
        "Register rows keep the source-selection no-gos proved against them"
        " as permanent evidence of what the current axioms do not supply."
        " The register is the anti-cheating surface of the V3 program."
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when the committed surface differs from the render",
    )
    args = parser.parse_args()

    register = load_json(REGISTER_PATH)
    rows = validate(register)
    surface = render(rows).encode("utf-8")
    if args.check:
        committed = SURFACE_PATH.read_bytes() if SURFACE_PATH.is_file() else b""
        if committed != surface:
            print(
                "premise register: docs/PREMISE_REGISTER_V3.md is stale; run "
                "python tools/build_premise_register.py",
                file=sys.stderr,
            )
            return 1
        print("premise register: surface is current")
        return 0
    SURFACE_PATH.write_bytes(surface)
    print(f"premise register: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
