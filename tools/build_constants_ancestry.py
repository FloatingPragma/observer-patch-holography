"""Build and validate the V3 constants ancestry surface (issue #736).

The machine-readable register is ``tracking/constants_ancestry.json``. This
tool validates it fail-closed and renders ``docs/CONSTANTS_ANCESTRY_V3.md``;
``--check`` fails when the committed page differs byte-for-byte from the
render.

One row per constants family in the masses-and-constants lane. Fail-closed
rules: the row ids, constants, statuses, and owning observation-ledger rows
are fixed program-wide and must equal the expected inventory exactly, in
order; every ancestry class comes from the six-value enum and every status
from the four-value enum; P and N entries cite the register row PR-17;
register-import entries cite an existing premise-register row or the literal
``proposed`` for a declared premise named descriptively; the register-row ids
cited on a row equal the premise list of its owning observation-ledger row;
rows with status ``postdiction_comparison`` carry a measured comparison
target in plain sight, ``diagnostic`` rows carry a measured comparison target
or a flagged empirical import, ``owed`` rows consume nothing, and
``conditional_determination`` rows carry no measured comparison target; every
evidence, value, and falsifier path exists in the repository. The rendered
page is a generated surface; edit the JSON, then regenerate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTER_PATH = ROOT / "tracking" / "constants_ancestry.json"
SURFACE_PATH = ROOT / "docs" / "CONSTANTS_ANCESTRY_V3.md"
PREMISE_REGISTER_PATH = ROOT / "tracking" / "premise_register.json"
OBSERVATION_LEDGER_PATH = ROOT / "tracking" / "observation_ledger.json"

SCHEMA = "oph.constants_ancestry.v1"
ISSUE = 736
ISSUE_URL = "https://github.com/muellerberndt/reverse-engineering-reality/issues"

CLASSES = (
    "P",
    "N",
    "register_import",
    "derived",
    "external_mathematics",
    "measured_comparison_target",
)
STATUSES = (
    "conditional_determination",
    "diagnostic",
    "postdiction_comparison",
    "owed",
)
PN_REGISTER_ROW = "PR-17"
PROPOSED = "proposed"

ROW_KEYS = {
    "id",
    "constant",
    "value_or_enclosure",
    "ancestry",
    "status",
    "ledger_row",
    "evidence",
    "falsifier",
}
POINTER_KEYS = {"artifact", "reading"}
FALSIFIER_KEYS = {"artifact", "rule"}

# Fixed program-wide inventory: (id, constant, status, ledger_row), in order.
EXPECTED_ROWS = (
    ("CA-01", "fine-structure constant", "diagnostic", "OL-H1"),
    ("CA-02", "cosmological constant via N", "diagnostic", "OL-H2"),
    ("CA-03", "charged-lepton masses", "postdiction_comparison", "OL-H3"),
    ("CA-04", "quark masses and mixings", "postdiction_comparison", "OL-H4"),
    ("CA-05", "W, Z, and Higgs sector", "postdiction_comparison", "OL-H5"),
    ("CA-06", "hadron spectrum", "diagnostic", "OL-H6"),
    ("CA-07", "neutrino stance", "owed", "OL-H7"),
    (
        "CA-08",
        "units and calibration: clock and energy anchors",
        "conditional_determination",
        "OL-H8",
    ),
)

STATUS_MEANING = {
    "conditional_determination": (
        "An exact statement certified conditional on the register rows and "
        "declared premises named on its row. The conditions are visible on "
        "the row, and the row presents no source-derived endpoint."
    ),
    "diagnostic": (
        "A certified numeric surface recorded for separation or transport "
        "against a measured value without a registered same-quantity "
        "bridge. The recorded numbers identify nothing and promote nothing."
    ),
    "postdiction_comparison": (
        "A row whose inputs include measured values as anchors or "
        "witnesses. Results on it are postdictions under the recorded "
        "decision rules, and enclosures and decision rules never move after "
        "comparison."
    ),
    "owed": (
        "No value, enclosure, or comparison is on record. The row consumes "
        "nothing and names the registration that must land first."
    ),
}


def fail(message: str) -> None:
    raise SystemExit(f"constants ancestry: {message}")


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        fail(f"missing input {path.relative_to(ROOT)}")
    except json.JSONDecodeError as error:
        fail(f"invalid JSON in {path.relative_to(ROOT)}: {error}")
    raise AssertionError("unreachable")


def load_premise_types() -> dict[str, str]:
    register = load_json(PREMISE_REGISTER_PATH)
    rows = register.get("rows")
    if not isinstance(rows, list) or not rows:
        fail("premise register rows must be a nonempty list")
    return {row["id"]: row["type"] for row in rows}


def load_ledger_premises() -> dict[str, list[str]]:
    ledger = load_json(OBSERVATION_LEDGER_PATH)
    rows = ledger.get("rows")
    if not isinstance(rows, list) or not rows:
        fail("observation ledger rows must be a nonempty list")
    return {row["id"]: list(row.get("premises", [])) for row in rows}


def validate_pointer(where: str, label: str, pointer: object, keys: set[str]) -> None:
    if not isinstance(pointer, dict):
        fail(f"{where}: {label} must be an object")
    if set(pointer) != keys:
        fail(f"{where}: {label} keys must be exactly {sorted(keys)}")
    for key in keys:
        value = pointer[key]
        if not isinstance(value, str) or not value.strip():
            fail(f"{where}: {label}.{key} must be a nonempty string")
    artifact = pointer["artifact"]
    if not (ROOT / artifact).exists():
        fail(f"{where}: {label} artifact path missing: {artifact}")


def validate_ancestry(
    where: str,
    status: str,
    ancestry: object,
    premise_types: dict[str, str],
) -> set[str]:
    if not isinstance(ancestry, list):
        fail(f"{where}: ancestry must be a list")
    if status == "owed":
        if ancestry:
            fail(f"{where}: an owed row consumes nothing; ancestry must be empty")
        return set()
    if not ancestry:
        fail(f"{where}: ancestry must be a nonempty list")

    cited: set[str] = set()
    classes_seen: list[str] = []
    inputs_seen: set[str] = set()
    empirical_import_seen = False
    for entry in ancestry:
        if not isinstance(entry, dict):
            fail(f"{where}: ancestry entries must be objects")
        keys = set(entry)
        if not {"input", "class"} <= keys or keys - {"input", "class", "register_row"}:
            fail(
                f"{where}: ancestry entry keys must be input, class, and "
                "optionally register_row"
            )
        text = entry["input"]
        if not isinstance(text, str) or not text.strip():
            fail(f"{where}: ancestry input must be a nonempty string")
        if text in inputs_seen:
            fail(f"{where}: ancestry inputs must be duplicate-free")
        inputs_seen.add(text)
        klass = entry["class"]
        if klass not in CLASSES:
            fail(f"{where}: class {klass!r} is not in the enum {CLASSES}")
        classes_seen.append(klass)

        if klass in ("P", "N"):
            if entry.get("register_row") != PN_REGISTER_ROW:
                fail(
                    f"{where}: a {klass} entry must cite register row "
                    f"{PN_REGISTER_ROW}"
                )
            cited.add(PN_REGISTER_ROW)
        elif klass == "register_import":
            if "register_row" not in entry:
                fail(f"{where}: a register_import entry must carry register_row")
            register_row = entry["register_row"]
            if register_row == PROPOSED:
                if "PR-" in text:
                    fail(
                        f"{where}: a proposed entry names its premise "
                        "descriptively and cites no PR id"
                    )
            elif register_row in premise_types:
                if register_row == PN_REGISTER_ROW:
                    fail(
                        f"{where}: {PN_REGISTER_ROW} enters through the P and "
                        "N classes, never as a register_import"
                    )
                cited.add(register_row)
                if premise_types[register_row] == "empirical_import":
                    empirical_import_seen = True
            else:
                fail(
                    f"{where}: register_row {register_row!r} is neither an "
                    f"existing premise-register row nor {PROPOSED!r}"
                )
        else:
            if "register_row" in entry:
                fail(
                    f"{where}: register_row is only legal on P, N, and "
                    "register_import entries"
                )

    measured = classes_seen.count("measured_comparison_target")
    if status == "postdiction_comparison" and measured == 0:
        fail(
            f"{where}: a postdiction_comparison row must carry a "
            "measured_comparison_target entry in plain sight"
        )
    if status == "diagnostic" and measured == 0 and not empirical_import_seen:
        fail(
            f"{where}: a diagnostic row must carry a measured comparison "
            "target or a flagged empirical import"
        )
    if status == "conditional_determination" and measured > 0:
        fail(
            f"{where}: a conditional_determination row carries no measured "
            "comparison target; comparisons live on the postdiction ledger"
        )
    return cited


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

    premise_types = load_premise_types()
    ledger_premises = load_ledger_premises()

    pn_classes_seen: set[str] = set()
    for index, row in enumerate(rows):
        expected_id, expected_constant, expected_status, expected_ledger = (
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
        if row["constant"] != expected_constant:
            fail(f"{where}: constant must equal {expected_constant!r}")
        if row["status"] not in STATUSES:
            fail(
                f"{where}: status {row['status']!r} is not in the enum "
                f"{STATUSES}"
            )
        if row["status"] != expected_status:
            fail(f"{where}: status must equal {expected_status!r}")
        if row["ledger_row"] != expected_ledger:
            fail(f"{where}: ledger_row must equal {expected_ledger!r}")
        if expected_ledger not in ledger_premises:
            fail(
                f"{where}: ledger row {expected_ledger} is absent from the "
                "observation ledger"
            )

        validate_pointer(where, "value_or_enclosure", row["value_or_enclosure"], POINTER_KEYS)
        validate_pointer(where, "falsifier", row["falsifier"], FALSIFIER_KEYS)

        cited = validate_ancestry(where, row["status"], row["ancestry"], premise_types)
        expected_premises = set(ledger_premises[expected_ledger])
        if cited != expected_premises:
            fail(
                f"{where}: cited register rows {sorted(cited)} must equal the "
                f"observation-ledger premises {sorted(expected_premises)}"
            )
        for entry in row["ancestry"]:
            pn_classes_seen.update(
                klass for klass in ("P", "N") if entry["class"] == klass
            )

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

    if pn_classes_seen != {"P", "N"}:
        fail("the surface must cite both numerical inputs P and N")
    return rows


def measured_display(row: dict, premise_types: dict[str, str]) -> str:
    flags: list[str] = []
    if any(
        entry["class"] == "measured_comparison_target" for entry in row["ancestry"]
    ):
        flags.append("measured_comparison_target")
    imports = sorted(
        entry["register_row"]
        for entry in row["ancestry"]
        if entry["class"] == "register_import"
        and entry.get("register_row") not in (PROPOSED, None)
        and premise_types.get(entry["register_row"]) == "empirical_import"
    )
    if imports:
        flags.append("empirical imports " + ", ".join(imports))
    return "; ".join(flags) if flags else "none"


def ancestry_lines(row: dict) -> list[str]:
    lines: list[str] = []
    if not row["ancestry"]:
        lines.append(
            "- none: the row consumes nothing until its registration lands."
        )
        return lines
    for entry in row["ancestry"]:
        klass = entry["class"]
        if klass == "measured_comparison_target":
            label = f"**`{klass}`**"
        else:
            label = f"`{klass}`"
        register_row = entry.get("register_row")
        if register_row == PROPOSED:
            label += " (declared premise; register row proposed)"
        elif register_row is not None:
            label += f" ({register_row})"
        lines.append(f"- {label}: {entry['input']}.")
    return lines


def render(rows: list[dict]) -> str:
    premise_types = load_premise_types()
    lines: list[str] = []
    lines.append("# OPH V3 Constants Ancestry")
    lines.append("")
    lines.append(
        "Generated by `tools/build_constants_ancestry.py` from"
        " `tracking/constants_ancestry.json`; edit the JSON, then regenerate."
        f" [Issue #{ISSUE}]({ISSUE_URL}/{ISSUE}) owns this surface."
    )
    lines.append("")
    lines.append(
        "One row per constants family in the masses-and-constants lane. Each"
        " row declares its complete input ancestry against the premise"
        " register: P and N are the only numerical inputs (register row"
        " PR-17), and every other number entering a row is classified in"
        " place as a register import, derived arithmetic, external"
        " mathematics, or a measured comparison target. A row that consumes"
        " a measured value carries that consumption in plain sight, through"
        " a **`measured_comparison_target`** entry or a flagged"
        " empirical-import register row. Declared premises without a register"
        " row are marked as proposed and named descriptively. The postdiction"
        " ledger and the claim scoreboard are the comparison authorities;"
        " this surface classifies ancestry, records falsifiers, and promotes"
        " nothing."
    )
    lines.append("")
    lines.append("| Row | Constants family | Status | Ledger row | Measured input |")
    lines.append("| --- | --- | --- | --- | --- |")
    for row in rows:
        lines.append(
            f"| {row['id']} | {row['constant']} | `{row['status']}` |"
            f" {row['ledger_row']} | {measured_display(row, premise_types)} |"
        )
    lines.append("")

    lines.append("## Rows")
    lines.append("")
    for row in rows:
        lines.append(f"### {row['id']} {row['constant']}")
        lines.append("")
        lines.append(
            f"Status `{row['status']}`; owning observation-ledger row"
            f" {row['ledger_row']}."
        )
        lines.append("")
        pointer = row["value_or_enclosure"]
        lines.append(
            f"- Value or enclosure: `{pointer['artifact']}`."
            f" {pointer['reading']}"
        )
        lines.append("- Ancestry:")
        for line in ancestry_lines(row):
            lines.append(f"  {line}")
        evidence = ", ".join(f"`{path}`" for path in row["evidence"])
        lines.append(f"- Evidence: {evidence}.")
        falsifier = row["falsifier"]
        lines.append(
            f"- Falsifier (`{falsifier['artifact']}`): {falsifier['rule']}"
        )
        lines.append("")

    lines.append("## What the statuses mean")
    lines.append("")
    for status in STATUSES:
        lines.append(f"### {status}")
        lines.append("")
        lines.append(STATUS_MEANING[status])
        lines.append("")

    lines.append(
        "Status changes on this surface are recorded register events, never"
        " silent edits. The ancestry surface is part of the anti-cheating"
        " boundary of the V3 program."
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
                "constants ancestry: docs/CONSTANTS_ANCESTRY_V3.md is stale; "
                "run python tools/build_constants_ancestry.py",
                file=sys.stderr,
            )
            return 1
        print("constants ancestry: surface is current")
        return 0
    SURFACE_PATH.write_bytes(surface)
    print(f"constants ancestry: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
