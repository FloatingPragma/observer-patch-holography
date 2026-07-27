"""Build and validate the canonical selection ledger surface.

The machine-readable ledger is ``claims/selection_ledger.json``. This tool
validates it fail-closed against the claim registry, the committed open-issue
snapshot, the physical-identification selector menus, the Lean corpus, and the
paper anchors, then renders ``docs/SELECTION_LEDGER.md`` from it. The rendered
page is a generated surface: ``--check`` fails when the committed page differs
from the render, and the mandatory suite runs that check.

Boundary rule enforced here (issue #554): every open selection names owner
issues that are live gates of its canonical or secondary claims (or live
blocking issues of a physical identification carrying one of those claims), so
the ledger, the claim registry, and the receipt registry expose identical
dependency boundaries.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "claims" / "selection_ledger.json"
REGISTRY_PATH = ROOT / "claims" / "claim_registry.yaml"
IDENTIFICATION_PATH = ROOT / "claims" / "physical_identification_registry.json"
SNAPSHOT_PATH = ROOT / "tracking" / "open_issues" / "open_problem_ledger.json"
SURFACE_PATH = ROOT / "docs" / "SELECTION_LEDGER.md"

SCHEMA = "oph.selection_ledger.v1"
CLASSES = ("forced", "exposed_premise", "open")
CLASS_TITLES = {
    "forced": "Forced by consensus",
    "exposed_premise": "Exposed premises",
    "open": "Open: the work-wanted list",
}

ROW_KEYS = {
    "row",
    "selector_id",
    "label",
    "class",
    "conditional_on",
    "canonical_claim_id",
    "secondary_claim_ids",
    "owner_issues",
    "where",
    "menu",
    "compression_input",
    "receipts",
    "lean",
    "negative_controls",
    "paper_anchors",
}

REQUIRED_SELECTOR_IDS = {
    "twelve_unit_port_splitting",
    "inverse_port_pairing_six_axes",
    "icosahedral_screen_selector_port_frame",
    "echosahedral_carrier_lineage",
    "compact_gauge_refinement_receipt",
    "declared_matter_completion_one_higgs_witness",
    "z6_reserve_pricing_input",
    "br0_cell_product_reading",
    "publicness_policy",
    "scheduler_selection_fair_block",
    "quantitative_closure_map_declaration",
    "physical_port_current_realization",
    "z6_global_form_descent",
    "physical_matter_spin_lift",
    "screen_to_family_attachment",
    "capacity_closure_selector",
    "br1_reserve_semantics",
    "br2_reserve_attachment",
    "br3_readback_count_effect",
    "br4_subfederation_marking",
    "br5_symmetry_quotient",
    "br6_cap_read_family",
}


def fail(message: str) -> None:
    raise SystemExit(f"selection ledger: {message}")


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        fail(f"missing input {path.relative_to(ROOT)}")
    except json.JSONDecodeError as error:
        fail(f"invalid JSON in {path.relative_to(ROOT)}: {error}")
    raise AssertionError("unreachable")


def positive_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


def validate(ledger: dict) -> list[dict]:
    if ledger.get("schema") != SCHEMA:
        fail(f"schema must equal {SCHEMA}")
    if ledger.get("issue") != 554:
        fail("issue must equal 554")
    rows = ledger.get("rows")
    if not isinstance(rows, list) or not rows:
        fail("rows must be a nonempty list")

    registry = load_json(REGISTRY_PATH)
    claim_ids = {claim["claim_id"] for claim in registry["claims"]}
    claim_gates = {claim["claim_id"]: set(claim["gates"]) for claim in registry["claims"]}

    identification = load_json(IDENTIFICATION_PATH)
    identification_blockers: dict[str, set[int]] = {}
    for entry in identification["physical_identifications"]:
        for claim_id in entry["claim_ids"]:
            identification_blockers.setdefault(claim_id, set()).update(
                blocker["number"] for blocker in entry["blocking_issues"]
            )
    selector_menus = {
        row["selector_id"]: row for row in identification["selector_menus"]
    }

    snapshot = load_json(SNAPSHOT_PATH)
    open_issues = {row["number"] for row in snapshot["rows"]}

    seen_ids: set[str] = set()
    undeclared_owned: set[int] = set()
    for index, row in enumerate(rows):
        where = f"rows[{index}]"
        if not isinstance(row, dict):
            fail(f"{where}: row must be an object")
        if set(row) != ROW_KEYS:
            missing = ROW_KEYS - set(row)
            extra = set(row) - ROW_KEYS
            fail(f"{where}: keys mismatch (missing {sorted(missing)}, extra {sorted(extra)})")
        if row["row"] != index + 1:
            fail(f"{where}: row numbers must be contiguous from 1")
        selector_id = row["selector_id"]
        if not isinstance(selector_id, str) or not selector_id:
            fail(f"{where}: selector_id must be a nonempty string")
        if selector_id in seen_ids:
            fail(f"duplicate selector_id {selector_id}")
        seen_ids.add(selector_id)
        if row["class"] not in CLASSES:
            fail(f"{where}: class must be one of {CLASSES}")
        conditional_on = row["conditional_on"]
        if (
            not isinstance(conditional_on, list)
            or any(not isinstance(item, str) or not item for item in conditional_on)
            or len(conditional_on) != len(set(conditional_on))
        ):
            fail(
                f"{where}: conditional_on must be a duplicate-free list of "
                "nonempty selector ids"
            )
        if not isinstance(row["label"], str) or not row["label"].strip():
            fail(f"{where}: label must be nonempty")
        if not isinstance(row["where"], str) or not row["where"].strip():
            fail(f"{where}: where must be nonempty")

        for claim_id in [row["canonical_claim_id"], *row["secondary_claim_ids"]]:
            if claim_id not in claim_ids:
                fail(f"{where}: unknown claim id {claim_id}")

        owners = row["owner_issues"]
        if not isinstance(owners, list) or any(not positive_int(n) for n in owners):
            fail(f"{where}: owner_issues must be positive integers")
        if row["class"] == "open":
            if not owners:
                fail(f"{where}: an open selection requires owner issues")
            allowed: set[int] = set()
            for claim_id in [row["canonical_claim_id"], *row["secondary_claim_ids"]]:
                allowed |= claim_gates[claim_id]
                allowed |= identification_blockers.get(claim_id, set())
            for number in owners:
                if number not in open_issues:
                    fail(f"{where}: owner issue #{number} is not open in the snapshot")
                if number not in allowed:
                    fail(
                        f"{where}: owner issue #{number} is not a live gate of the "
                        "row's claims; the ledger and registry boundaries diverge"
                    )
        else:
            if owners:
                fail(f"{where}: only open selections carry owner issues")

        menu = row["menu"]
        if set(menu) != {"size", "alternatives", "enumeration_source"}:
            fail(f"{where}: menu keys mismatch")
        size = menu["size"]
        if size is not None and not positive_int(size):
            fail(f"{where}: menu.size must be null or a positive integer")
        if not isinstance(menu["alternatives"], str) or not menu["alternatives"].strip():
            fail(f"{where}: menu.alternatives must describe the menu state")
        if row["class"] == "forced" and size is None:
            fail(f"{where}: a forced selection must record the menu it beat")
        if row["class"] == "forced" and not row["lean"] and not row["receipts"]:
            fail(f"{where}: a forced selection requires a Lean or receipt citation")

        if row["compression_input"]:
            entry = selector_menus.get(selector_id)
            if entry is None:
                fail(f"{where}: compression input missing from selector_menus")
            if size is None:
                if entry["status"] != "undeclared" or entry["menu_size"] is not None:
                    fail(f"{where}: selector_menus row must be undeclared with null size")
                if sorted(entry["blocking_issues"]) != sorted(owners):
                    fail(
                        f"{where}: selector_menus blockers {entry['blocking_issues']} "
                        f"must equal the owner issues {owners}"
                    )
                undeclared_owned.update(owners)
            else:
                if entry["status"] != "declared" or entry["menu_size"] != size:
                    fail(
                        f"{where}: selector_menus row must be declared with size {size}"
                    )
                if entry["blocking_issues"]:
                    fail(f"{where}: a declared selector cannot retain blockers")
        else:
            if size is not None:
                fail(f"{where}: a numeric menu must be a compression input")
            if selector_id in selector_menus:
                fail(f"{where}: non-compression row must stay out of selector_menus")

        for path in row["receipts"]:
            if not (ROOT / path).is_file():
                fail(f"{where}: receipt path missing: {path}")
        for ref in row["lean"]:
            module = ROOT / ref["module"]
            if not module.is_file():
                fail(f"{where}: Lean module missing: {ref['module']}")
            declaration = ref.get("declaration")
            if declaration is not None and declaration not in module.read_text(
                encoding="utf-8"
            ):
                fail(
                    f"{where}: declaration {declaration} absent from {ref['module']}"
                )
        for anchor in row["paper_anchors"]:
            paper = ROOT / anchor["file"]
            if not paper.is_file():
                fail(f"{where}: paper missing: {anchor['file']}")
            text = paper.read_text(encoding="utf-8")
            if (
                f"\\label{{{anchor['anchor']}}}" not in text
                and anchor["anchor"] not in text
            ):
                fail(
                    f"{where}: anchor {anchor['anchor']} absent from {anchor['file']}"
                )

    missing_required = REQUIRED_SELECTOR_IDS - seen_ids
    unexpected = seen_ids - REQUIRED_SELECTOR_IDS
    if missing_required or unexpected:
        fail(
            "selector inventory mismatch "
            f"(missing {sorted(missing_required)}, unexpected {sorted(unexpected)})"
        )

    dependencies = {
        row["selector_id"]: list(row["conditional_on"]) for row in rows
    }
    for selector_id, prerequisites in dependencies.items():
        unknown = set(prerequisites) - seen_ids
        if unknown:
            fail(
                f"selector {selector_id}: conditional_on names unknown selectors "
                f"{sorted(unknown)}"
            )
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(selector_id: str) -> None:
        if selector_id in visiting:
            fail(f"conditional_on cycle detected at selector {selector_id}")
        if selector_id in visited:
            return
        visiting.add(selector_id)
        for prerequisite in dependencies[selector_id]:
            visit(prerequisite)
        visiting.remove(selector_id)
        visited.add(selector_id)

    for selector_id in sorted(seen_ids):
        visit(selector_id)

    ledger_menu_ids = {
        row["selector_id"] for row in rows if row["compression_input"]
    }
    stray = set(selector_menus) - ledger_menu_ids
    if stray:
        fail(f"selector_menus rows without a ledger row: {sorted(stray)}")

    bound = identification["compression_bound"]
    if undeclared_owned:
        if set(bound["blocking_issues"]) != undeclared_owned:
            fail(
                "compression_bound blockers must equal the owner issues of "
                f"undeclared menus: {sorted(undeclared_owned)}"
            )
    return rows


def render(rows: list[dict]) -> str:
    by_class = {cls: [row for row in rows if row["class"] == cls] for cls in CLASSES}
    counts = {cls: len(by_class[cls]) for cls in CLASSES}

    def menu_cell(row: dict) -> str:
        size = row["menu"]["size"]
        if size is None:
            return row["menu"]["alternatives"]
        return f"{size} ({row['menu']['alternatives']})"

    lines: list[str] = []
    lines.append("# OPH Selection Ledger")
    lines.append("")
    lines.append(
        "Generated by `tools/build_selection_ledger.py` from"
        " `claims/selection_ledger.json`; edit the JSON, then regenerate."
    )
    lines.append("")
    lines.append(
        "OPH builds public reality from finite observers that read patches of a shared"
    )
    lines.append(
        "screen, compare overlaps, and repair disagreement until every observer agrees."
    )
    lines.append(
        "On the declared carrier lineage (row 4, an exposed premise) that consensus"
    )
    lines.append(
        "requirement forces the unit split, pairing, and frame of rows 1-3,"
    )
    lines.append(
        "with machine-checked proofs, and every remaining"
    )
    lines.append(
        "structural choice is finite: a declared premise with a stated menu, or an open"
    )
    lines.append(
        "item with a named owner. There is no tunable layer anywhere in the list. This"
    )
    lines.append(
        "ledger records every discrete structural selection in the framework exactly"
    )
    lines.append(
        "once, classified as FORCED (derived from the axioms and source objects), EXPOSED"
    )
    lines.append(
        "PREMISE (chosen and declared in theorem antecedents), or OPEN (neither, with the"
    )
    lines.append(
        "owning GitHub issue). It shows exactly what is proved and exactly what is left:"
    )
    lines.append(
        "pick an open row and work on it. It is the canonical free-versus-forced surface"
    )
    lines.append("for")
    lines.append(
        "[issue #554](https://github.com/FloatingPragma/observer-patch-holography/issues/554)."
    )
    lines.append("")

    lines.append(f"## {CLASS_TITLES['forced']}")
    lines.append("")
    lines.append("| # | Selection | Where it is forced | Menu beaten |")
    lines.append("| --- | --- | --- | --- |")
    for row in by_class["forced"]:
        lines.append(
            f"| {row['row']} | {row['label']} | {row['where']} | {menu_cell(row)} |"
        )
    lines.append("")

    lines.append(f"## {CLASS_TITLES['exposed_premise']}")
    lines.append("")
    lines.append("| # | Selection | Where it is declared | Menu |")
    lines.append("| --- | --- | --- | --- |")
    for row in by_class["exposed_premise"]:
        lines.append(
            f"| {row['row']} | {row['label']} | {row['where']} | {menu_cell(row)} |"
        )
    lines.append("")

    lines.append(f"## {CLASS_TITLES['open']}")
    lines.append("")
    lines.append("| # | Selection | Owner and state | Menu state |")
    lines.append("| --- | --- | --- | --- |")
    for row in by_class["open"]:
        owners = ", ".join(
            f"[#{n}](https://github.com/FloatingPragma/observer-patch-holography/issues/{n})"
            for n in row["owner_issues"]
        )
        lines.append(
            f"| {row['row']} | {row['label']} | {row['where']} Owner: {owners}. |"
            f" {menu_cell(row)} |"
        )
    lines.append("")

    lines.append(
        f"Totals: {len(rows)} selections. {counts['forced']} FORCED,"
        f" {counts['exposed_premise']} EXPOSED PREMISE, {counts['open']} OPEN."
    )
    lines.append(
        "Every open row belongs to the particle-realization and capacity programs;"
    )
    lines.append(
        "the recovery of relativity and quantum structure from observer consensus"
    )
    lines.append("consumes only the forced rows and declared premises above.")
    lines.append("")
    lines.append("## How to read a row")
    lines.append("")
    lines.append(
        "A FORCED row is a theorem: given the axioms and the cited premises, the"
    )
    lines.append(
        "selection has no alternative, the citation names the machine-checked or"
    )
    lines.append(
        "certified proof, and the menu column records the enumerated alternative"
    )
    lines.append(
        "set the theorem excluded. An EXPOSED PREMISE row is a declared choice: the"
    )
    lines.append(
        "theorem downstream of it states the choice in its antecedent, and the menu"
    )
    lines.append(
        "column records how many alternatives the declaration passed over. An OPEN"
    )
    lines.append(
        "row is a selection that is neither derived nor stated as a clean"
    )
    lines.append(
        "antecedent; the owning issue tracks the work that moves it into one of the"
    )
    lines.append(
        "other two classes. A row can move up (OPEN to EXPOSED PREMISE to FORCED)"
    )
    lines.append(
        "and never silently disappears: the consensus argument is exactly as strong"
    )
    lines.append(
        "as the FORCED column plus the declared premises it names. Rows with a"
    )
    lines.append(
        "numeric menu feed the compression accounting in"
    )
    lines.append(
        "`claims/physical_identification_registry.json`; rows whose declared object"
    )
    lines.append(
        "has no finite alternative enumeration carry their antecedent location"
    )
    lines.append(
        "instead and are excluded from compression accounting by construction."
    )
    lines.append(
        "Every open row is a well-posed problem with an issue link; solving one"
    )
    lines.append("moves it up the ledger.")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when the committed surface differs from the render",
    )
    args = parser.parse_args()

    ledger = load_json(LEDGER_PATH)
    rows = validate(ledger)
    surface = render(rows)
    if args.check:
        committed = (
            SURFACE_PATH.read_text(encoding="utf-8")
            if SURFACE_PATH.is_file()
            else ""
        )
        if committed != surface:
            print(
                "selection ledger: docs/SELECTION_LEDGER.md is stale; run "
                "python tools/build_selection_ledger.py",
                file=sys.stderr,
            )
            return 1
        print("selection ledger: surface is current")
        return 0
    SURFACE_PATH.write_text(surface, encoding="utf-8", newline="\n")
    print(f"selection ledger: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
