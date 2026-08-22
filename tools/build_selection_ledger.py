"""Build and validate the canonical selection ledger surface.

The machine-readable ledger is ``claims/selection_ledger.json``. This tool
validates it fail-closed against the claim registry, the physical-identification
selector menus, the Lean corpus, and the paper anchors, then renders
``docs/SELECTION_LEDGER.md`` from it. The rendered
page is a generated surface: ``--check`` fails when the committed page differs
from the render, and the mandatory suite runs that check.

Boundary rule enforced here (issue #554): every open selection names owner
issues that are declared gates of its canonical or secondary claims (or
blocking issues of a physical identification carrying one of those claims), so
the scientific registries expose identical dependency boundaries. Live issue
state belongs to GitHub and the workspace DAG, not to this repository.
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
    "br0_cell_product_reading",
    "publicness_policy",
    "scheduler_selection_fair_block",
    "quantitative_closure_map_declaration",
    "physical_port_current_realization",
    "oriented_face_minimum_repair_rule",
    "post_hoc_repair_count_phase_pairing",
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


def require_exact_keys(value: object, expected: set[str], where: str) -> dict:
    if not isinstance(value, dict):
        fail(f"{where}: value must be an object")
    actual = set(value)
    if actual != expected:
        fail(
            f"{where}: keys mismatch "
            f"(missing {sorted(expected - actual)}, extra {sorted(actual - expected)})"
        )
    return value


def validate_physical_identification_registry(
    identification: object,
    claim_ids: set[str],
    root: Path = ROOT,
) -> tuple[dict[str, set[int]], dict[str, dict]]:
    """Validate scientific identification and compression bookkeeping.

    This deliberately validates declarations only. Whether an issue is open,
    closed, or superseded is live project state owned by GitHub and the
    workspace DAG, not by this repository.
    """

    payload = require_exact_keys(
        identification,
        {
            "schema_version",
            "physical_identifications",
            "selector_menus",
            "compression_bound",
        },
        "physical identification registry",
    )
    if payload["schema_version"] != 1:
        fail("physical identification registry: schema_version must equal 1")

    identifications = payload["physical_identifications"]
    if not isinstance(identifications, list) or not identifications:
        fail("physical identification registry: physical_identifications must be nonempty")

    seen_identifications: set[str] = set()
    identification_blockers: dict[str, set[int]] = {}
    for index, raw_row in enumerate(identifications):
        where = f"physical_identifications[{index}]"
        row = require_exact_keys(
            raw_row,
            {
                "identification_id",
                "label",
                "status",
                "claim_ids",
                "source_anchors",
                "blocking_issues",
            },
            where,
        )
        identification_id = row["identification_id"]
        if not isinstance(identification_id, str) or not identification_id:
            fail(f"{where}: identification_id must be a nonempty string")
        if identification_id in seen_identifications:
            fail(f"duplicate physical identification id {identification_id}")
        seen_identifications.add(identification_id)
        if not isinstance(row["label"], str) or not row["label"].strip():
            fail(f"{where}: label must be a nonempty string")
        if row["status"] != "undischarged":
            fail(f"{where}: status must equal undischarged")

        row_claim_ids = row["claim_ids"]
        if (
            not isinstance(row_claim_ids, list)
            or not row_claim_ids
            or any(not isinstance(claim_id, str) or not claim_id for claim_id in row_claim_ids)
            or len(set(row_claim_ids)) != len(row_claim_ids)
        ):
            fail(f"{where}: claim_ids must be a nonempty unique string list")
        unknown_claims = sorted(set(row_claim_ids) - claim_ids)
        if unknown_claims:
            fail(f"{where}: unknown canonical claim ids {unknown_claims}")

        anchors = row["source_anchors"]
        if (
            not isinstance(anchors, list)
            or not anchors
            or any(not isinstance(anchor, str) or not anchor.strip() for anchor in anchors)
        ):
            fail(f"{where}: source_anchors must be nonempty strings")

        blockers = row["blocking_issues"]
        if not isinstance(blockers, list) or not blockers:
            fail(f"{where}: blocking_issues must be nonempty")
        blocker_numbers: list[int] = []
        for blocker_index, raw_blocker in enumerate(blockers):
            blocker_where = f"{where}.blocking_issues[{blocker_index}]"
            blocker = require_exact_keys(
                raw_blocker,
                {"number", "scope", "role"},
                blocker_where,
            )
            number = blocker["number"]
            if not positive_int(number):
                fail(f"{blocker_where}: number must be a positive integer")
            if blocker["scope"] not in {
                "unresolved_attachment",
                "resource_deferred",
                "bounded_out_of_scope",
            }:
                fail(f"{blocker_where}: scope is outside the controlled vocabulary")
            if not isinstance(blocker["role"], str) or not blocker["role"].strip():
                fail(f"{blocker_where}: role must be a nonempty string")
            blocker_numbers.append(number)
        if len(set(blocker_numbers)) != len(blocker_numbers):
            fail(f"{where}: blocking issue numbers must be unique")
        for claim_id in row_claim_ids:
            identification_blockers.setdefault(claim_id, set()).update(blocker_numbers)

    raw_menus = payload["selector_menus"]
    if not isinstance(raw_menus, list) or not raw_menus:
        fail("physical identification registry: selector_menus must be nonempty")
    selector_menus: dict[str, dict] = {}
    undeclared_blockers: set[int] = set()
    for index, raw_menu in enumerate(raw_menus):
        where = f"selector_menus[{index}]"
        menu = require_exact_keys(
            raw_menu,
            {"selector_id", "label", "source", "status", "menu_size", "blocking_issues"},
            where,
        )
        selector_id = menu["selector_id"]
        if not isinstance(selector_id, str) or not selector_id:
            fail(f"{where}: selector_id must be a nonempty string")
        if selector_id in selector_menus:
            fail(f"duplicate physical selector id {selector_id}")
        selector_menus[selector_id] = menu
        if not isinstance(menu["label"], str) or not menu["label"].strip():
            fail(f"{where}: label must be a nonempty string")
        source = menu["source"]
        if not isinstance(source, str) or not source:
            fail(f"{where}: source must be a repository-relative path")
        relative_source = Path(source)
        if relative_source.is_absolute():
            fail(f"{where}: source must be a repository-relative path")
        try:
            source_path = (root / relative_source).resolve()
            source_path.relative_to(root.resolve())
        except ValueError:
            fail(f"{where}: source escapes the repository")
        if not source_path.is_file():
            fail(f"{where}: source does not exist: {source}")

        blockers = menu["blocking_issues"]
        if (
            not isinstance(blockers, list)
            or any(not positive_int(number) for number in blockers)
            or len(set(blockers)) != len(blockers)
        ):
            fail(f"{where}: blocking_issues must be unique positive integers")
        if menu["status"] == "declared":
            if not positive_int(menu["menu_size"]):
                fail(f"{where}: a declared selector requires a positive integer menu_size")
            if blockers:
                fail(f"{where}: a declared selector cannot retain blocking_issues")
        elif menu["status"] == "undeclared":
            if menu["menu_size"] is not None:
                fail(f"{where}: an undeclared selector cannot carry a numeric menu_size")
            if not blockers:
                fail(f"{where}: an undeclared selector requires blocking issues")
            undeclared_blockers.update(blockers)
        else:
            fail(f"{where}: status must equal declared or undeclared")

    bound = require_exact_keys(
        payload["compression_bound"],
        {"status", "reason", "blocking_issues", "p_acc_upper_bound"},
        "compression_bound",
    )
    bound_blockers = bound["blocking_issues"]
    if (
        not isinstance(bound_blockers, list)
        or any(not positive_int(number) for number in bound_blockers)
        or len(set(bound_blockers)) != len(bound_blockers)
    ):
        fail("compression_bound: blocking_issues must be unique positive integers")
    if undeclared_blockers:
        if bound["p_acc_upper_bound"] is not None:
            fail("compression_bound: numeric P_acc is forbidden while menu sizes are undeclared")
        if bound["status"] != "not_computable" or bound["reason"] != "menu_sizes_undeclared":
            fail(
                "compression_bound: undeclared menus require "
                "not_computable / menu_sizes_undeclared"
            )
        if set(bound_blockers) != undeclared_blockers:
            fail("compression_bound: blockers must exactly match undeclared selector menus")
    else:
        if bound_blockers:
            fail("compression_bound: declared selector menus cannot retain blockers")
        if bound["status"] == "not_computed":
            if bound["reason"] != "calculation_not_run" or bound["p_acc_upper_bound"] is not None:
                fail(
                    "compression_bound: not_computed requires calculation_not_run "
                    "and no numeric P_acc"
                )
        elif bound["status"] == "computed":
            value = bound["p_acc_upper_bound"]
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not 0 <= value <= 1
                or bound["reason"] != "all_menu_sizes_declared"
            ):
                fail(
                    "compression_bound: computed requires an all-menu numeric P_acc "
                    "between zero and one"
                )
        else:
            fail("compression_bound: fully declared menus require computed or not_computed")

    return identification_blockers, selector_menus


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
    identification_blockers, selector_menus = (
        validate_physical_identification_registry(identification, claim_ids)
    )

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
                if number not in allowed:
                    fail(
                        f"{where}: owner issue #{number} is not a declared gate of the "
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
        "requirement forces the pairing and frame of rows 2-3. The unit split in row 1"
    )
    lines.append(
        "is exact inside a separately declared integer counting and cost realization."
    )
    lines.append(
        "These statements have machine-checked proofs. Finite menus are enumerated"
    )
    lines.append(
        "where they exist. Other exposed premises include continuous data or"
    )
    lines.append(
        "unenumerated maps and are excluded from finite compression accounting."
    )
    lines.append(
        "Across the premise register and this ledger, declared inputs and unresolved"
    )
    lines.append(
        "attachments are named. This ledger is the bounded discrete-selection and"
    )
    lines.append(
        "compression inventory defined in `claims/selection_ledger.json`; it is not an"
    )
    lines.append(
        "exhaustive inventory of every theorem premise. Declared effects, evolution laws,"
    )
    lines.append(
        "scheduler clauses, and finite action packets belong in the premise register."
    )
    lines.append(
        "Rows here are classified as FORCED (derived from the axioms and source objects), EXPOSED"
    )
    lines.append(
        "PREMISE (chosen and declared in theorem antecedents), or OPEN (neither, with the"
    )
    lines.append(
        "owning GitHub issue). It shows what is proved and open inside this finite menu:"
    )
    lines.append(
        "pick an open row and work on it. It is the free-versus-forced surface"
    )
    lines.append("for")
    lines.append(
        "[historical issue #554](https://github.com/FloatingPragma/observer-patch-holography/issues/554);"
    )
    lines.append(
        "the premise register was established under [#727](https://github.com/FloatingPragma/observer-patch-holography/issues/727) and remains a maintained scientific register,"
    )
    lines.append(
        "and source discharge is [#739](https://github.com/FloatingPragma/observer-patch-holography/issues/739)."
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
        "Open rows belong to bracket, particle-realization, and capacity programs;"
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
