#!/usr/bin/env python3
"""Generate the single-page claims scoreboard (Phase 0 of the completion plan).

The scoreboard is one generated artifact: physical-identification blockers,
claims x class x status x live gates, and compression computability. Its
fail-closed sources are the claim registry, the physical-identification
registry, and the committed GitHub issue snapshot. It is never edited by hand;
``--check`` fails when the committed page drifts from its sources, and the
mandatory suite runs that check on every push.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

PROJECTION_START = "<!-- CLAIM-STATUS-PROJECTION:BEGIN -->"
PROJECTION_END = "<!-- CLAIM-STATUS-PROJECTION:END -->"

COMPRESSION_TEMPLATE = """# OPH Compression Scorecard

This scorecard deliberately does not collapse heterogeneous mathematical,
artifact, empirical, and physical-establishment rows into one scalar evidence
score. The canonical claim classes, descriptive statuses, and live gates are
the only fail-closed accounting on this surface. A quantitative compression
claim must name its input menu and target exposure separately.

<!-- CLAIM-STATUS-PROJECTION:BEGIN -->
<!-- Regenerate with: python tools/build_scoreboard.py -->
<!-- CLAIM-STATUS-PROJECTION:END -->

## Reading rule

Zero adjustable continuous parameters in a declared map is not, by itself,
physical establishment or statistical compression. The claim class and live
gates in the generated projection remain binding.
"""

CLASS_ORDER = [
    "physical_establishment",
    "empirical_implementation",
    "emitted_artifact",
    "branch_entry",
    "conditional_implication",
    "declared_structure",
]


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_documents(root: Path = ROOT) -> tuple[dict, dict, dict]:
    return (
        load(root / "claims" / "claim_registry.yaml"),
        load(root / "tracking" / "open_issues" / "open_problem_ledger.json"),
        load(root / "claims" / "physical_identification_registry.json"),
    )


def _exact_keys(row: dict, expected: set[str], where: str) -> None:
    actual = set(row)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise SystemExit(f"{where}: closed schema mismatch, missing={missing}, extra={extra}")


def _positive_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


def _issue_link(repo: str, number: int) -> str:
    return f"[#{number}](https://github.com/{repo}/issues/{number})"


def validate_accounting(
    registry: dict,
    snapshot: dict,
    accounting: dict,
    root: Path = ROOT,
) -> None:
    """Validate physical blockers and refuse premature compression numbers."""

    if not isinstance(accounting, dict):
        raise SystemExit("physical identification registry: root must be an object")
    _exact_keys(
        accounting,
        {
            "schema_version",
            "physical_identifications",
            "selector_menus",
            "compression_bound",
        },
        "physical identification registry",
    )
    if accounting["schema_version"] != 1:
        raise SystemExit("physical identification registry: schema_version must equal 1")

    claim_ids = {claim["claim_id"] for claim in registry.get("claims", [])}
    open_issues = {row["number"] for row in snapshot.get("rows", [])}
    closed_out_of_scope_issues = {
        row["number"] for row in snapshot.get("closed_out_of_scope_records", [])
    }

    identifications = accounting["physical_identifications"]
    if not isinstance(identifications, list) or not identifications:
        raise SystemExit(
            "physical identification registry: physical_identifications must be nonempty"
        )
    seen_identifications: set[str] = set()
    for index, row in enumerate(identifications):
        where = f"physical_identifications[{index}]"
        if not isinstance(row, dict):
            raise SystemExit(f"{where}: row must be an object")
        _exact_keys(
            row,
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
            raise SystemExit(f"{where}: identification_id must be a nonempty string")
        if identification_id in seen_identifications:
            raise SystemExit(
                f"physical identification registry: duplicate identification_id "
                f"{identification_id}"
            )
        seen_identifications.add(identification_id)
        if not isinstance(row["label"], str) or not row["label"].strip():
            raise SystemExit(f"{where}: label must be a nonempty string")
        if row["status"] != "undischarged":
            raise SystemExit(f"{where}: status must equal undischarged")

        row_claim_ids = row["claim_ids"]
        if (
            not isinstance(row_claim_ids, list)
            or not row_claim_ids
            or any(not isinstance(claim_id, str) for claim_id in row_claim_ids)
            or len(set(row_claim_ids)) != len(row_claim_ids)
        ):
            raise SystemExit(f"{where}: claim_ids must be a nonempty unique list")
        unknown_claims = sorted(set(row_claim_ids) - claim_ids)
        if unknown_claims:
            raise SystemExit(f"{where}: unknown canonical claim ids {unknown_claims}")

        anchors = row["source_anchors"]
        if (
            not isinstance(anchors, list)
            or not anchors
            or any(not isinstance(anchor, str) or not anchor.strip() for anchor in anchors)
        ):
            raise SystemExit(f"{where}: source_anchors must be nonempty strings")

        blockers = row["blocking_issues"]
        if not isinstance(blockers, list) or not blockers:
            raise SystemExit(f"{where}: blocking_issues must be nonempty")
        blocker_numbers: list[int] = []
        for blocker_index, blocker in enumerate(blockers):
            blocker_where = f"{where}.blocking_issues[{blocker_index}]"
            if not isinstance(blocker, dict):
                raise SystemExit(f"{blocker_where}: blocker must be an object")
            _exact_keys(
                blocker,
                {"number", "disposition", "role"},
                blocker_where,
            )
            if not _positive_int(blocker["number"]):
                raise SystemExit(f"{blocker_where}: number must be a positive integer")
            disposition = blocker["disposition"]
            if disposition not in {
                "open_work_item",
                "parked_computational_blocker",
                "closed_out_of_scope_blocker",
                "resource_deferred_blocker",
            }:
                raise SystemExit(
                    f"{blocker_where}: disposition is outside the controlled vocabulary"
                )
            if not isinstance(blocker["role"], str) or not blocker["role"].strip():
                raise SystemExit(f"{blocker_where}: role must be a nonempty string")
            if disposition == "closed_out_of_scope_blocker":
                if blocker["number"] not in closed_out_of_scope_issues:
                    raise SystemExit(
                        f"{blocker_where}: closed blocker is absent from the canonical "
                        "closed-out-of-scope ledger"
                    )
            elif (
                disposition != "resource_deferred_blocker"
                and blocker["number"] not in open_issues
            ):
                raise SystemExit(
                    f"{blocker_where}: open or parked blocker is absent from the "
                    "canonical open-issue ledger"
                )
            blocker_numbers.append(blocker["number"])
        if len(set(blocker_numbers)) != len(blocker_numbers):
            raise SystemExit(f"{where}: blocking issue numbers must be unique")

    selector_menus = accounting["selector_menus"]
    if not isinstance(selector_menus, list) or not selector_menus:
        raise SystemExit("physical identification registry: selector_menus must be nonempty")
    seen_selectors: set[str] = set()
    undeclared_blockers: set[int] = set()
    for index, menu in enumerate(selector_menus):
        where = f"selector_menus[{index}]"
        if not isinstance(menu, dict):
            raise SystemExit(f"{where}: row must be an object")
        _exact_keys(
            menu,
            {
                "selector_id",
                "label",
                "source",
                "status",
                "menu_size",
                "blocking_issues",
            },
            where,
        )
        selector_id = menu["selector_id"]
        if not isinstance(selector_id, str) or not selector_id:
            raise SystemExit(f"{where}: selector_id must be a nonempty string")
        if selector_id in seen_selectors:
            raise SystemExit(
                f"physical identification registry: duplicate selector_id {selector_id}"
            )
        seen_selectors.add(selector_id)
        if not isinstance(menu["label"], str) or not menu["label"].strip():
            raise SystemExit(f"{where}: label must be a nonempty string")
        source = menu["source"]
        if not isinstance(source, str) or not source:
            raise SystemExit(f"{where}: source must be a repository-relative path")
        relative_source = Path(source)
        if relative_source.is_absolute():
            raise SystemExit(f"{where}: source must be a repository-relative path")
        source_path = root / relative_source
        try:
            source_path.resolve().relative_to(root.resolve())
        except ValueError:
            raise SystemExit(f"{where}: source escapes the repository") from None
        if not source_path.is_file():
            raise SystemExit(f"{where}: source does not exist: {source}")

        blockers = menu["blocking_issues"]
        if (
            not isinstance(blockers, list)
            or any(not _positive_int(number) for number in blockers)
            or len(set(blockers)) != len(blockers)
        ):
            raise SystemExit(f"{where}: blocking_issues must be unique positive integers")
        missing_issues = sorted(set(blockers) - open_issues)
        if missing_issues:
            raise SystemExit(
                f"{where}: blocking issues are not open in the canonical snapshot "
                f"{missing_issues}"
            )

        status = menu["status"]
        if status == "declared":
            if not _positive_int(menu["menu_size"]):
                raise SystemExit(
                    f"{where}: a declared selector requires a positive integer menu_size"
                )
            if blockers:
                raise SystemExit(
                    f"{where}: a declared selector cannot retain blocking_issues"
                )
        elif status == "undeclared":
            if menu["menu_size"] is not None:
                raise SystemExit(
                    f"{where}: an undeclared selector cannot carry a numeric menu_size"
                )
            if not blockers:
                raise SystemExit(
                    f"{where}: an undeclared selector requires a live blocking issue"
                )
            undeclared_blockers.update(blockers)
        else:
            raise SystemExit(f"{where}: status must equal declared or undeclared")

    bound = accounting["compression_bound"]
    if not isinstance(bound, dict):
        raise SystemExit("compression_bound: value must be an object")
    _exact_keys(
        bound,
        {"status", "reason", "blocking_issues", "p_acc_upper_bound"},
        "compression_bound",
    )
    blockers = bound["blocking_issues"]
    if (
        not isinstance(blockers, list)
        or any(not _positive_int(number) for number in blockers)
        or len(set(blockers)) != len(blockers)
    ):
        raise SystemExit(
            "compression_bound: blocking_issues must be unique positive integers"
        )
    if set(blockers) - open_issues:
        raise SystemExit(
            "compression_bound: every blocking issue must be open in the canonical snapshot"
        )

    if undeclared_blockers:
        if bound["p_acc_upper_bound"] is not None:
            raise SystemExit(
                "compression_bound: numeric P_acc is forbidden while selector menu "
                "sizes are undeclared"
            )
        if bound["status"] != "not_computable" or bound["reason"] != "menu_sizes_undeclared":
            raise SystemExit(
                "compression_bound: undeclared selector menus require "
                "not_computable / menu_sizes_undeclared"
            )
        if set(blockers) != undeclared_blockers:
            raise SystemExit(
                "compression_bound: blockers must exactly match undeclared selector menus"
            )
    else:
        if blockers:
            raise SystemExit(
                "compression_bound: declared selector menus cannot retain blockers"
            )
        if bound["status"] == "not_computed":
            if (
                bound["reason"] != "calculation_not_run"
                or bound["p_acc_upper_bound"] is not None
            ):
                raise SystemExit(
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
                raise SystemExit(
                    "compression_bound: computed requires an all-menu numeric P_acc "
                    "between zero and one"
                )
        else:
            raise SystemExit(
                "compression_bound: fully declared menus require computed or not_computed"
            )


def compression_state(accounting: dict) -> str:
    bound = accounting["compression_bound"]
    if bound["status"] == "not_computable":
        issues = ", ".join(f"#{number}" for number in bound["blocking_issues"])
        return f"compression_bound: not_computable - menu sizes undeclared ({issues})"
    if bound["status"] == "not_computed":
        return "compression_bound: not_computed - calculation not run"
    return (
        "compression_bound: computed - "
        f"P_acc <= {bound['p_acc_upper_bound']}"
    )


def canonical_projection_payload(
    registry: dict,
    snapshot: dict,
    accounting: dict,
    root: Path = ROOT,
) -> dict:
    """Return the normalized status payload shared by all claim surfaces.

    The full open-issue-number set is included. It is broader than the issue
    subset referenced by current claims. Thus opening or closing a GitHub gate
    makes the proof spine and compression scorecard stale until the ledger and
    projection are regenerated together.
    """
    validate_accounting(registry, snapshot, accounting, root)
    return {
        "schema": 2,
        "release_id": registry["release_id"],
        "claims": registry["claims"],
        "open_issue_numbers": sorted(row["number"] for row in snapshot["rows"]),
        "physical_identification_accounting": accounting,
    }


def projection_digest(payload: dict) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def render_projection(
    registry: dict,
    snapshot: dict,
    accounting: dict,
    root: Path = ROOT,
) -> str:
    payload = canonical_projection_payload(registry, snapshot, accounting, root)
    claims = registry["claims"]
    by_class = {
        cls: sum(claim["claim_class"] == cls for claim in claims)
        for cls in CLASS_ORDER
    }
    distinct_gates = sorted(
        {gate for claim in claims for gate in claim["gates"]}
    )
    class_census = ", ".join(
        f"`{cls}` {by_class[cls]}" for cls in CLASS_ORDER
    )
    return "\n".join(
        [
            PROJECTION_START,
            "## Canonical claims and live-gates projection",
            "",
            "Generated by `tools/build_scoreboard.py`; do not edit this block. "
            "Its SHA-256 commits to every complete ordered canonical claim row "
            "(including class, status, and gates)",
            "and to the complete open-issue-number set in the committed "
            "GitHub ledger. The generated "
            "[claims scoreboard](../tracking/claims_scoreboard.md) is the "
            "row-by-row view.",
            "",
            f"- Release: `{registry['release_id']}`",
            f"- Projection SHA-256: `{projection_digest(payload)}`",
            f"- Registry: {len(claims)} canonical claims ({class_census})",
            f"- Claim gates: {len(distinct_gates)} distinct open issues across "
            f"{sum(bool(claim['gates']) for claim in claims)} gated claims",
            f"- Open-problem ledger: {len(snapshot['rows'])} open issues",
            f"- Undischarged physical identifications: "
            f"{len(accounting['physical_identifications'])}",
            f"- `{compression_state(accounting)}`",
            PROJECTION_END,
        ]
    )


def replace_projection(text: str, projection: str, relative_path: Path) -> str:
    if text.count(PROJECTION_START) != 1 or text.count(PROJECTION_END) != 1:
        raise SystemExit(
            f"{relative_path}: expected exactly one generated claim-status block; "
            "restore its BEGIN/END markers"
        )
    before, remainder = text.split(PROJECTION_START, 1)
    _, after = remainder.split(PROJECTION_END, 1)
    return before + projection + after


def render(
    registry: dict | None = None,
    snapshot: dict | None = None,
    accounting: dict | None = None,
    root: Path = ROOT,
) -> str:
    if registry is None or snapshot is None or accounting is None:
        registry, snapshot, accounting = source_documents(root)
    validate_accounting(registry, snapshot, accounting, root)
    claims = registry["claims"]
    repo = snapshot.get("repo", "")
    open_titles = {row["number"]: row.get("title", "") for row in snapshot.get("rows", [])}

    by_class: dict[str, int] = {}
    for claim in claims:
        by_class[claim["claim_class"]] = by_class.get(claim["claim_class"], 0) + 1
    gated = [claim for claim in claims if claim["gates"]]
    gate_numbers = sorted({gate for claim in gated for gate in claim["gates"]})

    lines = [
        "# OPH claims scoreboard",
        "",
        "Generated by `tools/build_scoreboard.py` from `claims/claim_registry.yaml`,",
        "`claims/physical_identification_registry.json`, and",
        "`tracking/open_issues/open_problem_ledger.json`. Do not edit by hand;",
        "`tools/run_mandatory_suite.py` fails when this page drifts from its sources.",
        "",
        "## Undischarged physical identifications",
        "",
        "Generated from `claims/physical_identification_registry.json`. Each row "
        "names the canonical claim surface, each ledger-tracked blocker, and the "
        "blocker's execution disposition.",
        "",
    ]
    for index, row in enumerate(accounting["physical_identifications"], start=1):
        claims_text = ", ".join(f"`{claim_id}`" for claim_id in row["claim_ids"])
        blockers = "; ".join(
            f"{_issue_link(repo, blocker['number'])} "
            f"[{blocker['disposition'].replace('_', ' ')}] {blocker['role']}"
            for blocker in row["blocking_issues"]
        )
        anchors = ", ".join(f"`{anchor}`" for anchor in row["source_anchors"])
        lines.append(
            f"{index}. **`{row['identification_id']}`: {row['label']}.** "
            f"Status: `{row['status']}`. Canonical claims: {claims_text}. "
            f"Source anchors: {anchors}. Blocking issues: {blockers}."
        )
    lines += [
        "",
        "## Compression-bound state",
        "",
        f"`{compression_state(accounting)}`",
        "",
        "A numeric `P_acc` or compression bound is inadmissible until every "
        "selector menu in the canonical registry has a declared finite size.",
        "",
        f"Release: `{registry['release_id']}`. {len(claims)} claims, "
        f"{len(gate_numbers)} distinct live gates across {len(gated)} gated claims.",
        "",
        "| Class | Claims |",
        "|---|---|",
    ]
    for cls in CLASS_ORDER:
        lines.append(f"| `{cls}` | {by_class.get(cls, 0)} |")
    lines += [
        "",
        "| Claim | Class | Status | Live gates |",
        "|---|---|---|---|",
    ]
    for claim in claims:
        gates = ", ".join(
            f"[#{gate}](https://github.com/{repo}/issues/{gate})" for gate in claim["gates"]
        ) or "none"
        lines.append(
            f"| `{claim['claim_id']}` | `{claim['claim_class']}` "
            f"| `{claim['status']}` | {gates} |"
        )
    lines += [
        "",
        "## Live gates",
        "",
        "| Issue | Title | Gated claims |",
        "|---|---|---|",
    ]
    for gate in gate_numbers:
        holders = ", ".join(
            f"`{claim['claim_id']}`" for claim in claims if gate in claim["gates"]
        )
        title = open_titles.get(gate, "").replace("|", "\\|")
        lines.append(
            f"| [#{gate}](https://github.com/{repo}/issues/{gate}) | {title} | {holders} |"
        )
    lines.append("")
    return "\n".join(lines)


def check_generated_surfaces(root: Path = ROOT) -> None:
    """Fail when any generated scoreboard surface differs from its sources."""

    registry, snapshot, accounting = source_documents(root)
    page = render(registry, snapshot, accounting, root)
    projection = render_projection(registry, snapshot, accounting, root)
    scoreboard_path = root / "tracking" / "claims_scoreboard.md"
    proof_spine_path = root / "docs" / "PROOF_SPINE.md"
    compression_path = root / "docs" / "COMPRESSION_SCORECARD.md"
    on_disk = (
        scoreboard_path.read_text(encoding="utf-8")
        if scoreboard_path.exists()
        else ""
    )
    if on_disk != page:
        raise SystemExit(
            f"{scoreboard_path.relative_to(root)} has drifted from its sources; "
            "regenerate with: python tools/build_scoreboard.py"
        )
    for path in [proof_spine_path, compression_path]:
        if not path.exists():
            raise SystemExit(
                f"{path.relative_to(root)} is missing; regenerate with: "
                "python tools/build_scoreboard.py"
            )
        current = path.read_text(encoding="utf-8")
        expected = replace_projection(current, projection, path.relative_to(root))
        if current != expected:
            raise SystemExit(
                f"{path.relative_to(root)} claim-status projection has drifted; "
                "regenerate with: python tools/build_scoreboard.py"
            )


def write_generated_surfaces(root: Path = ROOT) -> list[Path]:
    """Generate all scoreboard surfaces below an injectable repository root."""

    registry, snapshot, accounting = source_documents(root)
    page = render(registry, snapshot, accounting, root)
    projection = render_projection(registry, snapshot, accounting, root)
    scoreboard_path = root / "tracking" / "claims_scoreboard.md"
    proof_spine_path = root / "docs" / "PROOF_SPINE.md"
    compression_path = root / "docs" / "COMPRESSION_SCORECARD.md"
    projection_surfaces = [
        (proof_spine_path, None),
        (compression_path, COMPRESSION_TEMPLATE),
    ]
    scoreboard_path.write_text(page, encoding="utf-8")
    written = [scoreboard_path.relative_to(root)]
    for path, template in projection_surfaces:
        current = path.read_text(encoding="utf-8") if path.exists() else template
        if current is None:
            raise SystemExit(f"{path.relative_to(root)} is missing")
        path.write_text(
            replace_projection(current, projection, path.relative_to(root)),
            encoding="utf-8",
        )
        written.append(path.relative_to(root))
    return written


def main() -> None:
    if "--check" in sys.argv[1:]:
        check_generated_surfaces()
        print(
            "claims scoreboard, proof spine, and compression scorecard OK: "
            "in sync with registries and issue snapshot"
        )
        return
    written = write_generated_surfaces()
    print("wrote " + ", ".join(map(str, written)))


if __name__ == "__main__":
    main()
