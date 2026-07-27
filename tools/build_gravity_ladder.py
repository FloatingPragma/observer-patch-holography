"""Build and validate the gravity premise-elimination ladder (issue #618).

The machine-readable ladder is ``claims/gravity_premise_ladder.json``. This
tool validates it fail-closed and renders ``docs/GRAVITY_PREMISE_LADDER.md``;
``--check`` fails when the committed page differs from the render, and the
mandatory suite runs that check.

Fail-closed rules: the rungs are contiguous one through eleven; every rung
names at least one theorem-or-countermodel reference and at least one
deletion or mutation test, and every reference resolves (Lean declarations
present in their modules, code symbols present in their files, receipt
paths existing, paper anchors present in their TeX sources, test names
present in their test files); every premise interface id exists in the
axiom registry with the rung status drawn from its interface classes, so a
rung cannot claim a status stronger than its registry classification; every
assumption token resolves in the assumption dictionary; every claim id
resolves in the claim registry; and pending physical-identification rungs
name live owning issues.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
LADDER_PATH = ROOT / "claims" / "gravity_premise_ladder.json"
SURFACE_PATH = ROOT / "docs" / "GRAVITY_PREMISE_LADDER.md"
AXIOM_REGISTRY_PATH = ROOT / "claims" / "axiom_registry.yaml"
CLAIM_REGISTRY_PATH = ROOT / "claims" / "claim_registry.yaml"
ASSUMPTION_PATH = ROOT / "claims" / "assumption_dictionary.md"
SNAPSHOT_PATH = ROOT / "tracking" / "open_issues" / "open_problem_ledger.json"

SCHEMA = "oph.gravity_premise_ladder.v1"

RUNG_KEYS = {
    "rung",
    "title",
    "status",
    "premise_interfaces",
    "assumption_tokens",
    "theorem_or_countermodel",
    "mutation_tests",
    "paper_anchors",
    "claim_ids",
    "owner_issues",
}


def fail(message: str) -> None:
    raise SystemExit(f"gravity premise ladder: {message}")


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        fail(f"missing input {path.relative_to(ROOT)}")
    except json.JSONDecodeError as error:
        fail(f"invalid JSON in {path.relative_to(ROOT)}: {error}")
    raise AssertionError("unreachable")


def validate(ladder: dict) -> list[dict]:
    if ladder.get("schema") != SCHEMA:
        fail(f"schema must equal {SCHEMA}")
    if ladder.get("issue") != 618:
        fail("issue must equal 618")
    rungs = ladder.get("rungs")
    if not isinstance(rungs, list) or [row.get("rung") for row in rungs] != list(
        range(1, 12)
    ):
        fail("rungs must be contiguous from one through eleven")

    axiom_registry = yaml.safe_load(AXIOM_REGISTRY_PATH.read_text(encoding="utf-8"))
    status_enum = set(axiom_registry["status_enum"])
    interface_classes = {
        row["id"]: row["class"] for row in axiom_registry["premise_interfaces"]
    }
    claim_registry = load_json(CLAIM_REGISTRY_PATH)
    claim_ids = {claim["claim_id"] for claim in claim_registry["claims"]}
    assumptions = ASSUMPTION_PATH.read_text(encoding="utf-8")
    snapshot = load_json(SNAPSHOT_PATH)
    open_issues = {row["number"] for row in snapshot["rows"]}

    for row in rungs:
        where = f"rung {row.get('rung')}"
        if set(row) != RUNG_KEYS:
            fail(f"{where}: keys mismatch")
        if not isinstance(row["title"], str) or not row["title"].strip():
            fail(f"{where}: title must be nonempty")
        if row["status"] not in status_enum:
            fail(f"{where}: status {row['status']} is not in the canonical enum")
        interfaces = row["premise_interfaces"]
        if not isinstance(interfaces, list) or not interfaces:
            fail(f"{where}: at least one premise interface is required")
        for interface in interfaces:
            if interface not in interface_classes:
                fail(f"{where}: unknown premise interface {interface}")
        allowed = {interface_classes[i] for i in interfaces}
        if row["status"] not in allowed:
            fail(
                f"{where}: status {row['status']} is not among its registry "
                f"interface classes {sorted(allowed)}"
            )
        for token in row["assumption_tokens"]:
            if f"`{token}`" not in assumptions:
                fail(f"{where}: assumption token {token} is not in the dictionary")
        refs = row["theorem_or_countermodel"]
        if not isinstance(refs, list) or not refs:
            fail(f"{where}: a theorem or countermodel reference is required")
        for ref in refs:
            kind = ref.get("kind")
            if kind == "lean":
                module = ROOT / ref["module"]
                if not module.is_file():
                    fail(f"{where}: Lean module missing: {ref['module']}")
                if ref["declaration"] not in module.read_text(encoding="utf-8"):
                    fail(
                        f"{where}: declaration {ref['declaration']} absent from "
                        f"{ref['module']}"
                    )
            elif kind == "code":
                path = ROOT / ref["path"]
                if not path.is_file():
                    fail(f"{where}: code path missing: {ref['path']}")
                if ref["symbol"] not in path.read_text(encoding="utf-8"):
                    fail(
                        f"{where}: symbol {ref['symbol']} absent from {ref['path']}"
                    )
            elif kind == "receipt":
                if not (ROOT / ref["path"]).is_file():
                    fail(f"{where}: receipt missing: {ref['path']}")
            elif kind == "paper":
                paper = ROOT / ref["file"]
                if not paper.is_file():
                    fail(f"{where}: paper missing: {ref['file']}")
                text = paper.read_text(encoding="utf-8")
                if (
                    f"\\label{{{ref['anchor']}}}" not in text
                    and ref["anchor"] not in text
                ):
                    fail(
                        f"{where}: anchor {ref['anchor']} absent from {ref['file']}"
                    )
            else:
                fail(f"{where}: unknown reference kind {kind}")
        tests = row["mutation_tests"]
        if not isinstance(tests, list) or not tests:
            fail(f"{where}: a deletion or mutation test is required")
        for test in tests:
            path = ROOT / test["path"]
            if not path.is_file():
                fail(f"{where}: test file missing: {test['path']}")
            if test["name"] not in path.read_text(encoding="utf-8"):
                fail(f"{where}: test {test['name']} absent from {test['path']}")
        for anchor in row["paper_anchors"]:
            paper = ROOT / anchor["file"]
            if not paper.is_file():
                fail(f"{where}: paper missing: {anchor['file']}")
            text = paper.read_text(encoding="utf-8")
            if (
                f"\\label{{{anchor['anchor']}}}" not in text
                and anchor["anchor"] not in text
            ):
                fail(f"{where}: anchor {anchor['anchor']} absent from {anchor['file']}")
        for claim_id in row["claim_ids"]:
            if claim_id not in claim_ids:
                fail(f"{where}: unknown claim id {claim_id}")
        for issue in row["owner_issues"]:
            if issue not in open_issues:
                fail(f"{where}: owner issue #{issue} is not open in the snapshot")
        if row["status"] == "physical_identification" and not row["owner_issues"]:
            fail(f"{where}: a pending physical identification requires owner issues")
    return rungs


def render(ladder: dict, rungs: list[dict]) -> str:
    lines: list[str] = []
    lines.append("# The gravity premise-elimination ladder")
    lines.append("")
    lines.append(
        "Generated by `tools/build_gravity_ladder.py` from"
        " `claims/gravity_premise_ladder.json`; edit the JSON, then"
        " regenerate. Issue #618 owns this ladder."
    )
    lines.append("")
    lines.append(ladder["policy"])
    lines.append("")
    lines.append("| Rung | Premise | Status | Theorem or countermodel | Deletion/mutation test | Owner |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for row in rungs:
        refs = "; ".join(
            ref.get("declaration") or ref.get("symbol") or ref.get("anchor") or Path(ref.get("path", "")).name
            for ref in row["theorem_or_countermodel"]
        )
        tests = "; ".join(
            f"`{Path(test['path']).name}`" for test in row["mutation_tests"]
        )
        owner = (
            ", ".join(
                f"[#{n}](https://github.com/FloatingPragma/observer-patch-holography/issues/{n})"
                for n in row["owner_issues"]
            )
            if row["owner_issues"]
            else "ladder"
        )
        lines.append(
            f"| {row['rung']} | {row['title']} | {row['status']} | {refs} |"
            f" {tests} | {owner} |"
        )
    lines.append("")
    lines.append("## Premise interfaces and assumption tokens per rung")
    lines.append("")
    for row in rungs:
        interfaces = ", ".join(f"`{i}`" for i in row["premise_interfaces"])
        tokens = (
            ", ".join(f"`{t}`" for t in row["assumption_tokens"])
            if row["assumption_tokens"]
            else "none beyond the interfaces"
        )
        claims = (
            ", ".join(f"`{c}`" for c in row["claim_ids"])
            if row["claim_ids"]
            else "carried by the interface registry"
        )
        lines.append(
            f"- **Rung {row['rung']}**: interfaces {interfaces}; assumption"
            f" tokens {tokens}; claim rows {claims}."
        )
    lines.append("")
    lines.append(
        "Rungs one and two are delivered by the equal-state-weights and"
    )
    lines.append(
        "optimizer-pushforward receipts, with the cross-refinement tail a"
    )
    lines.append(
        "named interface. Every rung status equals an axiom-registry"
    )
    lines.append(
        "interface class of that rung, and the ladder"
    )
    lines.append(
        "validation fails closed when any reference, test, token, claim, or"
    )
    lines.append("owner issue stops resolving.")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    ladder = load_json(LADDER_PATH)
    rungs = validate(ladder)
    surface = render(ladder, rungs)
    if args.check:
        committed = (
            SURFACE_PATH.read_text(encoding="utf-8") if SURFACE_PATH.is_file() else ""
        )
        if committed != surface:
            print(
                "gravity premise ladder: docs/GRAVITY_PREMISE_LADDER.md is stale; "
                "run python tools/build_gravity_ladder.py",
                file=sys.stderr,
            )
            return 1
        print("gravity premise ladder: surface is current")
        return 0
    SURFACE_PATH.write_text(surface, encoding="utf-8", newline="\n")
    print(f"gravity premise ladder: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
