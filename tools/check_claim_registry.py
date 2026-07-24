#!/usr/bin/env python3
"""Validate the OPH claim-registry seed files.

The registry file is JSON-compatible YAML, so this validator avoids an external
YAML dependency while keeping the requested `.yaml` public path.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

PAPER_EXTERNAL_REGISTRY_PATTERNS = [
    "claims/claim_registry",
    "\\ophid{claims/",
]

REQUIRED_CLAIM_FIELDS = {
    "claim_id",
    "statement",
    "owner_paper",
    "tier",
    "assumptions",
    "imported_results",
    "oph_specific_delta",
    "novelty_type",
    "evidence",
    "falsifier",
    "scope_if_false",
    "status",
    "claim_class",
    "gates",
}

# Controlled claim classification (issue #512). `status` stays descriptive
# free text; `claim_class` is the machine-checked epistemic class.
CLAIM_CLASS_VOCABULARY = {
    # Definitions, contracts, programs, interpretive frames; no theorem or
    # executed artifact is asserted.
    "declared_structure",
    # A proved implication whose premises include declared (not source-derived)
    # inputs; no claim that the physical world instantiates the branch.
    "conditional_implication",
    # The assertion that a declared physical branch is entered or nonempty.
    "branch_entry",
    # An executed pipeline artifact (audit, archive, selector run, protocol
    # record); the claim is the artifact's existence and content.
    "emitted_artifact",
    # A quantitative pipeline landing compared against measurement, with a
    # declared endpoint boundary or error budget.
    "empirical_implementation",
    # Physically established end to end; requires every live gate closed.
    "physical_establishment",
}

# Classes whose public wording may assert physical establishment. Every other
# class fails the wording gate if its statement claims establishment.
PROMOTED_CLASSES = {"physical_establishment"}

ESTABLISHMENT_WORDING = re.compile(
    r"physically established|experimentally confirmed|physically confirmed"
    r"|empirically established|is established as physical",
    re.IGNORECASE,
)

ISSUE_SNAPSHOT_RELATIVE = Path("tracking") / "open_issues" / "open_problem_ledger.json"
CLAIM_ID_TOKEN = re.compile(r"\bOPH-[A-Z0-9][A-Z0-9-]{2,}\b")


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SystemExit(f"{path}: invalid JSON-compatible YAML: {exc}") from exc


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def release_id_from_tex(root: Path) -> str:
    text = (root / "paper" / "release_info.tex").read_text(encoding="utf-8")
    match = re.search(r"\\newcommand\{\\OPHPaperReleaseID\}\{([^}]+)\}", text)
    require(match is not None, "paper/release_info.tex does not define OPHPaperReleaseID")
    return match.group(1)


def check_standalone_papers(root: Path) -> None:
    for folder in ["paper", "extra"]:
        for path in (root / folder).glob("*.tex"):
            text = path.read_text(encoding="utf-8")
            for pattern in PAPER_EXTERNAL_REGISTRY_PATTERNS:
                require(
                    pattern not in text,
                    f"{path.relative_to(root)} references the external claim registry; papers must remain standalone",
                )


def dictionary_tokens(root: Path) -> set[str]:
    """Assumption tokens with a canonical dictionary row (backtick-quoted)."""
    text = (root / "claims" / "assumption_dictionary.md").read_text(encoding="utf-8")
    return set(re.findall(r"^\|\s*`([^`]+)`", text, re.MULTILINE))


def load_issue_snapshot(root: Path) -> tuple[set[int], set[int], dict]:
    """Open and closed issue numbers from the committed GitHub issue snapshot.

    The snapshot is the fail-closed stand-in for live GitHub state: gates are
    validated against it, and CI regenerating it against the live repository is
    what keeps it honest. The snapshot itself must be internally consistent —
    a hard-coded or stale `open_issue_count` is rejected here.
    """
    snapshot = load_json(root / ISSUE_SNAPSHOT_RELATIVE)
    rows = snapshot.get("rows", [])
    numbers = [row.get("number") for row in rows]
    require(
        all(isinstance(n, int) for n in numbers),
        f"{ISSUE_SNAPSHOT_RELATIVE}: issue rows with non-integer numbers",
    )
    require(
        len(set(numbers)) == len(numbers),
        f"{ISSUE_SNAPSHOT_RELATIVE}: duplicate issue numbers in rows",
    )
    require(
        snapshot.get("open_issue_count") == len(rows),
        f"{ISSUE_SNAPSHOT_RELATIVE}: open_issue_count "
        f"{snapshot.get('open_issue_count')!r} does not equal the computed row "
        f"count {len(rows)} (hard-coded or stale count)",
    )
    closed = {row.get("number") for row in snapshot.get("closed_out_of_scope_records", [])}
    return set(numbers), closed, snapshot


def check_gates(claim: dict, open_issues: set[int], closed_issues: set[int]) -> None:
    """Fail closed on stale, missing, or promotion-violating live gates."""
    claim_id = claim["claim_id"]
    gates = claim["gates"]
    require(
        isinstance(gates, list) and all(isinstance(g, int) for g in gates),
        f"{claim_id}: gates must be a list of GitHub issue numbers",
    )
    for gate in gates:
        require(
            gate not in closed_issues,
            f"{claim_id}: gate #{gate} is closed on GitHub but still referenced "
            "as an open dependency",
        )
        require(
            gate in open_issues,
            f"{claim_id}: gate #{gate} is missing from the GitHub issue snapshot "
            "(open theorem work not tracked on GitHub)",
        )
    if claim["claim_class"] in PROMOTED_CLASSES:
        require(
            not gates,
            f"{claim_id}: claim_class {claim['claim_class']!r} asserts physical "
            f"establishment while gates {gates} are still open",
        )


def check_wording(claim: dict) -> None:
    """Public wording must not be stronger than the claim's class."""
    if claim["claim_class"] in PROMOTED_CLASSES:
        return
    match = ESTABLISHMENT_WORDING.search(claim["statement"])
    if match is not None:
        raise SystemExit(
            f"{claim['claim_id']}: statement wording {match.group(0)!r} asserts "
            f"establishment but claim_class is {claim['claim_class']!r}"
        )


def check_snapshot_claim_tokens(snapshot: dict, known_claims: set[str]) -> None:
    """Reverse direction: claim IDs cited by the issue snapshot must exist."""
    for row in snapshot.get("rows", []) + snapshot.get("closed_out_of_scope_records", []):
        for value in row.values():
            if not isinstance(value, str):
                continue
            for token in CLAIM_ID_TOKEN.findall(value):
                require(
                    token in known_claims,
                    f"{ISSUE_SNAPSHOT_RELATIVE}: issue #{row.get('number')} cites "
                    f"unknown claim id {token}",
                )


def main(root: Path = ROOT) -> None:
    registry = load_json(root / "claims" / "claim_registry.yaml")
    require(
        registry.get("release_id") == release_id_from_tex(root),
        f"registry release_id {registry.get('release_id')!r} does not match paper/release_info.tex",
    )
    claims = registry.get("claims", [])
    require(isinstance(claims, list) and claims, "claim registry has no claims")
    check_standalone_papers(root)
    defined_tokens = dictionary_tokens(root)
    open_issues, closed_issues, snapshot = load_issue_snapshot(root)

    seen: set[str] = set()
    owner_paths: set[str] = set()
    for claim in claims:
        missing = REQUIRED_CLAIM_FIELDS - set(claim)
        require(not missing, f"{claim.get('claim_id', '<missing>')}: missing fields {sorted(missing)}")
        claim_id = claim["claim_id"]
        require(claim_id not in seen, f"duplicate claim_id {claim_id}")
        seen.add(claim_id)
        require(claim["statement"].strip(), f"{claim_id}: empty statement")
        require(claim["assumptions"], f"{claim_id}: empty assumptions")
        require(claim["imported_results"], f"{claim_id}: empty imported_results")
        require(claim["oph_specific_delta"].strip(), f"{claim_id}: empty OPH delta")
        require(claim["falsifier"].strip(), f"{claim_id}: empty falsifier")
        require(
            claim["claim_class"] in CLAIM_CLASS_VOCABULARY,
            f"{claim_id}: claim_class {claim.get('claim_class')!r} is not in the "
            f"controlled vocabulary {sorted(CLAIM_CLASS_VOCABULARY)}",
        )
        check_gates(claim, open_issues, closed_issues)
        check_wording(claim)
        owner = root / claim["owner_paper"]
        require(owner.exists(), f"{claim_id}: owner paper does not exist: {claim['owner_paper']}")
        owner_paths.add(claim["owner_paper"])
        undefined = sorted(set(claim["assumptions"]) - defined_tokens)
        require(
            not undefined,
            f"{claim_id}: assumption tokens without a canonical dictionary row: {undefined}",
        )
        for evidence in claim["evidence"]:
            require(
                (root / evidence).exists(),
                f"{claim_id}: evidence path does not exist: {evidence}",
            )

    check_snapshot_claim_tokens(snapshot, seen)

    for matrix_name, required_columns in [
        ("novelty_matrix.csv", {"claim_id", "closest_prior_work", "oph_specific_delta", "novelty_type", "falsifier"}),
        ("falsification_matrix.csv", {"claim_id", "mathematical_falsifier", "physical_identification_falsifier", "phenomenological_falsifier", "scope_if_false"}),
    ]:
        matrix_path = root / "claims" / matrix_name
        rows = load_csv(matrix_path)
        require(rows, f"{matrix_path}: no rows")
        require(required_columns.issubset(rows[0].keys()), f"{matrix_path}: missing required columns")
        for row in rows:
            claim_id = row["claim_id"]
            require(claim_id in seen, f"{matrix_path}: unknown claim_id {claim_id}")

    graph = load_json(root / "claims" / "dependency_graph.json")
    nodes = set(graph.get("nodes", []))
    require(nodes <= seen, f"dependency graph has unknown nodes: {sorted(nodes - seen)}")
    for edge in graph.get("edges", []):
        require(edge.get("from") in seen, f"dependency graph has unknown edge source: {edge}")
        require(edge.get("to") in seen, f"dependency graph has unknown edge target: {edge}")
        require(edge.get("role"), f"dependency graph edge lacks role: {edge}")

    gated = [claim for claim in claims if claim["gates"]]
    gate_count = len({gate for claim in gated for gate in claim["gates"]})
    print(
        f"claim registry OK: {len(seen)} claims, {len(owner_paths)} owner papers, "
        f"{gate_count} live gates across {len(gated)} gated claims"
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT)
