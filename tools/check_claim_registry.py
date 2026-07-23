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
}


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

    print(f"claim registry OK: {len(seen)} claims, {len(owner_paths)} owner papers")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT)
