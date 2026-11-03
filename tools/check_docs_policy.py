"""Fail when a document under docs/ carries dated or status-record material.

The rule set mirrors AGENTS.md and docs/README.md: no dates in file names or
headings, no progress or audit-record prose, no pull-request handoffs, and no
root document outside the declared set of programs, references and ledgers.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DOCS = "docs"

# Root documents the repository declares as programs, references and ledgers.
ALLOWED_ROOT_DOCUMENTS = frozenset({
    "README.md",
    "AGREEMENT_RESEARCH_FOUNDATIONS.md",
    "AXIOM_REFERENCE.md",
    "BASIS.md",
    "COMMON_OBJECTIONS.md",
    "COMPLETION_PROGRAM.md",
    "EMPIRICAL_STATUS.md",
    "FROZEN_PREDICTION_LADDER.md",
    "INTERPRETATION_AND_TERMS.md",
    "OPH_EXECUTABLE_HYPOTHESIS.md",
    "OPH_FALSIFICATION_PROGRAM.md",
    "PAPER_SET.md",
    "POSTDICTION_LEDGER.md",
})

# The documented exception: retrospective public-data comparisons that no
# single code package reproduces.
EXEMPT_DOCUMENTS = {
    "docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md":
        "public-data comparison exception",
}

# Violations that open corrections remove; drop each entry with its correction.
KNOWN_VIOLATIONS = {
    "docs/comms/": {
        "dated-file-name": "archive tree leaves docs/ with its relocation",
        "dated-directory-name": "archive tree leaves docs/ with its relocation",
        "dated-heading": "archive tree leaves docs/ with its relocation",
        "status-prose": "archive tree leaves docs/ with its relocation",
    },
    "docs/policies/NUMERICAL_AUDIT.md": {
        "dated-heading": "review history leaves the policy with its rewrite",
        "pull-request-reference": "review history leaves the policy with its rewrite",
    },
    "docs/research/MARKED_TRIPLET_OPERATIONAL_ATTACHMENT_AUDIT.md": {
        "pull-request-reference": "audit history leaves the note with its rewrite",
        "status-prose": "audit history leaves the note with its rewrite",
    },
    "docs/EXECUTABLE_PROGRAM_CLASS_V2_PROPOSAL.md": {
        "root-document": "proposal moves beside its consumer",
    },
}

ISO_DATE = re.compile(
    r"(?<!\d)(?:19|20)\d{2}[-_](?:0[1-9]|1[0-2])[-_](?:0[1-9]|[12]\d|3[01])(?!\d)"
)
MONTHS = (
    "January|February|March|April|May|June|July|August|September|October|"
    "November|December"
)
MONTH_DATE = re.compile(rf"\b(?:\d{{1,2}}\s+)?(?:{MONTHS})\s+\d{{4}}\b")
HEADING = re.compile(r"^#{1,6}\s")
PULL_REQUEST = re.compile(r"\bPR #\d+|/pull/\d+")
STATUS_MARKERS = (
    "**audit date:**",
    "baseline commit",
    "generated_utc",
)
CHRONOLOGY_HEADING = re.compile(
    r"^#{1,6}\s+(?:\d+\.\s*)?(?:Chronology|Timeline|History)\b", re.IGNORECASE
)


@dataclass(frozen=True)
class Violation:
    path: str
    line: int | None
    rule: str
    text: str


def _document_paths(root: Path) -> list[Path]:
    docs = root / DOCS
    if not docs.is_dir():
        return []
    return sorted(path for path in docs.rglob("*") if path.is_file())


def _text_violations(relative: str, lines: list[str]) -> list[Violation]:
    found: list[Violation] = []
    for number, line in enumerate(lines, start=1):
        stripped = line.strip()
        lowered = line.lower()
        if HEADING.match(line) and (ISO_DATE.search(line) or MONTH_DATE.search(line)):
            found.append(Violation(relative, number, "dated-heading", stripped))
        if PULL_REQUEST.search(line):
            found.append(Violation(relative, number, "pull-request-reference", stripped))
        if any(marker in lowered for marker in STATUS_MARKERS):
            found.append(Violation(relative, number, "status-prose", stripped))
        if CHRONOLOGY_HEADING.match(line):
            found.append(Violation(relative, number, "chronology-heading", stripped))
    return found


def check(root: Path = ROOT) -> tuple[list[Violation], list[Violation]]:
    """Return (unexempted violations, violations covered by a known correction)."""
    root = root.resolve()
    found: list[Violation] = []
    for path in _document_paths(root):
        relative = path.relative_to(root).as_posix()
        if path.suffix == ".md" and relative in EXEMPT_DOCUMENTS:
            continue
        if ISO_DATE.search(path.name):
            found.append(Violation(relative, None, "dated-file-name", path.name))
        if path.suffix == ".md":
            found.extend(
                _text_violations(
                    relative, path.read_text(encoding="utf-8").splitlines()
                )
            )
            if path.parent == root / DOCS and path.name not in ALLOWED_ROOT_DOCUMENTS:
                found.append(Violation(relative, None, "root-document", path.name))
    for directory in sorted(
        path for path in (root / DOCS).rglob("*") if path.is_dir()
    ):
        if ISO_DATE.search(directory.name):
            found.append(
                Violation(
                    directory.relative_to(root).as_posix(),
                    None,
                    "dated-directory-name",
                    directory.name,
                )
            )
    failures: list[Violation] = []
    exempted: list[Violation] = []
    for violation in found:
        known = KNOWN_VIOLATIONS.get(violation.path) or next(
            (
                rules
                for prefix, rules in KNOWN_VIOLATIONS.items()
                if violation.path.startswith(prefix)
            ),
            {},
        )
        if violation.rule in known:
            exempted.append(violation)
        else:
            failures.append(violation)
    return failures, exempted


def main(root: Path = ROOT) -> None:
    failures, exempted = check(root)
    for violation in failures:
        location = violation.path + (f":{violation.line}" if violation.line else "")
        print(f"{location}: {violation.rule}: {violation.text}")
    for violation in exempted:
        location = violation.path + (f":{violation.line}" if violation.line else "")
        print(f"exempt: {location}: {violation.rule}")
    documents = len(_document_paths(root.resolve()))
    if failures:
        raise SystemExit(
            f"docs policy: {len(failures)} violation(s) across {documents} files"
        )
    print(
        f"docs policy OK: {documents} files checked, "
        f"{len(exempted)} known violation(s) exempted"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    main(parser.parse_args().root)
