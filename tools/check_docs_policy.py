"""Fail when a document under docs/ carries dated or status-record material.

The rule set mirrors AGENTS.md and docs/README.md: no dates in file names or
headings, no progress or audit-record prose, no pull-request handoffs, and no
root document outside the declared set of programs, references and ledgers.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DOCS = "docs"

# Root documents are the ones docs/README.md declares, plus the index itself.
ROOT_INDEX_LINK = re.compile(r"\]\(([A-Za-z0-9_.-]+\.md)\)")

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


def changed_documents(root: Path, base: str) -> list[Path]:
    """Documents added, copied, modified or renamed against a git revision."""
    result = subprocess.run(
        [
            "git", "-C", str(root), "diff", "--name-only", "--diff-filter=ACMR",
            f"{base}...HEAD",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise SystemExit(f"docs policy: cannot read the diff against {base}: "
                         f"{result.stderr.strip()}")
    documents = []
    for line in result.stdout.splitlines():
        if " => " in line:
            continue
        path = (root / line).resolve()
        if path.is_file() and path.is_relative_to((root / DOCS).resolve()):
            documents.append(path)
    return sorted(documents)


def declared_root_documents(root: Path) -> set[str]:
    """Root documents that the docs index declares, plus the index itself."""
    index = root / DOCS / "README.md"
    names = {"README.md"}
    if index.is_file():
        names.update(
            ROOT_INDEX_LINK.findall(index.read_text(encoding="utf-8"))
        )
    return names


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


def check(
    root: Path = ROOT, documents: list[Path] | None = None
) -> tuple[list[Violation], list[Violation]]:
    """Return (unexempted violations, violations covered by a known correction).

    With ``documents`` given, only those files and their ``docs/`` ancestors
    are scanned; that is the changed-file mode a pull request runs.
    """
    root = root.resolve()
    docs = root / DOCS
    targets = (
        _document_paths(root)
        if documents is None
        else [path.resolve() for path in documents if path.is_file()]
    )
    declared_root = declared_root_documents(root)
    found: list[Violation] = []
    for path in targets:
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
            if path.parent == docs and path.name not in declared_root:
                found.append(Violation(relative, None, "root-document", path.name))
    directories = (
        sorted(path for path in docs.rglob("*") if path.is_dir())
        if documents is None
        else sorted({
            parent
            for path in targets
            for parent in path.parents
            if parent.is_relative_to(docs) and parent != docs
        })
    )
    for directory in directories:
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


def _report(violations: list[Violation], prefix: str = "") -> None:
    for violation in violations:
        location = violation.path + (f":{violation.line}" if violation.line else "")
        print(f"{prefix}{location}: {violation.rule}: {violation.text}")


def main(root: Path = ROOT, base: str | None = None) -> None:
    documents = None if base is None else changed_documents(root, base)
    failures, exempted = check(root, documents)
    _report(failures)
    _report(exempted, prefix="exempt: ")
    scope = (
        f"{len(documents)} changed file(s) against {base}"
        if documents is not None
        else f"{len(_document_paths(root.resolve()))} files"
    )
    failures_by_rule = Counter(v.rule for v in failures)
    exempted_by_rule = Counter(v.rule for v in exempted)
    if failures:
        counts = ", ".join(
            f"{rule}={count}" for rule, count in sorted(failures_by_rule.items())
        )
        raise SystemExit(f"docs policy: {len(failures)} violation(s) in {scope}: {counts}")
    known = ", ".join(
        f"{rule}={count}" for rule, count in sorted(exempted_by_rule.items())
    )
    print(
        f"docs policy OK: {scope}, "
        f"{len(exempted)} known violation(s) exempted"
        + (f" ({known})" if known else "")
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="scan the whole tree (the default)",
    )
    parser.add_argument(
        "--changed",
        metavar="BASE",
        help="scan only documents added, copied, modified or renamed against BASE",
    )
    arguments = parser.parse_args()
    main(arguments.root, arguments.changed)
