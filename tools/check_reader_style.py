#!/usr/bin/env python3
"""Check reader-facing OPH style constraints.

The check is intentionally narrow. It enforces the durable rules that should
not depend on taste: no changelog/progress framing, no common AI-prose tells
called out in the style guide for this release, no banned h-word, and no internal
implementation identifiers in paper abstracts.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

READER_GLOBS = [
    "README.md",
    "README_FR.md",
    "REPRODUCE.md",
    "Lean/README.md",
    "Lean/docs/LIBRARY_GUIDE.md",
    "code/README.md",
    "docs/README.md",
    "docs/APPLICATIONS.md",
    "docs/COMMON_OBJECTIONS.md",
    "docs/ENTANGLEMENT_GEOMETRY_PROBLEM_OPH.md",
    "docs/STRANGE_LOOP.md",
    "assets/README.md",
    "paper/README.md",
    "extra/README.md",
    "book/**/*.md",
    "claims/*.md",
    "cosmology/**/*.md",
    "physics-problems/**/*.md",
]

STATUS_GLOBS = [
    "docs/PROOF_SPINE.md",
    "docs/BOOK_CHAPTER_LEDGER.md",
    "docs/CONSISTENCY_STACK.md",
    "docs/CLOSURE_LEDGER.md",
    "docs/SURVIVAL_PROOF_FINAL_INTEGRATION_*.md",
    "docs/UNIFIED_CARRIER_COHERENCE_AUDIT_*.md",
    "code/geometry/*STATUS*.md",
    "claims/claim_registry.yaml",
    "claims/dependency_graph.json",
    "claims/falsification_matrix.csv",
    "claims/novelty_matrix.csv",
    "code/geometry/runs/*status.json",
]

INFORMAL_GLOBS = [
    "README.md",
    "README_FR.md",
    "book/**/*.md",
]

PAPER_GLOBS = [
    "paper/**/*.tex",
    "extra/**/*.tex",
    "cosmology/**/*.tex",
]

PROGRESS_PATTERNS = [
    (re.compile(r"\bnow\b", re.IGNORECASE), "progress word: now"),
    (re.compile(r"\balready\b", re.IGNORECASE), "progress word: already"),
    (re.compile(r"\bcurrently\b", re.IGNORECASE), "progress word: currently"),
    (re.compile(r"\bpresently\b", re.IGNORECASE), "progress word: presently"),
    (re.compile(r"\blatest\b", re.IGNORECASE), "progress word: latest"),
    (re.compile(r"\brecent(?:ly)?\b", re.IGNORECASE), "progress word: recent/recently"),
    (re.compile(r"\bnewly\b", re.IGNORECASE), "progress word: newly"),
    (re.compile(r"\bpreviously\b", re.IGNORECASE), "progress word: previously"),
    (re.compile(r"\bformerly\b", re.IGNORECASE), "progress word: formerly"),
    (re.compile(r"\bno longer\b", re.IGNORECASE), "progress phrase: no longer"),
    (re.compile(r"\bstill\b", re.IGNORECASE), "progress word: still"),
    (re.compile(r"\b(?:has|have|had)?\s*not yet\b", re.IGNORECASE), "progress phrase: not yet"),
    (re.compile(r"\bremains? open\b", re.IGNORECASE), "progress phrase: remains open"),
    (re.compile(r"\bso far\b", re.IGNORECASE), "progress phrase: so far"),
    (re.compile(r"\bfuture work\b", re.IGNORECASE), "progress phrase: future work"),
    (re.compile(r"\bnext step\b", re.IGNORECASE), "progress phrase: next step"),
    (re.compile(r"\bcurrent state\b", re.IGNORECASE), "progress phrase: current state"),
    (re.compile(r"\bwhat changed\b", re.IGNORECASE), "progress phrase: what changed"),
    (re.compile(r"\bwhat did not change\b", re.IGNORECASE), "progress phrase: what did not change"),
]

AI_TELL_PATTERNS = [
    (re.compile(r"\u2014|&mdash;"), "em dash"),
    (re.compile(r"\bnot\s+only\b", re.IGNORECASE), "not-only contrast"),
    (re.compile(r"\bnot\s+just\b", re.IGNORECASE), "not-just contrast"),
    (re.compile(r"\bnot\s+merely\b", re.IGNORECASE), "not-merely contrast"),
    (re.compile(r"\bnot\s+simply\b", re.IGNORECASE), "not-simply contrast"),
    (re.compile(r"\bnot\b[^.!?\n]{0,120}\bbut\b", re.IGNORECASE), "not-X-but-Y contrast"),
    (re.compile(r"\bdelve\b", re.IGNORECASE), "AI-tell word: delve"),
    (re.compile(r"\btapestry\b", re.IGNORECASE), "AI-tell word: tapestry"),
    (re.compile(r"\bseamless(?:ly)?\b", re.IGNORECASE), "AI-tell word: seamless"),
    (re.compile(r"\bIt is important to note\b", re.IGNORECASE), "boilerplate intro"),
    (re.compile(r"\bAt its core\b", re.IGNORECASE), "boilerplate intro"),
    (re.compile(r"\b(?:In essence|Simply put|Put simply)\b", re.IGNORECASE), "boilerplate intro"),
    (re.compile(r"\bIt(?: is|'s) worth noting\b", re.IGNORECASE), "boilerplate intro"),
    (
        re.compile(
            r"(?:^|[.!?]\s+)(?:Moreover|Furthermore|Crucially|Importantly|Notably|"
            r"In conclusion|Ultimately|That said|With this in mind|Building on this idea),",
            re.IGNORECASE | re.MULTILINE,
        ),
        "formulaic sentence opener",
    ),
    (re.compile(r"\bthis (?:paper|section|chapter) (?:aims|seeks) to\b", re.IGNORECASE), "anthropomorphized document"),
    (re.compile(r"\bunpack\b", re.IGNORECASE), "AI-tell word: unpack"),
]

READER_IDENTIFIER_PATTERNS = [
    (re.compile(r"\bD\d+[a-z]?\b"), "internal D-lane identifier in reader prose"),
    (re.compile(r"\bUD\d+\b"), "internal UD identifier in reader prose"),
    (re.compile(r"\bRP-A\d+\b"), "internal RP identifier in reader prose"),
    (re.compile(r"\b(?:GAP|CL|DK)-[A-Z0-9-]+\b"), "internal tracker identifier in reader prose"),
    (re.compile(r"\bCP-\d+\b"), "internal premise identifier in reader prose"),
    (re.compile(r"\bHIERARCHY-SCREEN-READOUT\b"), "internal theorem label in reader prose"),
]

ABSTRACT_IDENTIFIER_PATTERNS = [
    (re.compile(r"\bD\d+\b"), "internal D-lane identifier in abstract"),
    (re.compile(r"\bOPH\d+\b"), "internal OPH checklist identifier in abstract"),
    (re.compile(r"\bsigma_ref\b", re.IGNORECASE), "internal code identifier in abstract"),
    (re.compile(r"\bsigma_[A-Za-z0-9_]+\b"), "internal sigma identifier in abstract"),
    (re.compile(r"\bmanifest\b", re.IGNORECASE), "internal manifest reference in abstract"),
    (re.compile(r"\btheorem_contract\b", re.IGNORECASE), "internal code identifier in abstract"),
    (re.compile(r"\bQFT-Q\d+\b", re.IGNORECASE), "internal QFT-stage identifier in abstract"),
    (re.compile(r"\bMGNS-1\b", re.IGNORECASE), "internal modular-state identifier in abstract"),
]

INFORMAL_IDENTIFIER_PATTERNS = [
    (re.compile(r"\bQFT-Q\d+\b", re.IGNORECASE), "internal QFT-stage identifier in informal prose"),
    (re.compile(r"\bMGNS-1\b", re.IGNORECASE), "internal modular-state identifier in informal prose"),
    (re.compile(r"\bFiniteCapBWCertificate\b"), "internal certificate identifier in informal prose"),
    (
        re.compile(r"\bCONDITIONAL_[A-Z0-9_]+\b"),
        "machine status identifier in informal prose",
    ),
]


def iter_paths(globs: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in globs:
        paths.extend(ROOT.glob(pattern))
    return sorted({path for path in paths if path.is_file()})


def line_col(text: str, index: int) -> tuple[int, int]:
    line = text.count("\n", 0, index) + 1
    line_start = text.rfind("\n", 0, index) + 1
    return line, index - line_start + 1


def add_matches(issues: list[str], path: Path, text: str, patterns: list[tuple[re.Pattern[str], str]]) -> None:
    rel = path.relative_to(ROOT)
    for pattern, label in patterns:
        for match in pattern.finditer(text):
            line, col = line_col(text, match.start())
            issues.append(f"{rel}:{line}:{col}: {label}")


def abstracts(text: str) -> list[tuple[int, str]]:
    out: list[tuple[int, str]] = []
    pattern = re.compile(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", re.DOTALL)
    for match in pattern.finditer(text):
        out.append((match.start(1), match.group(1)))
    return out


def main() -> int:
    issues: list[str] = []

    all_paths = iter_paths(READER_GLOBS + STATUS_GLOBS + PAPER_GLOBS)
    for path in all_paths:
        text = path.read_text(encoding="utf-8", errors="ignore")
        add_matches(
            issues,
            path,
            text,
            [(re.compile(r"\b" + "hon" + "est" + r"\b", re.IGNORECASE), "banned h-word")],
        )
        add_matches(issues, path, text, PROGRESS_PATTERNS)
        add_matches(issues, path, text, AI_TELL_PATTERNS)

    for path in iter_paths(READER_GLOBS):
        text = path.read_text(encoding="utf-8", errors="ignore")
        add_matches(issues, path, text, READER_IDENTIFIER_PATTERNS)

    for path in iter_paths(INFORMAL_GLOBS):
        text = path.read_text(encoding="utf-8", errors="ignore")
        add_matches(issues, path, text, INFORMAL_IDENTIFIER_PATTERNS)

    for path in iter_paths(PAPER_GLOBS):
        text = path.read_text(encoding="utf-8", errors="ignore")
        for offset, abstract in abstracts(text):
            for pattern, label in ABSTRACT_IDENTIFIER_PATTERNS:
                for match in pattern.finditer(abstract):
                    line, col = line_col(text, offset + match.start())
                    issues.append(f"{path.relative_to(ROOT)}:{line}:{col}: {label}")
            mar = re.search(r"\bMAR\b", abstract)
            if mar and "Minimal Admissible Realization" not in abstract[: mar.start()]:
                line, col = line_col(text, offset + mar.start())
                issues.append(f"{path.relative_to(ROOT)}:{line}:{col}: MAR used before definition in abstract")

    if issues:
        print("reader style check failed:")
        for issue in issues:
            print(f"- {issue}")
        return 1

    print("reader style check OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
