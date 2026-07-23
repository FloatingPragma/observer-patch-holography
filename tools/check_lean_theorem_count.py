#!/usr/bin/env python3
"""Check the stated Lean theorem-count floor against the library.

The public prose states the size of the Lean library as a lower bound
("more than 700 theorems and lemmas" / "plus de 700 théorèmes et lemmes")
so that the sentence stays true as proofs land. This check counts the
public ``theorem``/``lemma`` declarations under ``Lean/`` and verifies
every stated floor:

  * a floor above the actual count overstates the library and fails hard;
  * a floor more than one round hundred behind the count is stale and
    fails with the expected value.

The expected floor is the largest multiple of 100 strictly below the
count, so "more than N" is always literally true.

Usage:
  python3 tools/check_lean_theorem_count.py          # check, exit 1 on drift
  python3 tools/check_lean_theorem_count.py --fix    # rewrite stale floors
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LEAN_ROOT = REPO_ROOT / "Lean"

# Prose surfaces that state the floor.
CLAIM_FILES = [
    REPO_ROOT / "README.md",
    REPO_ROOT / "README_FR.md",
    REPO_ROOT / "docs" / "HACKER_NEWS_START_HERE.md",
    REPO_ROOT / "extra" / "compact_proof_of_oph.tex",
]

# Public declarations only: a leading `private` keeps the line from matching.
DECL_RE = re.compile(r"^\s*(?:theorem|lemma)\s", re.MULTILINE)

# "more than 700 ... theorems" / "plus de 700 ... théorèmes"; the noun must
# follow within a few words so unrelated numbers never match.
FLOOR_RE = re.compile(
    r"(?:[Mm]ore than|[Pp]lus de)\s+(\d+)(?=[^.\n]{0,80}(?:theorem|théorème))"
)


def count_declarations() -> int:
    total = 0
    for path in sorted(LEAN_ROOT.rglob("*.lean")):
        if any(part.startswith(".") or part == "lake-packages" for part in path.parts):
            continue
        total += len(DECL_RE.findall(path.read_text(encoding="utf-8")))
    return total


def expected_floor(count: int) -> int:
    return max(0, (count - 1) // 100 * 100)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fix", action="store_true", help="rewrite stale floors to the expected value"
    )
    args = parser.parse_args()

    count = count_declarations()
    target = expected_floor(count)
    print(f"Lean public theorem/lemma declarations: {count} (expected floor: {target})")

    problems: list[str] = []
    for path in CLAIM_FILES:
        text = path.read_text(encoding="utf-8")
        matches = list(FLOOR_RE.finditer(text))
        if not matches:
            problems.append(f"{path.relative_to(REPO_ROOT)}: no theorem-count floor found")
            continue
        mismatch = False
        for m in matches:
            floor = int(m.group(1))
            line = text.count("\n", 0, m.start()) + 1
            where = f"{path.relative_to(REPO_ROOT)}:{line}"
            if floor >= count:
                mismatch = True
                problems.append(
                    f"{where}: floor {floor} overstates the library (count {count})"
                )
            elif floor != target:
                mismatch = True
                problems.append(
                    f"{where}: floor {floor} is stale (expected {target}, count {count})"
                )
        if args.fix and mismatch:
            fixed = FLOOR_RE.sub(
                lambda m: m.group(0).replace(m.group(1), str(target)), text
            )
            path.write_text(fixed, encoding="utf-8")
            print(f"fixed: {path.relative_to(REPO_ROOT)}")

    if problems:
        for problem in problems:
            print(f"FAIL {problem}", file=sys.stderr)
        if args.fix:
            print("floors rewritten; re-run without --fix to confirm", file=sys.stderr)
        else:
            print(
                "run `python3 tools/check_lean_theorem_count.py --fix` to rewrite stale floors",
                file=sys.stderr,
            )
        return 1
    print("OK: every stated floor matches the expected value")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
