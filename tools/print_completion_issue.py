#!/usr/bin/env python3
"""Print the canonical GitHub issue body for one completion gate.

``tools/build_completion_program.py`` renders the same body. That module is
pinned as a source by replay receipts, so this wrapper corrects the dependency
checklist at print time: a predecessor whose stage is ``ACCEPTED`` is ticked.
It never writes a file.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_completion_program_current as builder


def issue_body(stage: dict, program: dict) -> str:
    """Return the generated issue body with accepted predecessors ticked."""
    by_id = {row["id"]: row for row in program["stages"]}
    body = builder.render_issue(stage, program)
    for dependency in stage["requires"]:
        predecessor = by_id[dependency]
        if predecessor["status"] != "ACCEPTED":
            continue
        link = builder.issue_link(predecessor["issue"])
        unticked = f"- [ ] {link}: {predecessor['title']}"
        ticked = f"- [x] {link}: {predecessor['title']}"
        if unticked not in body:
            raise ValueError(f"{stage['id']}: missing dependency line for {dependency}")
        body = body.replace(unticked, ticked)
    return body


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gate", metavar="CP-XX", help="completion gate to print")
    args = parser.parse_args()
    program, _ = builder.load_program()
    by_id = {row["id"]: row for row in program["stages"]}
    if args.gate not in by_id:
        raise SystemExit(f"unknown completion gate: {args.gate}")
    print(issue_body(by_id[args.gate], program), end="")


if __name__ == "__main__":
    main()
