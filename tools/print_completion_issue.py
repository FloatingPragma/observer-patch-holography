#!/usr/bin/env python3
"""Print the canonical GitHub issue body for one completion gate.

``tools/build_completion_program_current.py`` uses the retained renderer. Its base module is
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
import build_completion_work_queue as work_queue


def issue_body(stage: dict, program: dict, *, include_work_queue: bool = True,
               root: Path = builder.ROOT) -> str:
    """Return the canonical contract, checked predecessors and owned tasks."""
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
    if include_work_queue:
        queue = work_queue.load_queue(root, program)
        section = work_queue.render_section(queue, stage['id'])
        if section:
            anchor = '## Acceptance\n'
            if body.count(anchor) != 1:
                raise ValueError(f"{stage['id']}: task insertion anchor changed")
            body = body.replace(anchor, section + '\n' + anchor, 1)
    return body


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gate", metavar="CP-XX", help="completion gate to print")
    parser.add_argument('--without-work-queue', action='store_true',
                        help='print only the retained gate contract and predecessor checkboxes')
    args = parser.parse_args()
    program, _ = builder.load_program()
    by_id = {row["id"]: row for row in program["stages"]}
    if args.gate not in by_id:
        raise SystemExit(f"unknown completion gate: {args.gate}")
    print(issue_body(by_id[args.gate], program,
                     include_work_queue=not args.without_work_queue), end="")


if __name__ == "__main__":
    main()
