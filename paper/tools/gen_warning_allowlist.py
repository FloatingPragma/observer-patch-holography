#!/usr/bin/env python3
"""Regenerate paper/build_warning_allowlist.json from a clean build (issue #542).

Each underfull box is anchored to its stable source line-range and observed
badness band, so the allowlist is a per-location budget, not a per-file blanket
count: a new/unseen underfull (different location, or badness outside the band)
is reported by check_build_warnings.py, not silently absorbed.

Excerpt text is deliberately NOT used as the anchor: TeX discretionary hyphens
shift with engine/font drift and some boxes have empty excerpts, whereas the
source line-range is stable across rebuilds. Overfull boxes are never budgeted.

Usage:
  1. Build each release root with logs, e.g.:
       tectonic -X compile paper/<root>.tex --keep-logs -o <logdir>
  2. python3 paper/tools/gen_warning_allowlist.py <logdir> [<logdir> ...]

The logs whose basenames match the existing allowlist's ``log`` fields are used;
pass the directory (or directories) containing the freshly built .log files.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
from check_build_warnings import parse_log  # noqa: E402

ALLOWLIST = TOOLS.parent / "build_warning_allowlist.json"


def build(log_paths: list[Path]) -> dict:
    groups: dict[tuple, list[int]] = defaultdict(list)
    for log in log_paths:
        warnings, _ = parse_log(log)
        for w in warnings:
            if not w.kind.startswith("underfull"):
                continue
            groups[(log.name, w.source_file, w.kind, w.lines)].append(w.badness or 0)

    entries = []
    for (log_name, src, kind, lines), badnesses in sorted(groups.items()):
        entries.append(
            {
                "id": f"underfull::{log_name[:-4]}::{src}::L{lines}",
                "log": log_name,
                "source_file": src,
                "kind": kind,
                "lines": lines,
                "badness_min": min(badnesses),
                "badness_max": max(badnesses),
                "max_count": len(badnesses),
                "reason": "Layout-only microtype/hyphenation underfull at a fixed source location; "
                "no claim change. Anchored to source line-range + badness band so a new/unseen "
                "underfull (different location or badness) is reported, not absorbed.",
            }
        )
    return {
        "_comment": "Underfull-box warning budget for issue #542. Each entry anchors ONE known "
        "underfull to its stable source line-range and observed badness band (not a per-file blanket "
        "count). A new underfull at a different location or badness is UNEXPLAINED and fails the gate. "
        "Regenerate with paper/tools/gen_warning_allowlist.py after a clean tectonic build "
        "(microtype + emergencystretch). Overfull boxes are NOT budgeted (the checker fails on any overfull).",
        "allow": entries,
    }


def main() -> int:
    dirs = [Path(a) for a in sys.argv[1:]] or [Path.cwd()]
    logs: list[Path] = []
    for d in dirs:
        logs.extend(sorted(d.glob("*.log")) if d.is_dir() else [d])
    if not logs:
        print("no .log files found; build the roots with tectonic --keep-logs first", file=sys.stderr)
        return 2
    doc = build(logs)
    ALLOWLIST.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    boxes = sum(e["max_count"] for e in doc["allow"])
    print(f"wrote {len(doc['allow'])} anchored entries covering {boxes} underfull boxes -> {ALLOWLIST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
