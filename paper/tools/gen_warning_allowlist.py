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
  1. Build all source-derived papers and the book with the canonical builders.
  2. Run this script with no positional arguments. It obtains the paper set
     directly from ``build_tex_papers.ALL_PAPERS`` and adds the canonical
     reader-facing book log:

       python3 paper/tools/gen_warning_allowlist.py

     Positional files or directories remain available for explicit runs, but
     their expanded log set must equal the registered set exactly. Missing,
     duplicate, and unregistered logs are rejected before the allowlist is
     written. This prevents a hand-typed directory list from silently dropping
     the cosmology or flagship logs and shrinking the warning budget.

The canonical warning budget therefore includes the flagship, core,
supplemental, cosmology, and reader-facing book logs. Compare the regenerated
file against the committed one before staging it; a content change should now
mean warning-content drift, not accidental registry-membership drift.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
from check_build_warnings import parse_log  # noqa: E402

REPO_ROOT = TOOLS.parents[1]
ROOT_TOOLS = REPO_ROOT / "tools"
sys.path.insert(0, str(ROOT_TOOLS))
import build_tex_papers  # noqa: E402

ALLOWLIST = TOOLS.parent / "build_warning_allowlist.json"
BOOK_LOG = REPO_ROOT.parent / "temp" / "book_pdf_build" / "book_manuscript.log"


def registered_logs() -> dict[Path, str]:
    """Return the exact publication-log registry as canonical path -> label.

    Paper membership comes from the same ``ALL_PAPERS`` mapping consumed by
    the publication build/warning gate. The book is registered separately
    because its generated TeX root is owned by ``build_book_pdf.py`` rather
    than ``ALL_PAPERS``.
    """

    logs: dict[Path, str] = {}
    for paper_id, tex_path in build_tex_papers.ALL_PAPERS.items():
        log_path = tex_path.with_suffix(".log").resolve()
        if log_path in logs:
            raise RuntimeError(
                "duplicate registered paper log path: "
                f"{log_path} ({logs[log_path]}, paper:{paper_id})"
            )
        logs[log_path] = f"paper:{paper_id}"

    book_path = BOOK_LOG.resolve()
    if book_path in logs:
        raise RuntimeError(f"book log collides with a registered paper log: {book_path}")
    logs[book_path] = "book:reader-facing"

    # The checker keys allowlist entries by basename. A basename collision
    # would make two distinct publication roots indistinguishable there.
    basenames = Counter(path.name for path in logs)
    collisions = sorted(name for name, count in basenames.items() if count > 1)
    if collisions:
        raise RuntimeError(
            "registered log basenames are not unique: " + ", ".join(collisions)
        )
    return logs


def expand_log_inputs(inputs: list[Path]) -> list[Path]:
    """Expand explicit file/directory inputs without deduplicating them."""

    logs: list[Path] = []
    for candidate in inputs:
        if candidate.is_dir():
            logs.extend(sorted(candidate.glob("*.log")))
        else:
            logs.append(candidate)
    return logs


def validate_exact_coverage(
    log_paths: list[Path], expected: dict[Path, str]
) -> list[Path]:
    """Validate exact registered coverage and return canonical ordered paths.

    Validation is path-based, not basename-based: a copied log with the name
    of a registered paper is not the registered build artifact. Duplicate
    inputs are checked before converting to a set so overlapping directory and
    file arguments cannot be silently absorbed.
    """

    canonical = [path.resolve() for path in log_paths]
    counts = Counter(canonical)
    supplied = set(canonical)
    expected_paths = set(expected)

    duplicate = sorted(path for path, count in counts.items() if count > 1)
    missing = sorted(expected_paths - supplied)
    unregistered = sorted(supplied - expected_paths)
    not_files = sorted(path for path in supplied & expected_paths if not path.is_file())

    problems: list[str] = []
    if duplicate:
        problems.append(
            "duplicate logs:\n  "
            + "\n  ".join(f"{path} ({counts[path]} copies)" for path in duplicate)
        )
    if missing:
        problems.append(
            "missing registered logs:\n  "
            + "\n  ".join(f"{path} [{expected[path]}]" for path in missing)
        )
    if unregistered:
        problems.append(
            "unregistered logs:\n  " + "\n  ".join(str(path) for path in unregistered)
        )
    if not_files:
        problems.append(
            "registered logs not found on disk (build with --keep-logs first):\n  "
            + "\n  ".join(str(path) for path in not_files)
        )
    if problems:
        raise ValueError("\n".join(problems))

    return sorted(expected_paths)


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
                "the warning itself carries no semantic status. Anchored to source line-range + "
                "badness band so a new/unseen "
                "underfull (different location or badness) is reported, not absorbed.",
            }
        )
    return {
        "_comment": "Underfull-box warning budget for issue #542 across all source-derived "
        "paper roots and the reader-facing book. Each entry anchors ONE known "
        "underfull to its stable source line-range and observed badness band (not a per-file blanket "
        "count). A new underfull at a different location or badness is UNEXPLAINED and fails the gate. "
        "Regenerate with paper/tools/gen_warning_allowlist.py after a clean tectonic build "
        "(microtype + emergencystretch). Overfull boxes are NOT budgeted (the checker fails on any overfull).",
        "allow": entries,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate the anchored underfull-warning allowlist from the exact "
            "registered paper and reader-facing book log set."
        )
    )
    parser.add_argument(
        "logs",
        nargs="*",
        type=Path,
        help=(
            "Optional explicit .log files or directories. Their expanded set "
            "must exactly match the registry; default: use registered paths."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ALLOWLIST,
        help=f"allowlist output path (default: {ALLOWLIST})",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    expected = registered_logs()
    supplied = expand_log_inputs(args.logs) if args.logs else list(expected)
    try:
        logs = validate_exact_coverage(supplied, expected)
    except ValueError as exc:
        print(f"gen_warning_allowlist: {exc}", file=sys.stderr)
        return 2

    doc = build(logs)
    args.output.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    boxes = sum(e["max_count"] for e in doc["allow"])
    print(
        f"wrote {len(doc['allow'])} anchored entries covering {boxes} "
        f"underfull boxes from {len(logs)} registered logs -> {args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
