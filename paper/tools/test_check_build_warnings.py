#!/usr/bin/env python3
"""Regression tests for the underfull warning gate (issue #542).

The point of the anchored allowlist is that it is NOT a per-file blanket count:
an underfull box at a known source location + badness is budgeted, but a NEW /
unseen underfull in the same source (different line-range or badness) must be
reported, not silently absorbed. These tests pin that behavior.

Run with: python -m pytest paper/tools/test_check_build_warnings.py
"""
from __future__ import annotations

import sys
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
import check_build_warnings as chk  # noqa: E402

ALLOWLIST = TOOLS.parent / "build_warning_allowlist.json"


def _match(tmp_path: Path, log_name: str, body: str) -> dict[str, chk.BoxWarning]:
    log = tmp_path / log_name
    log.write_text(body, encoding="utf-8")
    warnings, _ = chk.parse_log(log)
    allowlist = chk.load_allowlist(ALLOWLIST)
    chk.match_allowlist(warnings, allowlist, log.name)
    return {w.lines: w for w in warnings if w.kind.startswith("underfull")}


def _pick_anchored_entry() -> dict:
    """A real allowlist entry to build a 'known' warning from."""
    entries = chk.load_allowlist(ALLOWLIST)
    assert entries, "allowlist is empty"
    return entries[0]


def test_known_underfull_is_allowed(tmp_path: Path) -> None:
    e = _pick_anchored_entry()
    body = (
        f"({e['source_file']}\n"
        f"Underfull \\hbox (badness {e['badness_min']}) in paragraph at lines {e['lines']}\n"
        " []\\T1/lmr/m/n/10 some justified prose that hyphenates\n\n)\n"
    )
    by_lines = _match(tmp_path, e["log"], body)
    w = by_lines[e["lines"]]
    assert w.allowed_by is not None, "a known anchored underfull must stay budgeted"


def test_unseen_underfull_in_allowlisted_source_is_unexplained(tmp_path: Path) -> None:
    # Same source file as a real allowlist entry, but a NOVEL line-range: the
    # anchored budget must not absorb it (this is exactly the regression the
    # blanket per-file count allowed).
    e = _pick_anchored_entry()
    body = (
        f"({e['source_file']}\n"
        "Underfull \\hbox (badness 4200) in paragraph at lines 99991--99999\n"
        " []\\T1/lmr/m/n/10 a brand new paragraph that never existed before\n\n)\n"
    )
    by_lines = _match(tmp_path, e["log"], body)
    w = by_lines["99991--99999"]
    assert w.allowed_by is None, "an unseen underfull at a new location must be unexplained"


def test_same_location_but_worse_badness_is_unexplained(tmp_path: Path) -> None:
    # A finite-band entry: a much worse badness at the same line-range falls
    # outside the observed band and is reported (tighter-bound behavior).
    entries = chk.load_allowlist(ALLOWLIST)
    finite = next((e for e in entries if e["badness_max"] < 10000), None)
    if finite is None:
        return  # no finite-band entry to exercise; skip cleanly
    worse = finite["badness_max"] + 5000
    body = (
        f"({finite['source_file']}\n"
        f"Underfull \\hbox (badness {worse}) in paragraph at lines {finite['lines']}\n"
        " []\\T1/lmr/m/n/10 same place, materially worse box\n\n)\n"
    )
    by_lines = _match(tmp_path, finite["log"], body)
    w = by_lines[finite["lines"]]
    assert w.allowed_by is None, "badness outside the anchored band must be unexplained"
