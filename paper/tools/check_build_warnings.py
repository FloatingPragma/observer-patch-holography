#!/usr/bin/env python3
"""Release check for paper build logs (issue #542).

Parses a LaTeX/tectonic .log file and extracts:
  - underfull \\hbox / \\vbox warnings (with source file, line range, badness),
  - overfull \\hbox / \\vbox warnings,
  - undefined references, undefined citations, multiply defined labels.

Underfull warnings are matched against paper/build_warning_allowlist.json.
The script exits nonzero when any warning falls outside the allowlist, when
any overfull box is present, or when any reference/citation/label problem is
present. The allowlist applies to underfull boxes only.

Usage (from repo root or anywhere):
  python3 paper/tools/check_build_warnings.py paper/observers_are_all_you_need.log \\
      paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.log
  python3 paper/tools/check_build_warnings.py --list <log>   # dump every box warning with source mapping

Build the logs first with: tectonic -X compile <root>.tex --keep-logs
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

PAPER_DIR = Path(__file__).resolve().parent.parent
DEFAULT_ALLOWLIST = PAPER_DIR / "build_warning_allowlist.json"

UNDERFULL_RE = re.compile(
    r"^Underfull \\([hv]box) \(badness (\d+)\) "
    r"(?:in paragraph at lines (\d+)--(\d+)|in alignment at lines (\d+)--(\d+)|detected at line (\d+)|has occurred while \\output is active)"
)
OVERFULL_RE = re.compile(
    r"^Overfull \\([hv]box) \(([\d.]+)pt too (?:wide|high)\) "
    r"(?:in paragraph at lines (\d+)--(\d+)|in alignment at lines (\d+)--(\d+)|detected at line (\d+)|has occurred while \\output is active)"
)
REF_WARNING_RES = (
    re.compile(r"LaTeX Warning: (Reference|Citation) `([^']+)' .*undefined"),
    re.compile(r"LaTeX Warning: Label `([^']+)' multiply defined"),
    re.compile(r"LaTeX Warning: There were undefined (references|citations)"),
    re.compile(r"LaTeX Warning: There were multiply-defined labels"),
)
FONT_SPEC_RE = re.compile(r"\\(?:T1|OT1|OML|OMS|OMX|U)/[\w-]+/[\w]+/[\w]+/[\d.]+ ?")


@dataclass
class BoxWarning:
    kind: str  # underfull_hbox, overfull_vbox, ...
    source_file: str
    lines: str
    badness: int | None
    overwidth: float | None
    excerpt: str
    log_line: int
    pages_before: int | None  # last shipped page number seen before the warning
    allowed_by: str | None = None

    def location(self) -> str:
        page = f", after page {self.pages_before}" if self.pages_before else ""
        return f"{self.source_file}:{self.lines}{page}"


def unwrap_log(text: str) -> list[str]:
    """Join lines hard-wrapped at the classic 79-column log width."""
    out: list[str] = []
    buf = ""
    for raw in text.splitlines():
        buf += raw
        if len(raw) == 79:
            continue
        out.append(buf)
        buf = ""
    if buf:
        out.append(buf)
    return out


def normalize_excerpt(excerpt: str) -> str:
    """Strip font-spec noise and collapse whitespace for stable matching."""
    text = FONT_SPEC_RE.sub("", excerpt)
    text = text.replace("[]", " ")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


class FileTracker:
    """Track the TeX file-open parenthesis stack across a log stream."""

    FILE_OPEN_RE = re.compile(r"\(((?:\./|\.\./|/|[\w.-]+/)?[\w./\\-]+\.(?:tex|sty|cls|aux|toc|def|fd|cfg|clo|bbl|out|ltx|dfu))")
    PAGE_RE = re.compile(r"\[(\d+)[\]{ ]")

    def __init__(self) -> None:
        self.stack: list[str] = []
        self.last_page: int | None = None

    def feed(self, line: str) -> None:
        i = 0
        for match in self.PAGE_RE.finditer(line):
            self.last_page = int(match.group(1))
        while i < len(line):
            ch = line[i]
            if ch == "(":
                m = self.FILE_OPEN_RE.match(line, i)
                if m:
                    self.stack.append(m.group(1))
                    i = m.end()
                    continue
                self.stack.append("<group>")
            elif ch == ")":
                if self.stack:
                    self.stack.pop()
            i += 1

    def current_tex(self) -> str:
        for name in reversed(self.stack):
            if name.endswith(".tex"):
                return name
        return "<root>"


def parse_log(log_path: Path) -> tuple[list[BoxWarning], list[str]]:
    lines = unwrap_log(log_path.read_text(errors="replace"))
    tracker = FileTracker()
    warnings: list[BoxWarning] = []
    ref_problems: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        m_under = UNDERFULL_RE.match(line)
        m_over = OVERFULL_RE.match(line)
        if m_under or m_over:
            m = m_under or m_over
            box = m.group(1)
            nums = [g for g in m.groups()[2:] if g]
            line_range = "--".join(nums) if nums else "output"
            excerpt_lines: list[str] = []
            j = i + 1
            while j < len(lines) and lines[j].strip() not in ("", "[]") and not lines[j].startswith(("Underfull", "Overfull", "LaTeX", "Package", "(", ")")):
                excerpt_lines.append(lines[j])
                j += 1
            excerpt = normalize_excerpt(" ".join(excerpt_lines))
            warnings.append(
                BoxWarning(
                    kind=("underfull_" if m_under else "overfull_") + box,
                    source_file=tracker.current_tex(),
                    lines=line_range,
                    badness=int(m.group(2)) if m_under else None,
                    overwidth=float(m.group(2)) if m_over else None,
                    excerpt=excerpt,
                    log_line=i + 1,
                    pages_before=tracker.last_page,
                )
            )
            i = j
            continue
        for ref_re in REF_WARNING_RES:
            if ref_re.search(line):
                ref_problems.append(line.strip())
                break
        tracker.feed(line)
        i += 1
    return warnings, ref_problems


def load_allowlist(path: Path) -> list[dict]:
    if not path.exists():
        return []
    data = json.loads(path.read_text())
    return data.get("allow", [])


def match_allowlist(warnings: list[BoxWarning], allowlist: list[dict], log_name: str) -> None:
    counts: dict[str, int] = {}
    for w in warnings:
        if not w.kind.startswith("underfull"):
            continue
        for entry in allowlist:
            if entry.get("log") and entry["log"] not in log_name:
                continue
            if entry.get("kind") and entry["kind"] != w.kind:
                continue
            if entry.get("source_file") and entry["source_file"] != w.source_file:
                continue
            if entry.get("lines") and entry["lines"] != w.lines:
                continue
            if entry.get("excerpt_contains") and entry["excerpt_contains"] not in w.excerpt:
                continue
            if entry.get("badness_min") and (w.badness or 0) < entry["badness_min"]:
                continue
            if entry.get("badness_max") and (w.badness or 0) > entry["badness_max"]:
                continue
            eid = entry["id"]
            if counts.get(eid, 0) >= entry.get("max_count", 1):
                continue
            counts[eid] = counts.get(eid, 0) + 1
            w.allowed_by = eid
            break


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("logs", nargs="+", type=Path, help=".log files to check")
    parser.add_argument("--allowlist", type=Path, default=DEFAULT_ALLOWLIST)
    parser.add_argument("--list", action="store_true", help="dump every box warning with its source mapping and exit 0")
    args = parser.parse_args()

    allowlist = load_allowlist(args.allowlist)
    failures = 0
    for log_path in args.logs:
        if not log_path.exists():
            print(f"FAIL {log_path}: log file missing (build with tectonic -X compile <root>.tex --keep-logs)")
            failures += 1
            continue
        warnings, ref_problems = parse_log(log_path)
        underfull = [w for w in warnings if w.kind.startswith("underfull")]
        overfull = [w for w in warnings if w.kind.startswith("overfull")]
        match_allowlist(warnings, allowlist, log_path.name)
        unexplained = [w for w in underfull if w.allowed_by is None]

        if args.list:
            print(f"== {log_path.name}: {len(underfull)} underfull, {len(overfull)} overfull, {len(ref_problems)} ref/cite/label problems")
            for w in warnings:
                tag = f" [allowed: {w.allowed_by}]" if w.allowed_by else ""
                extra = f"badness {w.badness}" if w.badness is not None else f"{w.overwidth}pt"
                print(f"  {w.kind:16s} {w.location():60s} {extra:14s}{tag}  | {w.excerpt[:90]}")
            for p in ref_problems:
                print(f"  REF: {p}")
            continue

        ok = not unexplained and not overfull and not ref_problems
        status = "OK  " if ok else "FAIL"
        print(f"{status} {log_path.name}: {len(underfull)} underfull ({len(unexplained)} unexplained), {len(overfull)} overfull, {len(ref_problems)} ref/cite/label problems")
        for w in unexplained:
            print(f"  UNEXPLAINED {w.kind} at {w.location()} (badness {w.badness}): {w.excerpt[:120]}")
        for w in overfull:
            print(f"  OVERFULL {w.kind} at {w.location()} ({w.overwidth}pt): {w.excerpt[:120]}")
        for p in ref_problems:
            print(f"  REF: {p}")
        if not ok:
            failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
