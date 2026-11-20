"""Refresh the OpenWave rows of the replay app from the research team's current reports.

Rows with a passing result on the OPH law are ported by hand into the simulator, with an
exact check and a replay; the sync keeps those and greys every other row.

The reports directory holds a status README with one table row per reported row and one
``ROWnn_RESULTS.md`` per row. The sync takes each reported row's status and latest result
from the table, and the title, headline and tables of its report. Links are dropped, so no
path into the research repository reaches the app. Rows without a report keep no result.
The problem explanations, animations and charts in ``openwave_rows.json`` are kept.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROWS = HERE / "openwave_rows.json"
DEFAULT_REPORTS = HERE.parent.parent / "openwave" / "simulation" / "reports"
# Rows whose passing result is ported into the simulator, with a check and a replay. The sync
# leaves them as they are; every other row is greyed until its result is ported.
PORTED = {4, 30}


def plain(text: str) -> str:
    """Markdown to plain text: drop links with their labels, images and emphasis."""
    text = re.sub(r"!\[[^\]]*\]\([^)]*\)", "", text)
    # Links in the reports cite research files; drop them with their labels.
    text = re.sub(r"(?:\s*·\s*)?\[[^\]]*\]\([^)]*\)", "", text)
    text = re.sub(r"(?:\s*·\s*)+$", "", text.strip())
    text = text.replace("**", "").replace("`", "")
    return re.sub(r"\s+", " ", text).strip()


def _cells(line: str) -> list[str]:
    return [plain(c) for c in line.strip().strip("|").split("|")]


def tables(markdown: str) -> list[dict]:
    out, heading, lines = [], "", markdown.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("#"):
            heading = plain(line.lstrip("#"))
        if line.startswith("|") and i + 1 < len(lines) and re.match(r"^\|\s*:?-{3}", lines[i + 1]):
            header, rows = _cells(line), []
            i += 2
            while i < len(lines) and lines[i].startswith("|"):
                rows.append(_cells(lines[i]))
                i += 1
            out.append({"heading": heading, "header": header, "rows": rows})
            continue
        i += 1
    return out


# The reports call the OPH law and the variants run alongside it the finalists or the maturation
# laws. A sentence or table counts when it names them; sections about older, earlier, separate or historical sources are
# skipped, and so is every sentence that does not name the law.
LAWNAME = re.compile(r"finalist|maturation law|maturing|three (?:repair |formation |candidate |maturation )?laws|"
                       r"age-8192|8,192-age|candidate law|minus (?:cosine|linear)|cosine[- ](?:32|128)|linear[- ]128", re.I)
OTHER = re.compile(r"separate|earlier|older|historical|exterior|previously specified|archived|sibling|"
                   r"optical source|colour source|weak-bank|radiative-repair|Dirac", re.I)


def sentences(text: str) -> list[str]:
    return [s for s in re.split(r"(?<=[.!?])\s+", plain(text)) if s]


def law_sentences(text: str) -> list[str]:
    return [s for s in sentences(text) if LAWNAME.search(s) and not OTHER.search(s)]


def report(path: Path) -> dict | None:
    text = path.read_text(encoding="utf-8")
    title = next((plain(l.lstrip("#")) for l in text.splitlines() if l.startswith("# ")), "")
    blocks, heading = [], ""
    for chunk in text.split("\n\n"):
        chunk = chunk.strip()
        if chunk.startswith("#"):
            heading = plain(chunk.splitlines()[0].lstrip("#"))
            continue
        blocks.append((heading, chunk))
    findings = []
    for heading, chunk in blocks:
        if OTHER.search(heading) or chunk.startswith(("|", "```", "$$", "![")):
            continue
        findings.extend(law_sentences(chunk))
    kept = []
    for table in tables(text):
        body = " ".join(" ".join(r) for r in table["rows"])
        if not OTHER.search(table["heading"]) and LAWNAME.search(body + " " + table["heading"]) and not OTHER.search(body):
            kept.append(table)
    if not findings and not kept:
        return None
    seen, unique = set(), []
    for s in findings:
        if s not in seen:
            seen.add(s); unique.append(s)
    return {"title": title, "findings": unique[:8], "tables": kept}


def status_class(status: str) -> str:
    """pass only for a pass that the status attributes to the law; open otherwise."""
    if "pass" in status.lower() and LAWNAME.search(status) and not OTHER.search(status):
        return "pass"
    return "open"


def reported_rows(readme: str) -> dict:
    """Rows named in the status README: a table row ('| **N · ...') or a list item ('Row N · ...')."""
    found = {}
    for line in readme.splitlines():
        m = re.match(r"^\|\s*\*{0,2}(\d+)\s*·", line) or re.search(r"\[Row (\d+) ·", line)
        if m:
            found[int(m.group(1))] = plain(line)
    return found


def sync(reports: Path = DEFAULT_REPORTS) -> dict:
    """Grey every row that is not ported; return the rows the reports cover, for review."""
    readme = (reports / "README.md").read_text(encoding="utf-8")
    current = reported_rows(readme)
    data = json.loads(ROWS.read_text(encoding="utf-8"))
    for row in data["rows"]:
        if row["number"] in PORTED:
            continue
        row.update(status="open", status_label="no result", result="", report=None, comparison=None, runs=[])
    data["source"]["research_status"] = "current reports of the research team"
    ROWS.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return {n: ("ported" if n in PORTED else "reported, not ported") for n in sorted(current)}


if __name__ == "__main__":
    import sys
    print(sync(Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_REPORTS))
