"""Refresh the OpenWave rows of the replay app from the research team's current reports.

Only results on the three frozen candidate laws are taken over.

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
KEEP_WITHOUT_REPORT = {30}   # supplied directly by the owner with its own report


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


# Only results on the three frozen candidate laws are documented. A sentence or table counts
# when it names those laws; sections about older, earlier, separate or historical sources are
# skipped, and so is every sentence that does not name the candidate laws.
CANDIDATE = re.compile(r"finalist|maturation law|maturing|three (?:repair |formation |candidate |maturation )?laws|"
                       r"age-8192|8,192-age|candidate law|minus (?:cosine|linear)|cosine[- ](?:32|128)|linear[- ]128", re.I)
OTHER = re.compile(r"separate|earlier|older|historical|exterior|previously specified|archived|sibling|"
                   r"optical source|colour source|weak-bank|radiative-repair|Dirac", re.I)


def sentences(text: str) -> list[str]:
    return [s for s in re.split(r"(?<=[.!?])\s+", plain(text)) if s]


def candidate_sentences(text: str) -> list[str]:
    return [s for s in sentences(text) if CANDIDATE.search(s) and not OTHER.search(s)]


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
        findings.extend(candidate_sentences(chunk))
    kept = []
    for table in tables(text):
        body = " ".join(" ".join(r) for r in table["rows"])
        if not OTHER.search(table["heading"]) and CANDIDATE.search(body + " " + table["heading"]) and not OTHER.search(body):
            kept.append(table)
    if not findings and not kept:
        return None
    seen, unique = set(), []
    for s in findings:
        if s not in seen:
            seen.add(s); unique.append(s)
    return {"title": title, "findings": unique[:8], "tables": kept}


def status_class(status: str) -> str:
    """pass only for a pass that the status attributes to the candidate laws; open otherwise."""
    if "pass" in status.lower() and CANDIDATE.search(status) and not OTHER.search(status):
        return "pass"
    return "open"


def sync(reports: Path = DEFAULT_REPORTS) -> dict:
    readme = (reports / "README.md").read_text(encoding="utf-8")
    current = {}
    for line in readme.splitlines():
        m = re.match(r"^\|\s*\*{0,2}(\d+)\s*·", line)
        if m:
            cells = _cells(line)
            current[int(m.group(1))] = {"status": cells[1], "result": cells[2]}
    data = json.loads(ROWS.read_text(encoding="utf-8"))
    for row in data["rows"]:
        n = row["number"]
        if n in KEEP_WITHOUT_REPORT and n not in current:
            continue
        found = None
        if n in current:
            path = reports / f"ROW{n:02d}_RESULTS.md"
            found = report(path) if path.exists() else None
            lines = candidate_sentences(current[n]["result"])
            if lines or found:
                row["status"] = status_class(current[n]["status"])
                row["status_label"] = "passed on the candidate laws" if row["status"] == "pass" else "open on the candidate laws"
                row["result"] = " ".join(lines or (found["findings"][:2] if found else []))
                row["report"] = found
                continue
        row.update(status="open", status_label="no result on the candidate laws", result="", report=None,
                   comparison=None, runs=[])
    data["source"]["research_status"] = "current reports of the research team"
    ROWS.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return {n: v["status"] for n, v in current.items()}


if __name__ == "__main__":
    import sys
    print(sync(Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_REPORTS))
