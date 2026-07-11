#!/usr/bin/env python3
"""Generate a reproducible author-contribution report from git history."""

from __future__ import annotations

import json
import subprocess
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TRACKING_DIR = Path(__file__).resolve().parent
CONFIG_PATH = TRACKING_DIR / "authors.json"
LATEST_JSON = TRACKING_DIR / "latest.json"
LATEST_MD = TRACKING_DIR / "latest.md"


@dataclass(frozen=True)
class Alias:
    name: str
    email: str

    @property
    def author_spec(self) -> str:
        return f"{self.name} <{self.email}>"


def git(*args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(ROOT), *args],
        text=True,
        errors="replace",
    )


def load_config() -> dict:
    with CONFIG_PATH.open() as fh:
        return json.load(fh)


def classify_path(path: str) -> str:
    if path.startswith("paper/"):
        return "paper"
    if path.startswith("book/"):
        return "book"
    if path in {"README.md", "README_FR.md"}:
        return "readme"
    if path.startswith("extra/"):
        return "extra"
    if path.startswith("code/"):
        return "code"
    if path.startswith("assets/"):
        return "assets"
    if path.lower().startswith("lean/"):
        return "lean"
    return "other"


def list_commits(alias: Alias) -> dict[str, dict]:
    output = git(
        "log",
        "--all",
        "--no-merges",
        "--author",
        alias.author_spec,
        "--date=short",
        "--format=%H\t%ad\t%s",
    )
    commits = {}
    for line in output.splitlines():
        if not line.strip():
            continue
        commit_hash, date, subject = line.split("\t", 2)
        commits[commit_hash] = {
            "hash": commit_hash,
            "date": date,
            "subject": subject,
        }
    return commits


def commit_numstat(commit_hash: str) -> list[tuple[str, str, str]]:
    rows = []
    output = git("show", "--format=", "--numstat", commit_hash)
    for line in output.splitlines():
        if "\t" not in line:
            continue
        added, deleted, path = line.split("\t", 2)
        rows.append((added, deleted, path))
    return rows


def build_report() -> dict:
    config = load_config()
    paper_targets = config["paper_targets"]
    report = {
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "repo_head": git("rev-parse", "--short", "HEAD").strip(),
        "repo_status_short": git("status", "--short").splitlines(),
        "authors": [],
        "paper_rankings": [],
    }

    per_author = []

    for author in config["authors"]:
        aliases = [Alias(**alias) for alias in author["aliases"]]
        commit_map: dict[str, dict] = {}
        for alias in aliases:
            commit_map.update(list_commits(alias))

        category_counts = Counter()
        file_touch_counts = Counter()
        paper_commit_counts = Counter()
        total_adds = 0
        total_dels = 0

        for commit_hash, commit in commit_map.items():
            touched_paths = set()
            for added, deleted, path in commit_numstat(commit_hash):
                touched_paths.add(path)
                file_touch_counts[path] += 1
                category_counts[classify_path(path)] += 1
                if added.isdigit() and deleted.isdigit():
                    total_adds += int(added)
                    total_dels += int(deleted)
            for target in paper_targets:
                if any(path in touched_paths for path in target["paths"]):
                    paper_commit_counts[target["key"]] += 1

        paper_commits = [
            commit
            for commit_hash, commit in commit_map.items()
            if any(path.startswith("paper/") for _, _, path in commit_numstat(commit_hash))
        ]
        paper_commits.sort(key=lambda item: (item["date"], item["hash"]), reverse=True)

        author_entry = {
            "display_name": author["display_name"],
            "notes": author.get("notes", ""),
            "aliases": [alias.author_spec for alias in aliases],
            "non_merge_commit_count": len(commit_map),
            "line_additions": total_adds,
            "line_deletions": total_dels,
            "touched_file_count": sum(file_touch_counts.values()),
            "unique_file_count": len(file_touch_counts),
            "category_touch_counts": dict(sorted(category_counts.items())),
            "top_files": [
                {"path": path, "touches": touches}
                for path, touches in file_touch_counts.most_common(15)
            ],
            "paper_commit_counts": {
                target["key"]: paper_commit_counts.get(target["key"], 0)
                for target in paper_targets
            },
            "recent_paper_commits": paper_commits[:15],
        }
        per_author.append(author_entry)

    report["authors"] = per_author

    for target in paper_targets:
        ranking = []
        for author in per_author:
            count = author["paper_commit_counts"][target["key"]]
            if count:
                ranking.append(
                    {
                        "display_name": author["display_name"],
                        "commit_count": count,
                    }
                )
        ranking.sort(key=lambda item: (-item["commit_count"], item["display_name"]))
        report["paper_rankings"].append(
            {
                "key": target["key"],
                "label": target["label"],
                "paths": target["paths"],
                "ranking": ranking,
            }
        )

    return report


def render_markdown(report: dict) -> str:
    lines = []
    lines.append("# Author Contribution Tracking")
    lines.append("")
    lines.append(f"- Generated (UTC): `{report['generated_at_utc']}`")
    lines.append(f"- Repo head: `{report['repo_head']}`")
    if report["repo_status_short"]:
        lines.append("- Repo status when generated: dirty")
    else:
        lines.append("- Repo status when generated: clean")
    lines.append("")
    lines.append("## Author Summary")
    lines.append("")
    lines.append("| Author | Non-merge commits | +lines | -lines | Category touches |")
    lines.append("| --- | ---: | ---: | ---: | --- |")
    for author in report["authors"]:
        categories = ", ".join(
            f"{name}={count}"
            for name, count in sorted(author["category_touch_counts"].items())
        )
        lines.append(
            f"| {author['display_name']} | {author['non_merge_commit_count']} | "
            f"{author['line_additions']} | {author['line_deletions']} | {categories} |"
        )
    lines.append("")
    lines.append("## Paper Rankings")
    lines.append("")
    for ranking in report["paper_rankings"]:
        lines.append(f"### {ranking['label']}")
        lines.append("")
        if not ranking["ranking"]:
            lines.append("No tracked co-author commits on this file.")
            lines.append("")
            continue
        lines.append("| Rank | Author | Commits touching target file |")
        lines.append("| --- | --- | ---: |")
        for index, row in enumerate(ranking["ranking"], start=1):
            lines.append(f"| {index} | {row['display_name']} | {row['commit_count']} |")
        lines.append("")

    lines.append("## Detail")
    lines.append("")
    for author in report["authors"]:
        lines.append(f"### {author['display_name']}")
        lines.append("")
        lines.append(f"- Aliases: `{', '.join(author['aliases'])}`")
        if author["notes"]:
            lines.append(f"- Notes: {author['notes']}")
        lines.append(f"- Non-merge commits: `{author['non_merge_commit_count']}`")
        lines.append(
            f"- Line delta: `+{author['line_additions']} / -{author['line_deletions']}`"
        )
        lines.append(
            f"- Files touched: `{author['touched_file_count']}` touches across `{author['unique_file_count']}` unique files"
        )
        lines.append("- Core paper file counts:")
        for key, count in author["paper_commit_counts"].items():
            if count:
                lines.append(f"  - `{key}`: `{count}`")
        lines.append("- Most-touched files:")
        for file_row in author["top_files"][:10]:
            lines.append(f"  - `{file_row['path']}`: `{file_row['touches']}`")
        lines.append("- Recent paper commits:")
        for commit in author["recent_paper_commits"][:10]:
            lines.append(
                f"  - `{commit['date']}` `{commit['hash'][:7]}` {commit['subject']}"
            )
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    report = build_report()
    LATEST_JSON.write_text(json.dumps(report, indent=2) + "\n")
    LATEST_MD.write_text(render_markdown(report) + "\n")
    print(f"Wrote {LATEST_JSON}")
    print(f"Wrote {LATEST_MD}")


if __name__ == "__main__":
    main()
