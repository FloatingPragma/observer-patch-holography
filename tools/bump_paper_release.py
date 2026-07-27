#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from datetime import datetime
from pathlib import Path
import json


RELEASE_INFO_RELATIVE = Path("paper/release_info.tex")
CLAIM_REGISTRY_RELATIVE = Path("claims/claim_registry.yaml")
RELEASE_ID_MACRO = "OPHPaperReleaseID"
RELEASE_DATE_MACRO = "OPHPaperReleaseDate"
RELEASE_NUMBER_PATTERN = re.compile(r"^(?P<prefix>.*?)(?P<number>\d+)$")
NEXT_STEP = (
    "Next: run python3 tools/refresh_paper_release.py --publication"
)


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    release_info_path = repo_root / RELEASE_INFO_RELATIVE
    claim_registry_path = repo_root / CLAIM_REGISTRY_RELATIVE
    text = release_info_path.read_text(encoding="utf-8")
    registry = json.loads(claim_registry_path.read_text(encoding="utf-8"))

    current_release_id = extract_macro(text, RELEASE_ID_MACRO)
    current_release_date = extract_macro(text, RELEASE_DATE_MACRO)
    next_release_id = args.release_id or increment_release_id(current_release_id)
    next_release_date = args.date or today_release_date()

    updated_text = replace_macro(text, RELEASE_ID_MACRO, next_release_id)
    updated_text = replace_macro(updated_text, RELEASE_DATE_MACRO, next_release_date)

    if args.dry_run:
        print(f"{release_info_path}")
        print(f"release_id: {current_release_id} -> {next_release_id}")
        print(f"released_at: {current_release_date} -> {next_release_date}")
        print(f"claim_registry.release_id: {registry.get('release_id')} -> {next_release_id}")
        return 0

    registry["release_id"] = next_release_id
    release_info_path.write_text(updated_text, encoding="utf-8")
    claim_registry_path.write_text(
        json.dumps(registry, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Updated {release_info_path}")
    print(f"Updated {claim_registry_path}")
    print(f"release_id: {current_release_id} -> {next_release_id}")
    print(f"released_at: {current_release_date} -> {next_release_date}")
    print(NEXT_STEP)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Bump the shared paper release ID and release date in paper/release_info.tex for all current papers.",
    )
    parser.add_argument(
        "--release-id",
        help="Set an explicit release ID instead of incrementing the numeric suffix.",
    )
    parser.add_argument(
        "--date",
        help="Set an explicit release date instead of using today's local date.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the proposed release change without writing the file.",
    )
    return parser.parse_args()


def extract_macro(text: str, macro_name: str) -> str:
    pattern = re.compile(r"\\newcommand\{\\%s\}\{([^}]*)\}" % re.escape(macro_name))
    match = pattern.search(text)
    if not match:
        raise SystemExit(f"missing macro {macro_name} in release info")
    return match.group(1).strip()


def replace_macro(text: str, macro_name: str, value: str) -> str:
    pattern = re.compile(r"(\\newcommand\{\\%s\}\{)([^}]*)(\})" % re.escape(macro_name))
    updated_text, count = pattern.subn(r"\g<1>%s\g<3>" % value, text, count=1)
    if count != 1:
        raise SystemExit(f"missing macro {macro_name} in release info")
    return updated_text


def increment_release_id(release_id: str) -> str:
    match = RELEASE_NUMBER_PATTERN.match(release_id)
    if not match:
        raise SystemExit(
            f"could not increment release id {release_id!r}; pass --release-id to set it explicitly"
        )
    prefix = match.group("prefix")
    number = int(match.group("number"))
    return f"{prefix}{number + 1}"


def today_release_date() -> str:
    today = datetime.now().astimezone().date()
    return f"{today.strftime('%B')} {today.day}, {today.year}"


if __name__ == "__main__":
    raise SystemExit(main())
