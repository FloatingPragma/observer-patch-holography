"""Shared deterministic environment for OPH publication builders."""

from __future__ import annotations

import os
import re
from datetime import datetime, timezone
from pathlib import Path


RELEASE_DATE_RE = re.compile(
    r"\\newcommand\{\\OPHPaperReleaseDate\}\{([^}]*)\}"
)
DATE_VALUE_RE = re.compile(r"([A-Za-z]+) ([0-9]{1,2}), ([0-9]{4})")
MONTHS = {
    name: month
    for month, name in enumerate(
        (
            "January",
            "February",
            "March",
            "April",
            "May",
            "June",
            "July",
            "August",
            "September",
            "October",
            "November",
            "December",
        ),
        start=1,
    )
}


def build_environment(repo_root: Path) -> dict[str, str]:
    """Return an environment whose PDF dates are fixed by the release record."""
    release_info = (repo_root / "paper" / "release_info.tex").read_text(
        encoding="utf-8"
    )
    match = RELEASE_DATE_RE.search(release_info)
    if not match:
        raise RuntimeError("paper/release_info.tex has no OPHPaperReleaseDate")
    date_match = DATE_VALUE_RE.fullmatch(match.group(1).strip())
    if not date_match or date_match.group(1) not in MONTHS:
        raise RuntimeError("OPHPaperReleaseDate must use 'Month D, YYYY'")
    month_name, day, year = date_match.groups()
    release_day = datetime(
        int(year),
        MONTHS[month_name],
        int(day),
        tzinfo=timezone.utc,
    )
    env = os.environ.copy()
    env["SOURCE_DATE_EPOCH"] = str(int(release_day.timestamp()))
    env["TZ"] = "UTC"
    return env
