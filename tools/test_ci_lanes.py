"""The lane table is complete, packs without loss, and names existing paths."""
from __future__ import annotations

import re
import shlex
from collections import Counter
from pathlib import Path

import pytest
import yaml

from tools import ci_lanes

ROOT = Path(__file__).resolve().parents[1]


def _matrix(workflow: str, job: str) -> list[dict]:
    data = yaml.safe_load((ROOT / ".github/workflows" / workflow).read_text(encoding="utf-8"))
    return data["jobs"][job]["strategy"]["matrix"]["include"]


def test_lane_names_are_unique_and_classified():
    names = Counter(lane.name for lane in ci_lanes.LANES)
    assert not [name for name, count in names.items() if count > 1]
    for lane in ci_lanes.LANES:
        assert lane.tier in {"push", "nightly"}
        assert lane.profile in {"root", "model"}
        assert lane.commands and lane.minutes > 0


@pytest.mark.parametrize(("workflow", "job", "tier"), [("ci.yml", "lanes", "push"), ("nightly.yml", "lanes", "nightly")])
def test_workflow_buckets_run_every_lane_once(workflow, job, tier):
    include = _matrix(workflow, job)
    by_profile: dict[str, set[int]] = {}
    counts: dict[str, int] = {}
    for row in include:
        by_profile.setdefault(row["profile"], set()).add(row["bucket"])
        counts[row["profile"]] = row["buckets"]
    expected = {lane.profile for lane in ci_lanes.LANES if lane.tier == tier}
    assert set(by_profile) == expected
    for profile, buckets in by_profile.items():
        assert buckets == set(range(counts[profile]))
        packed = ci_lanes.pack(ci_lanes.select(tier, profile), counts[profile])
        flat = [lane.name for group in packed for lane in group]
        assert sorted(flat) == sorted(lane.name for lane in ci_lanes.select(tier, profile))
        assert all(packed)


def test_every_referenced_repository_path_exists():
    missing = []
    for lane in ci_lanes.LANES:
        for command in lane.commands:
            for token in shlex.split(command):
                path = token.split("::", 1)[0].removeprefix("--ignore=")
                if re.match(r"^(code|model|tools|evidence|Lean)/", path) and "$" not in path and "{" not in path:
                    if not (ROOT / path).exists():
                        missing.append((lane.name, path))
    assert not missing
