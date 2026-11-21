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


def _job(job: str) -> dict:
    data = yaml.safe_load((ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8"))
    return data["jobs"][job]


def test_lane_names_are_unique_and_classified():
    names = Counter(lane.name for lane in ci_lanes.LANES)
    assert not [name for name, count in names.items() if count > 1]
    for lane in ci_lanes.LANES:
        assert lane.tier == "push"
        assert lane.profile in {"root", "model"}
        assert lane.commands and lane.minutes > 0


def test_ci_runs_every_lane_once():
    job = _job("certificates")
    shards = job["strategy"]["matrix"]["shard"]
    count = int(job["env"]["SHARDS"])
    assert shards == list(range(count))
    steps = "\n".join(step.get("run", "") for step in job["steps"])
    assert "ci_lanes.py run --tier push --profile root --buckets $SHARDS --bucket ${{ matrix.shard }}" in steps
    packed = ci_lanes.pack(ci_lanes.select("push", "root"), count)
    assert sorted(lane.name for group in packed for lane in group) == sorted(
        lane.name for lane in ci_lanes.select("push", "root"))
    assert all(packed)
    model = "\n".join(step.get("run", "") for step in _job("model")["steps"])
    assert "ci_lanes.py run --tier push --profile model --buckets 1 --bucket 0" in model


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
