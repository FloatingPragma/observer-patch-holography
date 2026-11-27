"""The lane table is complete, packs without loss, and names existing paths."""
from __future__ import annotations

import re
import shlex
import sys
from collections import Counter
from pathlib import Path

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
        assert lane.profile in {"root", "model", "einstein"}
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


def test_model_tests_and_replays_run_once_and_both_gate_ci():
    job = _job("model")
    suites = job["strategy"]["matrix"]["suite"]
    assert sorted(suites) == ["replays", "tests"] and job["strategy"]["fail-fast"] is False
    assert "model" in _job("ci-ok")["needs"]
    for step in job["steps"]:
        if "run" in step and not step["run"].startswith("python -m pip install"):
            expected = "tests" if step["run"] == "python -B -m pytest -q -W error model" else "replays"
            assert step["if"] == f"matrix.suite == '{expected}'"
    assert sum(s.get("if") == "matrix.suite == 'tests'" for s in job["steps"]) == 1


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


def test_einstein_jobs_own_every_matched_test_once_on_model_pins():
    job = _job("einstein")
    suites = job["strategy"]["matrix"]["suite"]
    commands = [step.get("run", "") for step in job["steps"]]
    assert "python -m pip install -r model/requirements.txt" in commands
    assert "python tools/ci_lanes.py run --lane einstein_${{ matrix.suite }}" in commands
    assert "einstein" in _job("ci-ok")["needs"]
    assert job["strategy"]["fail-fast"] is False
    lanes = ci_lanes.select("push", "einstein")
    assert {lane.name for lane in lanes} == {f"einstein_{name}" for name in suites}
    assert len(lanes) == len(suites)
    root = Path("code/matched_einstein_limit")
    existing = {p.relative_to(ROOT) for p in (ROOT / root).rglob("test_*.py")}
    # Previously unseen directories must also be selected, including one at the
    # repository package root and one under a carved-out expensive subpackage.
    future = {root / "future" / "test_future.py",
              root / "source_implication/directional_stress/regional_modular/future/test_future.py",
              root / "source_implication/directional_stress/regional_modular/station_memories/test_future.py"}
    owners = {p: [] for p in existing | future}
    for lane in lanes:
        assert lane.env == {**ci_lanes.CODE, **ci_lanes.BLAS1}
        (command,) = lane.commands
        tokens = shlex.split(command)
        assert tokens[:6] == ["python", "-m", "pytest", "-q", "-W", "error"]
        paths = [Path(t) for t in tokens[6:] if not t.startswith("--ignore=")]
        excluded = [Path(t.removeprefix("--ignore=")) for t in tokens[6:] if t.startswith("--ignore=")]
        for test in owners:
            hierarchy = {test, *test.parents}
            if hierarchy.intersection(paths) and not hierarchy.intersection(excluded):
                owners[test].append(lane.name)
    assert all(len(names) == 1 for names in owners.values()), owners
    assert set(n for names in owners.values() for n in names) == {lane.name for lane in lanes}
    assert not any("code/matched_einstein_limit" in command for lane in ci_lanes.LANES
                   if lane.profile != "einstein" for command in lane.commands)


def test_original_einstein_lane_replays_every_suite(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["ci_lanes.py", "list", "--lane", "matched_einstein_limit"])
    ci_lanes.main()
    listing = capsys.readouterr().out
    for lane in ci_lanes.select("push", "einstein"):
        assert listing.count(lane.name) == 1
