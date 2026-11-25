"""Regression checks for a renamed task, the public draw and failed controls."""

from copy import deepcopy
import subprocess

import pytest

from . import experiment as e


def test_official_draw_uses_first_matching_version(tmp_path):
    def git(*args):
        return (
            subprocess.check_output(["git", "-C", str(tmp_path), *args])
            .decode()
            .strip()
        )

    git("init", "-q")
    git("config", "user.name", "Fixture")
    git("config", "user.email", "fixture@invalid.example")
    git("commit", "-qm", "baseline", "--allow-empty")
    base = git("rev-parse", "HEAD")
    spec = dict(e.terms(), platform_base=base)
    task = tmp_path / spec["platform_task"]
    task.parent.mkdir(parents=True)

    def publish(content):
        task.write_bytes(content)
        git("add", str(task))
        git("commit", "-qm", "task")
        git("update-ref", "refs/remotes/origin/main", "HEAD")
        return git("rev-parse", "HEAD")

    publish(b"different task")
    assert e.official_anchor(tmp_path, spec) is None
    first = publish((e.HERE / "TASK.md").read_bytes())
    publish(b"another task")
    publish((e.HERE / "TASK.md").read_bytes())
    assert e.official_anchor(tmp_path, spec) == first
    assert e.selection(first, spec) == e.selection(
        e.official_anchor(tmp_path, spec), spec
    )


def test_failed_control_and_changed_address_are_detected():
    settings = dict(home=13, start=0, direction=1, exposure=37)
    result = {}
    spec = e.terms()
    e.replay(settings, result, spec)
    assert result["passed"]
    other = deepcopy(result["bins"])
    other["coherent"][:12] = reversed(other["coherent"][:12])
    assert not e.compare(result["bins"], other, result["floor_ppt"], spec["gates"])[
        "passed"
    ]
    unchanged = deepcopy(result["bins"])
    unchanged["stochastic_swap"] = unchanged["coherent"]
    assert not e.compare(unchanged, unchanged, result["floor_ppt"], spec["gates"])[
        "passed"
    ]
    with pytest.raises(ValueError, match="Invalid"):
        broken = deepcopy(result["bins"])
        broken["coherent"][13] += 0.1
        e.compare(broken, broken, result["floor_ppt"], spec["gates"])
