"""Exercise the runner's actual argument parsing and subprocess dispatch."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest

from tools import run_mandatory_suite as runner


def _standard_steps():
    return [step for step in runner.MANDATORY_STEPS if step[0] not in runner.HEAVY_STEP_TITLES]


@pytest.fixture
def dispatch(monkeypatch):
    commands = []

    def capture(command, *, cwd):
        assert cwd == runner.ROOT
        commands.append(list(command))
        return SimpleNamespace(returncode=0)

    # Keep main(), run_steps(), and their ordering/error checks in the path.
    # Only the child scientific programs are replaced by a dispatch recorder.
    monkeypatch.setattr(runner.subprocess, "run", capture)

    def invoke(*arguments):
        commands.clear()
        monkeypatch.setattr(sys, "argv", ["run_mandatory_suite.py", *arguments])
        runner.main()
        return list(commands)

    return invoke, commands


@pytest.mark.parametrize("count", [1, 2, 3, 7, 9, 10, len(_standard_steps())])
def test_all_standard_shards_dispatch_exact_ordered_coverage(dispatch, capsys, count):
    invoke, _ = dispatch
    expected = [command for _, command in _standard_steps()]
    partitions = []
    for index in range(count):
        partitions.append(invoke("--shard-index", str(index), "--shard-count", str(count)))
        output = capsys.readouterr().out
        assert f"standard shard {index + 1}/{count}:" in output
        assert f"mandatory shard OK ({index + 1}/{count}; standard mode)" in output
        assert "mandatory suite OK" not in output
    assert all(partitions)
    assert max(map(len, partitions)) - min(map(len, partitions)) <= 1
    # Exact concatenation checks both contiguous order and multiplicity.
    assert [command for partition in partitions for command in partition] == expected


def test_duplicate_titles_and_commands_survive_sharding(monkeypatch, dispatch):
    steps = [
        ("duplicate", ["command", "same"]),
        ("duplicate", ["command", "different"]),
        ("heavy", ["excluded", "only-from-standard"]),
        ("duplicate", ["command", "same"]),
        ("last", ["command", "last"]),
    ]
    monkeypatch.setattr(runner, "MANDATORY_STEPS", steps)
    monkeypatch.setattr(runner, "HEAVY_STEP_TITLES", frozenset({"heavy"}))
    invoke, _ = dispatch
    first = invoke("--shard-index", "0", "--shard-count", "2")
    second = invoke("--shard-index", "1", "--shard-count", "2")
    assert first == [["command", "same"], ["command", "different"]]
    assert second == [["command", "same"], ["command", "last"]]


@pytest.mark.parametrize(
    ("arguments", "mode"),
    [
        ([], "standard"),
        (["--full"], "full"),
        (["--certificates"], "standard-certificates"),
        (["--full", "--certificates"], "full-certificates"),
        (["--certificates-only"], "certificates"),
        (["--certificate-smoke-only"], "smoke"),
    ],
)
def test_unsharded_modes_keep_their_complete_command_lists(dispatch, capsys, arguments, mode):
    expected = {
        "standard": _standard_steps(),
        "full": runner.MANDATORY_STEPS,
        "standard-certificates": _standard_steps() + runner.CERTIFICATE_STEPS,
        "full-certificates": runner.MANDATORY_STEPS + runner.CERTIFICATE_STEPS,
        "certificates": runner.CERTIFICATE_STEPS,
        "smoke": runner.CERTIFICATE_SMOKE_STEPS,
    }[mode]
    invoke, _ = dispatch
    assert invoke(*arguments) == [command for _, command in expected]
    assert "shard OK" not in capsys.readouterr().out


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (["--shard-index", "0"], "must be supplied together"),
        (["--shard-count", "2"], "must be supplied together"),
        (["--shard-index", "-1", "--shard-count", "2"], "0 <= index"),
        (["--shard-index", "2", "--shard-count", "2"], "0 <= index"),
        (["--shard-index", "0", "--shard-count", "0"], "shards must be nonempty"),
        (["--shard-index", "0", "--shard-count", "-1"], "shards must be nonempty"),
        (["--shard-index", "0", "--shard-count", str(len(_standard_steps()) + 1)], "shards must be nonempty"),
        (["--shard-index", "0.5", "--shard-count", "2"], "invalid int value"),
        (["--shard-index", "0", "--shard-count", "two"], "invalid int value"),
        (["--shard-index", "0", "--shard-count", "2", "--full"], "only for the standard suite"),
        (["--shard-index", "0", "--shard-count", "2", "--certificates"], "only for the standard suite"),
        (["--shard-index", "0", "--shard-count", "2", "--certificates-only"], "only for the standard suite"),
        (["--shard-index", "0", "--shard-count", "2", "--certificate-smoke-only"], "only for the standard suite"),
    ],
)
def test_invalid_shard_arguments_fail_before_any_child_runs(dispatch, capsys, arguments, message):
    invoke, commands = dispatch
    with pytest.raises(SystemExit) as error:
        invoke(*arguments)
    assert error.value.code == 2
    assert message in capsys.readouterr().err
    assert commands == []


def test_empty_standard_list_cannot_report_a_successful_shard(monkeypatch, dispatch, capsys):
    monkeypatch.setattr(runner, "MANDATORY_STEPS", [])
    monkeypatch.setattr(runner, "HEAVY_STEP_TITLES", frozenset())
    invoke, commands = dispatch
    with pytest.raises(SystemExit) as error:
        invoke("--shard-index", "0", "--shard-count", "1")
    assert error.value.code == 2
    assert "shards must be nonempty" in capsys.readouterr().err
    assert commands == []


def test_sharding_keeps_unknown_heavy_step_guard(monkeypatch, dispatch):
    monkeypatch.setattr(runner, "HEAVY_STEP_TITLES", frozenset({"missing step"}))
    invoke, commands = dispatch
    with pytest.raises(SystemExit, match="HEAVY_STEP_TITLES entries missing"):
        invoke("--shard-index", "0", "--shard-count", "2")
    assert commands == []


def test_failed_child_stops_its_shard_and_cannot_report_success(monkeypatch, capsys):
    monkeypatch.setattr(runner, "MANDATORY_STEPS", [
        ("first shard", ["unused"]),
        ("first shard again", ["unused-too"]),
        ("failing child", ["fails"]),
        ("must not run", ["later"]),
    ])
    monkeypatch.setattr(runner, "HEAVY_STEP_TITLES", frozenset())
    commands = []

    def failure(command, *, cwd):
        assert cwd == runner.ROOT
        commands.append(command)
        return SimpleNamespace(returncode=7)

    monkeypatch.setattr(runner.subprocess, "run", failure)
    monkeypatch.setattr(sys, "argv", [
        "run_mandatory_suite.py", "--shard-index", "1", "--shard-count", "2",
    ])
    with pytest.raises(SystemExit, match=r"FAILED: failing child \(exit 7\)"):
        runner.main()
    assert commands == [["fails"]]
    assert "shard OK" not in capsys.readouterr().out


@pytest.mark.parametrize("count", [1, 10])
def test_long_steps_run_once_between_skipping_shards_and_their_own_job(dispatch, capsys, count):
    invoke, _ = dispatch
    standard = _standard_steps()
    long_commands = [command for title, command in standard if title in runner.LONG_STEP_TITLES]
    assert len(long_commands) == len(runner.LONG_STEP_TITLES) >= 1
    partitions = [invoke("--shard-index", str(index), "--shard-count", str(count), "--skip-long")
                  for index in range(count)]
    assert all(partitions)
    long_only = invoke("--long-only")
    assert "mandatory long steps OK" in capsys.readouterr().out
    assert long_only == long_commands
    shards = [command for partition in partitions for command in partition]
    assert not any(command in long_commands for command in shards)
    # The shards and the long job together are exactly the standard list.
    assert sorted(map(tuple, shards + long_only)) == sorted(tuple(command) for _, command in standard)


@pytest.mark.parametrize("arguments,message", [
    (["--skip-long", "--long-only"], "mutually exclusive"),
    (["--long-only", "--full"], "only to the standard suite"),
    (["--skip-long", "--certificates-only"], "only to the standard suite"),
    (["--long-only", "--shard-index", "0", "--shard-count", "2"], "cannot be sharded"),
])
def test_long_step_options_reject_other_modes(dispatch, capsys, arguments, message):
    invoke, commands = dispatch
    with pytest.raises(SystemExit):
        invoke(*arguments)
    assert message in capsys.readouterr().err
    assert commands == []


def test_unknown_long_step_title_fails_closed(monkeypatch, dispatch):
    monkeypatch.setattr(runner, "LONG_STEP_TITLES", frozenset({"missing step"}))
    invoke, commands = dispatch
    with pytest.raises(SystemExit, match="LONG_STEP_TITLES entries missing"):
        invoke("--long-only")
    assert commands == []
