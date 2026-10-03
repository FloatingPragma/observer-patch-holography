"""Replay the workflow's shard routing and its fail-closed status aggregate."""
from __future__ import annotations

from collections import Counter
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from types import SimpleNamespace

import pytest
import yaml

from tools import run_mandatory_suite as runner

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github/workflows/mandatory-suite.yml"
ADDITIONAL_EVIDENCE_TESTS = (
    "code/P_derivation/test_selection_accounting.py",
    "code/P_derivation/test_printed_pair_identity.py",
    "evidence/source_net_causal_poset/manifold_refinement/test_verify.py",
    "code/particles/flavor/test_full_entropy_w5_selector.py",
    "code/particles/flavor/test_conditional_quark_mass_replay.py",
    "code/particles/flavor/test_verify_conditional_quark_mass_replay.py",
    "code/particles/flavor/test_compare_conditional_quark_masses.py",
    "code/particles/calibration/test_conditional_top_fixed_p.py",
    "code/particles/flavor/test_source_w5_response_constraints.py",
    "code/particles/flavor/test_native_repair_flavor_constraints.py",
    "code/particles/calibration/test_top_positive_capacity_band.py",
    "code/particles/flavor/test_native_record_renewal.py",
    "code/particles/flavor/test_exterior_higgs_response.py",
    "code/particles/calibration/test_source_ew_vev_matching.py",
    "code/particles/calibration/test_verify_source_ew_vev_matching.py",
    "code/particles/calibration/test_compare_source_ew_vev_matching.py",
)


def _workflow():
    # BaseLoader preserves GitHub's `on` key rather than YAML 1.1's boolean key.
    return yaml.load(WORKFLOW.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)


def _step(job, name):
    return next(step for step in job["steps"] if step.get("name") == name)


def _shards(job, full):
    expression = job["strategy"]["matrix"]["shard"]
    match = re.fullmatch(
        r"\$\{\{ fromJSON\(inputs\.full && '([^']+)' \|\| '([^']+)'\) \}\}",
        expression,
    )
    assert match, "unsupported worker matrix expression"
    return json.loads(match[1 if full else 2])


def _arguments(job, full, shard):
    command = _step(job, "Run the mandatory suite")["run"]
    match = re.fullmatch(
        r"python tools/run_mandatory_suite\.py "
        r"\$\{\{ inputs\.full && '([^']+)' \|\| "
        r"format\('([^']+)', matrix\.shard\) \}\}",
        command,
    )
    assert match, "unsupported mandatory command expression"
    return shlex.split(match[1] if full else match[2].format(shard))


def _evidence_dispatch(job, shard=None):
    """Read actual workflow commands; reject collection-only or partial execution."""
    dispatched = []
    for step in job["steps"]:
        for command in step.get("run", "").splitlines():
            if not any(path in command for path in ADDITIONAL_EVIDENCE_TESTS):
                continue
            assert not step.get("continue-on-error")
            if step.get("if"):
                assert shard is not None, "nightly evidence must be unconditional"
                assert step["if"] == "matrix.shard == 0"
                if shard != 0:
                    continue
            arguments = shlex.split(command)
            assert arguments[:4] == ["python", "-m", "pytest", "-q"]
            assert arguments[4:]
            # Only whole test files are allowed: no -k selection, collect-only,
            # shell success fallback, or ignored mutation tests.
            assert all(path in ADDITIONAL_EVIDENCE_TESTS for path in arguments[4:])
            assert all((ROOT / path).is_file() for path in arguments[4:])
            dispatched.extend(arguments[4:])
    return dispatched


def _assert_evidence_coverage(workflow, full, operating_system):
    job = workflow["jobs"]["mandatory-shards"]
    assert operating_system in job["strategy"]["matrix"]["os"]
    dispatched = [
        path for shard in _shards(job, full)
        for path in _evidence_dispatch(job, shard)
    ]
    assert Counter(dispatched) == Counter(ADDITIONAL_EVIDENCE_TESTS)
    assert Counter(_evidence_dispatch(workflow["jobs"]["mandatory-full"])) == (
        Counter(ADDITIONAL_EVIDENCE_TESTS)
    )


@pytest.mark.parametrize("full", [False, True])
@pytest.mark.parametrize("operating_system", ["ubuntu-latest", "windows-latest"])
def test_new_evidence_executes_once_per_os_and_in_nightly(full, operating_system):
    _assert_evidence_coverage(_workflow(), full, operating_system)


@pytest.mark.parametrize("job_name", ["mandatory-shards", "mandatory-full"])
@pytest.mark.parametrize(
    "mutation", ["omit_test", "duplicate", "collect_only", "skip", "allow_failure"]
)
def test_new_evidence_coverage_rejects_false_green_routes(job_name, mutation):
    workflow = _workflow()
    job = workflow["jobs"][job_name]
    step = _step(job, "Verify additional scientific evidence")
    if mutation == "omit_test":
        step["run"] = step["run"].replace(ADDITIONAL_EVIDENCE_TESTS[-1], "")
    elif mutation == "duplicate":
        job["steps"].append(dict(step))
    elif mutation == "collect_only":
        step["run"] += " --collect-only"
    elif mutation == "skip":
        step["if"] = "false"
    else:
        step["continue-on-error"] = "true"
    with pytest.raises(AssertionError):
        _assert_evidence_coverage(workflow, False, "ubuntu-latest")


def test_new_evidence_cannot_be_repeated_on_every_push_shard():
    workflow = _workflow()
    step = _step(workflow["jobs"]["mandatory-shards"], "Verify additional scientific evidence")
    del step["if"]
    with pytest.raises(AssertionError):
        _assert_evidence_coverage(workflow, False, "ubuntu-latest")


@pytest.mark.parametrize("full", [False, True])
@pytest.mark.parametrize("operating_system", ["ubuntu-latest", "windows-latest"])
def test_workflow_matrix_dispatches_the_entire_selected_suite_once(
    monkeypatch, full, operating_system
):
    job = _workflow()["jobs"]["mandatory-shards"]
    assert operating_system in job["strategy"]["matrix"]["os"]
    assert _shards(job, full) == ([0] if full else list(range(9)))
    actual = []

    def record(command, *, cwd):
        assert cwd == ROOT
        actual.append(command)
        return SimpleNamespace(returncode=0)

    # Exercise the arguments from actual YAML through main/run_steps. Only the
    # expensive child scientific programs are replaced by a dispatch recorder.
    monkeypatch.setattr(runner.subprocess, "run", record)
    for shard in _shards(job, full):
        arguments = _arguments(job, full, shard)
        if full:
            assert arguments == ["--full"]
        else:
            assert arguments == ["--shard-index", str(shard), "--shard-count", "9"]
        monkeypatch.setattr(sys, "argv", ["run_mandatory_suite.py", *arguments])
        runner.main()
    expected = [
        command for title, command in runner.MANDATORY_STEPS
        if full or title not in runner.HEAVY_STEP_TITLES
    ]
    assert actual == expected  # Exact ordered multiplicity, including repeats.


def test_each_worker_has_a_fresh_complete_checkout_and_the_runtime_ceiling():
    job = _workflow()["jobs"]["mandatory-shards"]
    assert job["if"] == "github.event_name != 'schedule'"
    assert job["runs-on"] == "${{ matrix.os }}"
    assert job["timeout-minutes"] == "30"
    assert job["strategy"]["fail-fast"] == "false"
    assert job["strategy"]["matrix"]["os"] == ["ubuntu-latest", "windows-latest"]
    assert set(job["strategy"]["matrix"]) == {"os", "shard"}
    assert not job.get("continue-on-error")
    checkout = _step(job, "Checkout")
    assert checkout["uses"].startswith("actions/checkout@")
    assert checkout["with"]["fetch-depth"] == "0"
    assert _step(job, "Set up Python")["with"]["python-version"] == "3.12"
    install = _step(job, "Install pinned dependencies")["run"]
    assert "pip install -r requirements.txt" in install
    assert _step(job, "Regression-test the Lean CI runtime ceiling")["run"] == (
        "python -m pytest -q tools/test_check_lean_ci_budget.py"
    )
    assert _step(job, "Regression-test mandatory sharding and workflow coverage")["run"] == (
        "python -m pytest -q tools/test_mandatory_suite_shards.py "
        "tools/test_mandatory_suite_workflow.py"
    )
    for step in job["steps"]:
        assert not step.get("continue-on-error")
        if step.get("name") == "Verify additional scientific evidence":
            assert step["if"] == "matrix.shard == 0"
        else:
            assert not step.get("if")
    assert not any("download-artifact@" in step.get("uses", "") for step in job["steps"])


def test_historical_statuses_require_all_workers_even_after_failure_or_skip():
    job = _workflow()["jobs"]["mandatory"]
    assert job["name"] == "mandatory collection (${{ matrix.os }})"
    assert job["strategy"]["matrix"]["os"] == ["ubuntu-latest", "windows-latest"]
    assert job["if"] == "${{ always() && github.event_name != 'schedule' }}"
    assert job["needs"] == ["mandatory-shards"]
    assert job["runs-on"] == "ubuntu-latest"
    assert int(job["timeout-minutes"]) <= 30
    assert not job.get("continue-on-error")
    assert len(job["steps"]) == 1
    step = job["steps"][0]
    assert step["shell"] == "python"
    assert not step.get("continue-on-error")
    assert not step.get("if")
    assert step["env"] == {"MANDATORY_RESULT": "${{ needs.mandatory-shards.result }}"}


@pytest.mark.parametrize(
    "result", ["success", "failure", "cancelled", "skipped", "", None, "Success", "success\n"]
)
def test_actual_aggregate_python_only_accepts_exact_success(result):
    step = _workflow()["jobs"]["mandatory"]["steps"][0]
    assert step["shell"] == "python"
    env = os.environ.copy()
    env.pop("MANDATORY_RESULT", None)
    if result is not None:
        env["MANDATORY_RESULT"] = result
    # Run exactly the Python script executed by GitHub, including missing env.
    completed = subprocess.run(
        [sys.executable, "-c", step["run"]], env=env,
        text=True, capture_output=True, check=False,
    )
    assert (completed.returncode == 0) is (result == "success")
    if result != "success":
        assert "did not all succeed" in completed.stderr


def test_nightly_full_and_manual_full_are_retained_without_path_filters():
    workflow = _workflow()
    triggers = workflow["on"]
    assert triggers["push"] == {"branches": ["main"]}
    assert triggers["pull_request"] == ""
    assert triggers["schedule"] == [{"cron": "17 3 * * *"}]
    assert triggers["workflow_dispatch"]["inputs"]["full"]["type"] == "boolean"
    assert triggers["workflow_dispatch"]["inputs"]["full"]["default"] == "false"
    job = workflow["jobs"]["mandatory-full"]
    assert job["if"] == "github.event_name == 'schedule'"
    assert job["runs-on"] == "ubuntu-latest"
    assert "strategy" not in job
    assert _step(job, "Run the full mandatory suite")["run"] == (
        "python tools/run_mandatory_suite.py --full"
    )
    assert _step(job, "Checkout")["with"]["fetch-depth"] == "0"
