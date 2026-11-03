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
    "code/geometry/test_collar_scaling_boundaries.py",
    "code/thermodynamics/test_collar_stationary_regressions.py",
    "code/thermodynamics/test_exact_count_chain_controls.py",
    "evidence/exact_federation_tower_20260924/test_streamed_expectations.py",
    "evidence/observer_dynamics_20260925/lean/test_verify.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/boltzmann/test_forward_contract.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/boltzmann_data/test_audit_inputs.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/causal/test_experiment.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/dynamics/test_dynamics.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/inflation/native_history/test_native_history.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/necessity/boundary_bulk/test_boundary_bulk.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/necessity/event_records/test_event_records.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/necessity/temporal_filters/test_temporal_filters.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/observation/test_conditional_bridges.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/observation/test_observations.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/observer_cmb/test_measured_comparison.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/observer_cmb/test_observer_spectrum.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/primordial/test_sources.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/source_derivation/innovations/test_extract.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/source_derivation/innovations/test_nested.py",
    "evidence/observer_dynamics_20260925/sim-analysis/codex/source_derivation/native/test_native_source.py",
    "evidence/observer_dynamics_20260925/sim-analysis/scripts/test_codex_causal_cosmology_audit.py",
    "evidence/observer_dynamics_20260925/sim-analysis/tests/test_codex_golden_fourier.py",
    "evidence/observer_dynamics_20260925/sim-analysis/tests/test_codex_scale_ensemble.py",
    "evidence/source_net_causal_poset/manifold_sampled/tests/test_exact_source_net_manifold_sampled.py",
    "evidence/source_net_causal_poset/manifold_sampled/tests/test_retained_receipt_mutations.py",
)

# Bounded but heavy: the eight-test tower replay takes about a minute, so it runs
# in the nightly mandatory-full job only instead of on every push.
NIGHTLY_EVIDENCE_TESTS = (
    "evidence/exact_federation_tower_20260924/test_verify_tower.py",
)

# Evidence test files that another workflow already executes on every relevant
# change. They are not duplicated here.
OTHER_WORKFLOW_EVIDENCE_TESTS = {
    "evidence/federation_recovery/test_recovery.py": ".github/workflows/federation-recovery.yml",
}

# Evidence test files that cannot run from a clean checkout. Their disposition
# is tracked on the standing audit issue instead of being silently skipped.
DISPOSITIONED_EVIDENCE_TESTS = {
    "evidence/observer_dynamics_20260925/sim-analysis/codex/necessity/universality/test_universality.py": (
        "six tests need commit 14d1699 of the upstream oph-physics-sim history, "
        "which this repository does not carry (tracked on #1033)"
    ),
}


# CI runs on Ubuntu only, for every trigger.
PLATFORMS = ["ubuntu-latest"]


def _operating_systems(job, event, ref="refs/heads/main"):
    assert job["strategy"]["matrix"]["os"] == PLATFORMS, "unsupported OS matrix"
    return PLATFORMS


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
    gated = ADDITIONAL_EVIDENCE_TESTS + NIGHTLY_EVIDENCE_TESTS
    dispatched = []
    for step in job["steps"]:
        for command in step.get("run", "").splitlines():
            if not command.strip() or command.lstrip().startswith("#"):
                continue
            if not any(path in command for path in gated):
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
            assert all(path in gated for path in arguments[4:])
            assert all((ROOT / path).is_file() for path in arguments[4:])
            dispatched.extend(arguments[4:])
    return dispatched


def _assert_evidence_coverage(workflow, full, operating_system):
    job = workflow["jobs"]["mandatory-shards"]
    assert operating_system in _operating_systems(job, "workflow_dispatch")
    dispatched = [
        path for shard in _shards(job, full)
        for path in _evidence_dispatch(job, shard)
    ]
    assert Counter(dispatched) == Counter(ADDITIONAL_EVIDENCE_TESTS)
    assert Counter(_evidence_dispatch(workflow["jobs"]["mandatory-full"])) == (
        Counter(ADDITIONAL_EVIDENCE_TESTS + NIGHTLY_EVIDENCE_TESTS)
    )


@pytest.mark.parametrize("full", [False, True])
@pytest.mark.parametrize("operating_system", ["ubuntu-latest"])
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
        # Drop the last gated test wherever the workflow dispatches it, including
        # the sampled-manifold step.
        for candidate in job["steps"]:
            if "run" in candidate:
                candidate["run"] = candidate["run"].replace(
                    ADDITIONAL_EVIDENCE_TESTS[-1], ""
                )
    elif mutation == "duplicate":
        job["steps"].append(dict(step))
    elif mutation == "collect_only":
        step["run"] = step["run"].rstrip() + " --collect-only"
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


def test_heavy_evidence_runs_nightly_only():
    workflow = _workflow()
    job = workflow["jobs"]["mandatory-shards"]
    for shard in _shards(job, False):
        assert not set(_evidence_dispatch(job, shard)) & set(NIGHTLY_EVIDENCE_TESTS)
    assert set(NIGHTLY_EVIDENCE_TESTS) <= set(
        _evidence_dispatch(workflow["jobs"]["mandatory-full"])
    )


def _tracked_evidence_tests():
    output = subprocess.check_output(
        ["git", "ls-files", "evidence/**/test_*.py"], cwd=ROOT, text=True
    )
    return {line for line in output.splitlines() if line}


def test_no_evidence_test_file_is_silently_unreached():
    tracked = _tracked_evidence_tests()
    groups = {
        "gated": set(ADDITIONAL_EVIDENCE_TESTS) | set(NIGHTLY_EVIDENCE_TESTS),
        "another workflow": set(OTHER_WORKFLOW_EVIDENCE_TESTS),
        "dispositioned": set(DISPOSITIONED_EVIDENCE_TESTS),
    }
    claims = Counter()
    for name, paths in groups.items():
        claims.update(paths)
    assert [path for path, count in claims.items() if count > 1] == [], (
        "an evidence test must belong to exactly one group"
    )
    claimed_evidence = {path for path in claims if path.startswith("evidence/")}
    assert claimed_evidence == tracked, (
        "every tracked evidence test must be gated, owned by another workflow, "
        "or explicitly dispositioned instead of being skipped silently"
    )
    assert all((ROOT / path).is_file() for path in tracked)


@pytest.mark.parametrize("full", [False, True])
@pytest.mark.parametrize("operating_system", ["ubuntu-latest"])
def test_workflow_matrix_dispatches_the_entire_selected_suite_once(
    monkeypatch, full, operating_system
):
    job = _workflow()["jobs"]["mandatory-shards"]
    assert operating_system in _operating_systems(job, "workflow_dispatch")
    assert _shards(job, full) == ([0] if full else list(range(10)))
    actual = []
    # These complete command groups exhausted the 30-minute Windows budget
    # when the nine-shard layout put both in its last worker.
    expensive_markers = (
        "code/a5_fingerprint/test_a5_multipole_fixed_point_certificate.py",
        "code/a5_closure/tests/test_equal_state_weights_certificate.py",
    )
    assigned_workers = {marker: [] for marker in expensive_markers}

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
            assert arguments == ["--shard-index", str(shard), "--shard-count", "10"]
        monkeypatch.setattr(sys, "argv", ["run_mandatory_suite.py", *arguments])
        start = len(actual)
        runner.main()
        for command in actual[start:]:
            for marker in expensive_markers:
                if marker in command:
                    assigned_workers[marker].append(shard)
    expected = [
        command for title, command in runner.MANDATORY_STEPS
        if full or title not in runner.HEAVY_STEP_TITLES
    ]
    assert actual == expected  # Exact ordered multiplicity, including repeats.
    assert all(len(owners) == 1 for owners in assigned_workers.values())
    if not full:
        assert len({owners[0] for owners in assigned_workers.values()}) == 2


def test_each_worker_has_a_fresh_complete_checkout_and_the_runtime_ceiling():
    job = _workflow()["jobs"]["mandatory-shards"]
    assert "if" not in job
    assert job["runs-on"] == "${{ matrix.os }}"
    assert job["timeout-minutes"] == "30"
    assert job["strategy"]["fail-fast"] == "false"
    for event, ref in (("push", "refs/heads/main"), ("pull_request", "refs/heads/main"),
                       ("schedule", "refs/heads/main"), ("push", "refs/tags/r2057"),
                       ("workflow_dispatch", "refs/heads/main")):
        assert _operating_systems(job, event, ref) == ["ubuntu-latest"]
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
        if step.get("name") in {
            "Verify additional scientific evidence",
            "Verify sampled-manifold evidence",
            "Verify numerical recovery bounds without optional hardware SDKs",
        }:
            assert step["if"] == "matrix.shard == 0"
        else:
            assert not step.get("if")
    assert not any("download-artifact@" in step.get("uses", "") for step in job["steps"])


def test_historical_statuses_require_all_workers_even_after_failure_or_skip():
    job = _workflow()["jobs"]["mandatory"]
    assert job["name"] == "mandatory collection (${{ matrix.os }})"
    assert job["strategy"]["matrix"]["os"] == PLATFORMS
    assert job["if"] == "${{ always() }}"
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
    assert triggers["push"] == {"branches": ["main"], "tags": ["r*"]}
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
