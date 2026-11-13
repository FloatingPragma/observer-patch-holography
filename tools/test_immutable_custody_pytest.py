"""Execution-count and routing mutation tests; actual science is checked separately."""
from pathlib import Path
import types

import pytest

from tools import immutable_custody_pytest as runner


def collected(root, ids=None, path=None):
    return types.SimpleNamespace(items=[types.SimpleNamespace(nodeid=node, path=path or root / runner.TEST)
                                       for node in sorted(runner.EXPECTED_IDS if ids is None else ids)])


def complete(root):
    plugin = runner.ExactTests(root)
    plugin.pytest_collection_finish(collected(root))
    for node in sorted(runner.EXPECTED_IDS):
        for phase in ("setup", "call", "teardown"):
            plugin.pytest_runtest_logreport(types.SimpleNamespace(nodeid=node, when=phase, outcome="passed"))
    return plugin


def test_accepts_only_complete_original_test_phase_set(tmp_path):
    plugin = complete(tmp_path)
    plugin.validate(0)
    assert len(plugin.evidence()["phase_reports"]) == 36


@pytest.mark.parametrize("change", ["omission", "duplicate", "foreign", "different_path"])
def test_collection_cannot_omit_duplicate_or_substitute_tests(tmp_path, change):
    ids = sorted(runner.EXPECTED_IDS)
    path = None
    if change == "omission":
        ids.pop()
    elif change == "duplicate":
        ids[-1] = ids[0]
    elif change == "foreign":
        ids[-1] = "model/unrelated/test.py::test_other"
    else:
        path = tmp_path / "wrong.py"
    with pytest.raises(ValueError, match="collection differs"):
        runner.ExactTests(tmp_path).pytest_collection_finish(collected(tmp_path, ids, path))


@pytest.mark.parametrize("mutation", ["failed", "skipped", "xfail", "duplicate", "missing_teardown", "exit_failure"])
def test_incomplete_or_nonpassing_execution_never_becomes_acceptance(tmp_path, mutation):
    plugin = complete(tmp_path)
    node = sorted(runner.EXPECTED_IDS)[0]
    if mutation == "missing_teardown":
        del plugin.reports[(node, "teardown")]
    elif mutation != "exit_failure":
        if mutation != "duplicate":
            del plugin.reports[(node, "call")]
        report = types.SimpleNamespace(nodeid=node, when="call", outcome="passed")
        if mutation in ("failed", "skipped"):
            report.outcome = mutation
        if mutation == "xfail":
            report.wasxfail = "must not be accepted"
        plugin.pytest_runtest_logreport(report)
    with pytest.raises(ValueError, match="not accepted"):
        plugin.validate(1 if mutation == "exit_failure" else 0)


def test_collection_error_is_retained_even_with_complete_test_results(tmp_path):
    plugin = complete(tmp_path)
    plugin.pytest_collectreport(types.SimpleNamespace(outcome="failed"))
    with pytest.raises(ValueError):
        plugin.validate(0)


def test_original_test_bytes_are_pinned(tmp_path):
    target = tmp_path / runner.TEST
    target.parent.mkdir(parents=True)
    target.write_text("def test_replacement(): pass\n")
    with pytest.raises(ValueError, match="renewed review"):
        runner.authenticate_test(tmp_path)


def test_original_test_symlink_cannot_escape_root(tmp_path):
    root = tmp_path / "root"
    target = root / runner.TEST
    target.parent.mkdir(parents=True)
    outside = tmp_path / "outside.py"
    outside.write_text("pass\n")
    target.symlink_to(outside)
    with pytest.raises(ValueError, match="renewed review"):
        runner.authenticate_test(root)


def test_reviewed_guard_is_authenticated_before_loading():
    guard, path = runner.load_guard()
    assert runner.sha(path) == runner.GUARD_SHA256
    assert guard.OWNER == "model.particle_identification.completion_certificate"


@pytest.mark.parametrize("variable", ["PYTEST_ADDOPTS", "PYTEST_PLUGINS"])
def test_environment_cannot_add_skip_or_replace_test_routing(tmp_path, monkeypatch, variable):
    monkeypatch.setenv(variable, "untrusted-extra-routing")
    with pytest.raises(SystemExit) as error:
        runner.main(["--root", str(tmp_path), "--audit", str(tmp_path.parent / "external-audit.json")])
    assert error.value.code == 2


def test_audit_cannot_overwrite_original_test_module(tmp_path):
    with pytest.raises(SystemExit) as error:
        runner.main(["--root", str(tmp_path), "--audit", str(tmp_path / runner.TEST)])
    assert error.value.code == 2


def harness(root):
    repository = Path(runner.__file__).resolve().parents[1]
    (root / "pytest.ini").write_bytes((repository / "pytest.ini").read_bytes())


@pytest.mark.parametrize("name", runner.ABSENT_CONFTESTS)
def test_added_ancestor_conftest_requires_review(tmp_path, name):
    harness(tmp_path)
    runner.authenticate_harness(tmp_path)
    path = tmp_path / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("# A newly supplied fixture hook\n")
    with pytest.raises(ValueError, match="conftest"):
        runner.authenticate_harness(tmp_path)


def test_changed_repository_pytest_configuration_requires_review(tmp_path):
    (tmp_path / "pytest.ini").write_text("[pytest]\naddopts = -k skip_original\n")
    with pytest.raises(ValueError, match="harness"):
        runner.authenticate_harness(tmp_path)


@pytest.mark.parametrize("failure", ["missing_receipt", "escaping_source"])
def test_initial_validation_failure_does_not_leak_runtime_state(tmp_path, monkeypatch, failure):
    import os
    import subprocess
    import sys

    root = tmp_path / "repo"
    root.mkdir()
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    harness(root)
    target = root / runner.TEST
    target.parent.mkdir(parents=True)
    original = Path(runner.__file__).resolve().parents[1] / runner.TEST
    target.write_bytes(original.read_bytes())
    if failure == "escaping_source":
        outside = tmp_path / "outside.py"
        outside.write_text("pass\n")
        (root / "model/escape.py").symlink_to(outside)
    monkeypatch.delenv("PYTEST_ADDOPTS", raising=False)
    monkeypatch.delenv("PYTEST_PLUGINS", raising=False)
    monkeypatch.setenv("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "original-value")
    before = (Path.cwd(), list(sys.path), list(sys.argv), os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"])
    exception = FileNotFoundError if failure == "missing_receipt" else ValueError
    with pytest.raises(exception):
        runner.main(["--root", str(root), "--audit", str(tmp_path / "audit.json")])
    assert (Path.cwd(), list(sys.path), list(sys.argv), os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"]) == before
