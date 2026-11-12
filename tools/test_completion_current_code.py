"""The documentary view preserves current code and every scientific check."""
import hashlib
from pathlib import Path
import subprocess
import sys

import pytest
from tools import check_completion_current_code as current


def git(root, *args):
    return subprocess.run(["git", *args], cwd=root, check=True, capture_output=True, text=True)


@pytest.fixture
def source(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    git(root, "init", "-q")
    git(root, "config", "user.email", "fixture@example.invalid")
    git(root, "config", "user.name", "Fixture")
    (root / "test_probe.py").write_text("def test_current_code(): assert True\n")
    (root / "receipt.json").write_text('{"mass": 7, "accepted": false}\n')
    docs = current.authenticated_documents(current.ROOT)
    for name in docs:
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("Current hypothesis with distinct semantics.\n")
    git(root, "add", ".")
    git(root, "commit", "-qm", "fixture")
    return root, docs


def test_only_exact_two_git_authenticated_documentary_inputs():
    docs = current.authenticated_documents(current.ROOT)
    assert set(docs) == {
        "docs/OPH_EXECUTABLE_HYPOTHESIS.md", "docs/research/OPH_FIELD_RECONSTRUCTION.md"}
    assert all(hashlib.sha256(data).hexdigest() == current.DOCUMENTS[name]["sha256"]
               for name, data in docs.items())


def test_current_modified_code_is_executed_and_failure_is_not_hidden(source, tmp_path):
    root, docs = source
    (root / "test_probe.py").write_text("def test_current_code(): assert False, 'CURRENT CODE EXECUTED'\n")
    files = current.current_files(root)
    view = tmp_path / "view"
    current.install_current_inputs(root, view, files, docs)
    result = subprocess.run([sys.executable, "-m", "pytest", "-q", "test_probe.py"],
                            cwd=view, env=current.environment(view), capture_output=True, text=True)
    assert result.returncode == 1
    assert "CURRENT CODE EXECUTED" in result.stdout
    assert (view / "receipt.json").read_bytes() == (root / "receipt.json").read_bytes()
    assert (root / next(iter(docs))).read_text() == "Current hypothesis with distinct semantics.\n"
    current.verify_view(view, files)


@pytest.mark.parametrize("damage", ["added", "removed", "altered"])
def test_unlisted_or_modified_replacement_is_rejected(source, tmp_path, damage):
    root, docs = source
    if damage == "added":
        docs["test_probe.py"] = b"pass\n"
    elif damage == "removed":
        docs.pop(next(iter(docs)))
    else:
        docs[next(iter(docs))] += b"\n"
    with pytest.raises(ValueError, match="replacement inventory|unauthenticated"):
        current.install_current_inputs(root, tmp_path / "view", current.current_files(root), docs)


@pytest.mark.parametrize("name", ["test_probe.py", "receipt.json", *current.DOCUMENTS])
def test_private_view_mutation_of_code_science_or_document_is_rejected(source, tmp_path, name):
    root, docs = source
    files = current.current_files(root)
    view = tmp_path / "view"
    current.install_current_inputs(root, view, files, docs)
    path = view / name
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="view changed"):
        current.verify_view(view, files)


def test_concurrent_current_source_change_is_rejected_before_copy(source, tmp_path):
    root, docs = source
    files = current.current_files(root)
    (root / "test_probe.py").write_text("def test_changed(): pass\n")
    with pytest.raises(ValueError, match="current input changed"):
        current.install_current_inputs(root, tmp_path / "view", files, docs)


def test_plugin_changes_only_two_explicit_input_roots_and_restores_them(tmp_path):
    from types import SimpleNamespace
    from tools import completion_contract_inputs as plugin
    calls = []
    marker = lambda: calls.append("original function")
    modules = tuple(SimpleNamespace(__name__=name, ROOT=tmp_path / str(i), load_contract=marker)
                    for i, name in enumerate(plugin.CONSUMERS))
    view = tmp_path / "view"
    original = plugin.bind_consumers(modules, view)
    for module in modules:
        assert module.load_contract is marker
        module.load_contract()
    assert calls == ["original function", "original function"]
    plugin.restore_consumers(modules, original, view)
    assert tuple(module.ROOT for module in modules) == original
    assert not hasattr(plugin, "pytest_collection_modifyitems")
    assert len(plugin.COLD_TESTS) == 3


def test_plugin_rejects_unlisted_consumer_or_changed_root(tmp_path):
    from types import SimpleNamespace
    from tools import completion_contract_inputs as plugin
    modules = tuple(SimpleNamespace(__name__=name, ROOT=tmp_path / str(i))
                    for i, name in enumerate(plugin.CONSUMERS))
    with pytest.raises(ValueError, match="unreviewed"):
        plugin.bind_consumers(modules + (SimpleNamespace(__name__="other", ROOT=tmp_path),), tmp_path)
    original = plugin.bind_consumers(modules, tmp_path)
    modules[0].ROOT = tmp_path / "unreviewed"
    with pytest.raises(ValueError, match="root changed"):
        plugin.restore_consumers(modules, original, tmp_path)


def test_inherited_pytest_selection_cannot_hide_tests(monkeypatch, tmp_path):
    monkeypatch.setenv("PYTEST_ADDOPTS", "-k never_execute_anything --ignore=model")
    assert current.environment(tmp_path)["PYTEST_ADDOPTS"] == ""


def test_runtime_skip_cannot_be_reported_as_complete_execution(tmp_path):
    path = tmp_path / "test_skipped.py"
    path.write_text("import pytest\ndef test_runtime_skip(): pytest.skip('runtime skip')\n")
    output = tmp_path / "execution.json"
    environment = current.environment(current.ROOT)
    result = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "tools.check_completion_current_code",
                             "--oph-execution-report", str(output), path.name],
                            cwd=tmp_path, env=environment, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    with pytest.raises(ValueError, match="did not execute and pass"):
        current.verify_execution_report(output, ("test_skipped.py::test_runtime_skip",))


def test_deselection_cannot_be_reported_as_complete_execution(tmp_path):
    import json
    path = tmp_path / "execution.json"
    path.write_text(json.dumps({"collected": [], "reports": {}, "exitstatus": 0}))
    with pytest.raises(ValueError, match="inventory"):
        current.verify_execution_report(path, ("test_probe.py::test_current_code",))


def test_external_pytest_plugin_cannot_change_the_test_inventory(monkeypatch, tmp_path):
    from tools import completion_contract_inputs as plugin
    monkeypatch.setenv("PYTEST_PLUGINS", "external_selector")
    assert current.environment(tmp_path)["PYTEST_PLUGINS"] == ""
    with pytest.raises(ValueError, match="inherited pytest options"):
        plugin.validate_environment()


def test_cold_module_scope_rejects_non_model_entrypoint(tmp_path):
    from tools import finite_contract_command as command
    with pytest.raises(ValueError, match="explicit current model"):
        command.execute("external.module", [], tmp_path)


@pytest.fixture
def cold_view(tmp_path):
    import json
    import shutil
    root = current.ROOT
    family = json.loads((root / "model/completion/campaign.json").read_text())
    clock = json.loads((root / "model/clock_defect_energy/campaign.json").read_text())
    control = json.loads((root / "model/common_source/native_defect_decision_results.json").read_text())
    names = {"model/completion/campaign.json", "model/clock_defect_energy/campaign.json",
             "model/common_source/native_defect_decision_results.json"}
    names.update(row["path"] for group in ("source_bindings", "retained_controls") for row in family[group])
    names.update(clock["source_sha256"])
    names.update(control["source_sha256"])
    view = tmp_path / "cold-view"
    for name in names:
        target = view / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / name, target)
    for name, data in current.authenticated_documents(root).items():
        (view / name).write_bytes(data)
    return view


def test_cold_view_accepts_exact_documentary_inputs_with_current_science(cold_view):
    from tools import finite_contract_command as command
    command.authenticate_view(current.ROOT, cold_view)


def test_cold_view_rejects_a_changed_scientific_source(cold_view):
    from tools import finite_contract_command as command
    path = cold_view / "model/completion/observer.py"
    path.write_bytes(path.read_bytes() + b"\n# altered source\n")
    with pytest.raises(ValueError, match="scoped scientific or documentary input changed"):
        command.authenticate_view(current.ROOT, cold_view)


def test_cold_dispatch_accepts_only_three_original_node_ids(tmp_path):
    from tools import completion_contract_inputs as plugin
    assert plugin.COLD_TESTS == {
        "model/protected_defect_catalog/test_packet.py::test_fresh_full_producer_binds_runtime_imports",
        "model/clock_defect_energy/test_localized.py::test_retained_exact_receipt_and_mutation_rejection",
        "model/completion/test_defect_algebra.py::test_cli_rejects_scope_promotion_without_rewriting_input"}
    with pytest.raises(ValueError, match="unreviewed cold-test"):
        plugin.dispatch_cold_test("model/other/test_packet.py::test_result", tmp_path)


@pytest.mark.parametrize("body", ["assert False, 'original assertion failed'", "import pytest; pytest.skip('hidden skip')"])
def test_cold_dispatch_preserves_failure_and_rejects_skip(tmp_path, body):
    import shutil
    from tools import completion_contract_inputs as plugin
    node = "model/completion/test_defect_algebra.py::test_cli_rejects_scope_promotion_without_rewriting_input"
    view = tmp_path / "view"
    test = view / "model/completion/test_defect_algebra.py"
    test.parent.mkdir(parents=True)
    test.write_text("def test_cli_rejects_scope_promotion_without_rewriting_input():\n    " + body + "\n")
    (view / "pytest.ini").write_text("[pytest]\n")
    tool = view / "tools/check_completion_current_code.py"
    tool.parent.mkdir()
    shutil.copyfile(Path(current.__file__), tool)
    (tool.parent / "__init__.py").write_text("")
    if "skip" in body:
        with pytest.raises(ValueError, match="did not execute and pass"):
            plugin.dispatch_cold_test(node, view)
    else:
        reports = plugin.dispatch_cold_test(node, view)
        assert [report.outcome for report in reports] == ["passed", "failed", "passed"]
        assert "original assertion failed" in str(reports[1].longrepr)


def test_cold_report_cannot_forward_a_different_status_or_identity():
    from _pytest.reports import TestReport
    from tools import completion_contract_inputs as plugin
    node = next(iter(plugin.COLD_TESTS))
    report = TestReport(nodeid=node, location=(node.split("::")[0], 0, "test"), keywords={},
                        outcome="passed", longrepr=None, when="call", duration=0)
    row = {"phase": "call", "outcome": "passed", "xfail": False, "pytest_report": report._to_json()}
    assert plugin.decode_reports([row], node)[0].passed
    for key, value in [("phase", "setup"), ("outcome", "failed"), ("xfail", True)]:
        altered = dict(row, **{key: value})
        with pytest.raises(ValueError, match="identity or status"):
            plugin.decode_reports([altered], node)
    with pytest.raises(ValueError, match="identity or status"):
        plugin.decode_reports([row], node + "_other")


def test_fixed_owner_command_uses_original_cold_python_in_authenticated_view(tmp_path, monkeypatch):
    from tools import finite_contract_command as command
    calls = []
    class ReplacedProcess(Exception):
        pass
    monkeypatch.setattr(command, "authenticate_view", lambda root, view: calls.append(("authenticate", root, view)))
    monkeypatch.setattr(command.os, "chdir", lambda path: calls.append(("cwd", path)))
    def execute(executable, argv, environment):
        calls.append(("execve", executable, argv, environment))
        raise ReplacedProcess
    monkeypatch.setattr(command.os, "execve", execute)
    arguments = ["native_fair_local_readout", "--check", "model/native_fair_local_readout/results.json"]
    with pytest.raises(ReplacedProcess):
        command.execute("model.replay_native_fixed", arguments, tmp_path, root=tmp_path / "actual")
    assert calls[0] == ("authenticate", tmp_path / "actual", tmp_path)
    assert calls[1] == ("cwd", tmp_path.resolve())
    assert calls[2][1:3] == (command.sys.executable, [command.sys.executable, "-m", "model.replay_native_fixed", *arguments])
    assert calls[2][3]["PYTHONPATH"] == str(tmp_path.resolve())


@pytest.mark.parametrize("raising_finalizer", [False, True])
def test_cold_dispatch_preserves_parent_fixture_transition(tmp_path, raising_finalizer):
    """Exercise ordinary A, cold A, ordinary B in a real pytest process."""
    import json
    import shutil
    from tools import completion_contract_inputs as plugin

    view = tmp_path / "view"
    tools = view / "tools"
    tools.mkdir(parents=True)
    (tools / "__init__.py").write_text("")
    for module in (current, plugin):
        shutil.copyfile(Path(module.__file__), tools / Path(module.__file__).name)
    (view / "pytest.ini").write_text("[pytest]\n")
    (view / "conftest.py").write_text('''
import os
from pathlib import Path
import pytest
from tools import completion_contract_inputs as original

def pytest_configure(config):
    os.environ.setdefault("OPH_LIFECYCLE_FIXTURE_PARENT_PID", str(os.getpid()))
    config._oph_finite_view = Path(__file__).parent

@pytest.hookimpl(tryfirst=True)
def pytest_runtest_protocol(item, nextitem):
    if os.getpid() == int(os.environ["OPH_LIFECYCLE_FIXTURE_PARENT_PID"]):
        return original.pytest_runtest_protocol(item, nextitem)
''')
    first = view / "model/completion/test_defect_algebra.py"
    first.parent.mkdir(parents=True)
    first.write_text('''
import json
import os
from pathlib import Path
import pytest

def emit(event):
    with Path(os.environ["OPH_LIFECYCLE_FIXTURE_EVENTS"]).open("a") as stream:
        stream.write(json.dumps([os.getpid(), event]) + "\\n")

@pytest.fixture(scope="module", autouse=True)
def module_scope():
    emit("setup")
    yield
    emit("teardown")
    if (os.getpid() == int(os.environ["OPH_LIFECYCLE_FIXTURE_PARENT_PID"])
            and os.environ["OPH_LIFECYCLE_FIXTURE_RAISE"] == "1"):
        raise RuntimeError("original parent module finalizer failed")

def test_parent_before_cold():
    emit("ordinary")

def test_cli_rejects_scope_promotion_without_rewriting_input():
    assert os.getpid() != int(os.environ["OPH_LIFECYCLE_FIXTURE_PARENT_PID"])
    emit("cold")
''')
    second = first.with_name("test_field_current_measurement.py")
    second.write_text('''
import json
import os
from pathlib import Path

def test_next_module():
    events = Path(os.environ["OPH_LIFECYCLE_FIXTURE_EVENTS"])
    rows = [json.loads(line) for line in events.read_text().splitlines()]
    parent = int(os.environ["OPH_LIFECYCLE_FIXTURE_PARENT_PID"])
    assert [event for pid, event in rows if pid == parent] == ["setup", "ordinary", "teardown"]
    with events.open("a") as stream:
        stream.write(json.dumps([os.getpid(), "next"]) + "\\n")
''')
    report = tmp_path / "parent-report.json"
    events = tmp_path / "events.jsonl"
    environment = current.environment(view)
    environment.pop("OPH_LIFECYCLE_FIXTURE_PARENT_PID", None)
    environment["OPH_LIFECYCLE_FIXTURE_EVENTS"] = str(events)
    environment["OPH_LIFECYCLE_FIXTURE_RAISE"] = str(int(raising_finalizer))
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-W", "error", "-p",
         "tools.check_completion_current_code", "--oph-execution-report", str(report),
         str(first.relative_to(view)), str(second.relative_to(view))],
        cwd=view, env=environment, capture_output=True, text=True)
    assert result.returncode == int(raising_finalizer), result.stdout + result.stderr
    nodes = (
        "model/completion/test_defect_algebra.py::test_parent_before_cold",
        "model/completion/test_defect_algebra.py::test_cli_rejects_scope_promotion_without_rewriting_input",
        "model/completion/test_field_current_measurement.py::test_next_module",
    )
    value = json.loads(report.read_text())
    assert value["collected"] == list(nodes)
    if raising_finalizer:
        phases = value["reports"][nodes[1]]
        assert [(row["phase"], row["outcome"]) for row in phases] == [
            ("setup", "passed"), ("call", "passed"), ("teardown", "passed"),
            ("teardown", "failed")]
        assert "original parent module finalizer failed" in result.stdout
        assert [(row["phase"], row["outcome"]) for row in value["reports"][nodes[2]]] == [
            (phase, "passed") for phase in ("setup", "call", "teardown")]
    else:
        current.verify_execution_report(report, nodes)
    by_process = {}
    for pid, event in map(json.loads, events.read_text().splitlines()):
        by_process.setdefault(pid, []).append(event)
    assert sorted(by_process.values()) == sorted([
        ["setup", "ordinary", "teardown", "next"], ["setup", "cold", "teardown"]])
