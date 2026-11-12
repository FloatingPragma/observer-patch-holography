"""Supply two frozen documentary inputs to their original finite consumers.

Load explicitly with ``-p tools.completion_contract_inputs``. Current code and
tests retain their exact bytes. Only the filesystem roots of the two
named finite-contract readers use an authenticated copy. All other readers,
including current-hypothesis and live-receipt checks, retain current inputs.
No scientific function, receipt, assertion, or test-selection hook is replaced.
Three named tests with cold descendants execute in the authenticated file view;
their original node IDs and complete phase reports are retained.
"""
from pathlib import Path
from contextlib import contextmanager
import os
import json
import subprocess
import sys
import tempfile

import pytest
from tools import check_completion_current_code as views

CONSUMERS = ("model.completion.family", "model.clock_defect_energy.contract")
_REPORTS = {}
COLD_TESTS = frozenset((
    "model/protected_defect_catalog/test_packet.py::test_fresh_full_producer_binds_runtime_imports",
    "model/clock_defect_energy/test_localized.py::test_retained_exact_receipt_and_mutation_rejection",
    "model/completion/test_defect_algebra.py::test_cli_rejects_scope_promotion_without_rewriting_input",
))


def validate_environment():
    views.require(not os.environ.get("PYTEST_ADDOPTS", "") and not os.environ.get("PYTEST_PLUGINS", ""),
                  "inherited pytest options are not part of the finite-contract test command")


def pytest_runtest_logreport(report):
    if report.nodeid.startswith(("model/completion/", "model/clock_defect_energy/", "tools/test_completion_program.py::")):
        _REPORTS.setdefault(report.nodeid, []).append(
            (report.when, report.outcome, bool(getattr(report, "wasxfail", False))))


def bind_consumers(modules, view):
    views.require(tuple(module.__name__ for module in modules) == CONSUMERS,
                  "unreviewed finite-contract consumer")
    originals = tuple(module.ROOT for module in modules)
    for module in modules:
        module.ROOT = view
    return originals


def restore_consumers(modules, originals, view):
    views.require(all(module.ROOT == view for module in modules),
                  "finite-contract input root changed during tests")
    for module, original in zip(modules, originals, strict=True):
        module.ROOT = original


@contextmanager
def finite_inputs(root):
    """Execute current producers with only two named filesystem roots scoped."""
    root = Path(root).resolve()
    files = views.current_files(root)
    documents = views.authenticated_documents(root)
    live = subprocess.run([sys.executable, "tools/renew_native_defect_live.py", "--check"],
                          cwd=root, env=views.environment(root), capture_output=True, text=True)
    views.require(live.returncode == 0, "complete live native producer comparison failed:\n"
                  + live.stdout + live.stderr)
    temporary = tempfile.TemporaryDirectory(prefix="oph-finite-contract-inputs-")
    view = Path(temporary.name) / "view"
    modules, originals = (), ()
    try:
        views.install_current_inputs(root, view, files, documents)
        from model.completion import family
        from model.clock_defect_energy import contract
        modules = (family, contract)
        originals = bind_consumers(modules, view)
        yield view
        views.verify_view(view, files)
        views.require(views.current_files(root) == files,
                      "current executable or input changed during finite-contract execution")
        restore_consumers(modules, originals, view)
    finally:
        for module, original in zip(modules, originals, strict=True):
            module.ROOT = original
        temporary.cleanup()


def pytest_sessionstart(session):
    validate_environment()
    context = finite_inputs(session.config.rootpath)
    view = context.__enter__()
    session.config._oph_finite_view = view
    session.config._oph_finite_inputs = context
    reporter = session.config.pluginmanager.getplugin("terminalreporter")
    if reporter is not None:
        reporter.write_line("Finite contracts use two authenticated documentary inputs; all current code and tests execute. "
                            "This does not validate or identify the current hypothesis prose with the frozen premise.")


def dispatch_cold_test(nodeid, view):
    """Run one exact original node in a cold process, including its fixtures."""
    views.require(nodeid in COLD_TESTS, "unreviewed cold-test dispatch")
    view = Path(view)
    report_path = view.parent / ("test-" + views.sha(nodeid.encode()) + ".json")
    command = [sys.executable, "-m", "pytest", "-q", "-W", "error",
               "-p", "tools.check_completion_current_code", "--oph-execution-report", str(report_path), nodeid]
    completed = subprocess.run(command, cwd=view, env=views.environment(view), capture_output=True, text=True)
    views.require(report_path.is_file(), "cold test did not produce execution evidence:\n"
                  + completed.stdout + completed.stderr)
    value = json.loads(report_path.read_text())
    views.require(value["collected"] == [nodeid] and set(value["reports"]) == {nodeid},
                  "cold test execution inventory changed")
    if completed.returncode == 0:
        views.verify_execution_report(report_path, (nodeid,))
    else:
        views.require(any(row["outcome"] == "failed" for row in value["reports"][nodeid]),
                      "cold test process failed outside its original assertions:\n" + completed.stdout + completed.stderr)
    return decode_reports(value["reports"][nodeid], nodeid)


def decode_reports(rows, nodeid):
    from _pytest.reports import TestReport
    reports = []
    for row in rows:
        report = TestReport._from_json(row["pytest_report"])
        views.require((report.nodeid, report.when, report.outcome, bool(getattr(report, "wasxfail", False)))
                      == (nodeid, row["phase"], row["outcome"], row["xfail"]),
                      "cold phase report identity or status changed")
        reports.append(report)
    return reports


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_protocol(item, nextitem):
    if item.nodeid not in COLD_TESTS:
        return None
    item.ihook.pytest_runtest_logstart(nodeid=item.nodeid, location=item.location)
    try:
        reports = dispatch_cold_test(item.nodeid, item.config._oph_finite_view)
    except Exception as error:
        from _pytest.reports import TestReport
        reports = [TestReport(nodeid=item.nodeid, location=item.location, keywords={},
                              outcome="failed", longrepr=str(error), when="call", duration=0)]
    for report in reports:
        item.ihook.pytest_runtest_logreport(report=report)
    # The preceding ordinary item can retain collectors shared with this
    # item. Its fixtures still need the ordinary transition to the next item,
    # even though this item's own setup/call/teardown ran in the child.
    # Do not invoke the whole parent teardown hook: setup-dependent hooks
    # (including pytest's logging stash) were never entered for this item.
    from _pytest.outcomes import Exit
    from _pytest.runner import CallInfo
    if item.session.shouldfail or item.session.shouldstop:
        nextitem = None
    boundary = CallInfo.from_call(
        lambda: item.session._setupstate.teardown_exact(nextitem),
        when="teardown", reraise=(Exit,) if item.config.getoption("usepdb", False)
        else (Exit, KeyboardInterrupt))
    if boundary.excinfo is not None:
        report = item.ihook.pytest_runtest_makereport(item=item, call=boundary)
        item.ihook.pytest_runtest_logreport(report=report)
    item.ihook.pytest_runtest_logfinish(nodeid=item.nodeid, location=item.location)
    return True


@pytest.hookimpl(trylast=True)
def pytest_sessionfinish(session, exitstatus):
    context = getattr(session.config, "_oph_finite_inputs", None)
    if context is None:
        return
    try:
        if int(exitstatus) == 0 and not session.config.getoption("collectonly", default=False):
            expected = [item.nodeid for item in session.items if item.nodeid.startswith(
                ("model/completion/", "model/clock_defect_energy/", "tools/test_completion_program.py::"))]
            views.require(set(_REPORTS) == set(expected), "finite-contract test execution inventory changed")
            views.require(all(_REPORTS[name] == [(phase, "passed", False) for phase in ("setup", "call", "teardown")]
                              for name in expected), "a finite-contract test was skipped or did not pass all phases")
    finally:
        try:
            context.__exit__(None, None, None)
        finally:
            del session.config._oph_finite_inputs
            del session.config._oph_finite_view
