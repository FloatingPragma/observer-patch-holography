"""Run exactly the unchanged particle-interface tests with guarded metadata reuse.

This is an adapted pytest cold replay, not an unmodified CLI replay. Only the
four reviewed path inventories are cached; the original tests, scientific
report and strict receipt comparison execute unchanged. Use an immutable,
isolated checkout with a current receipt and an external audit destination.

This isolated module uses only pytest's built-in tmp_path/monkeypatch fixtures
and its own parametrization. It has no repository conftest or plugin hooks.
That reviewed harness is pinned below; it is not a routing policy for the
general model suite or its retained-documentary adapters.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import tempfile
import time

TEST = "model/particle_identification/test_completion_certificate.py"
TEST_SHA256 = "522b78c013be67a508aba3b48e5ee32fea8b0fc95b1120105fefd6147a1fad05"
GUARD_SHA256 = "19f030d5c246efeb8a4556154a3ac0789ea2974d1483813f52f8bfb6d9c2ad61"
PYTEST_INI_SHA256 = "4f40f83f8d1457725aba1ca61b9d9f251ec33e561316395dfb793c7ef72a7448"
ABSENT_CONFTESTS = ("conftest.py", "model/conftest.py", "model/particle_identification/conftest.py")
TEST_NAMES = (
    "test_complete_interface_receipt_replays",
    "test_receipt_equality_preserves_scalar_types[False]",
    "test_receipt_equality_preserves_scalar_types[0.0]",
    "test_receipt_equality_rejects_nonfinite_values",
    "test_fair_interface_receipts_are_freshly_computed",
    "test_fair_interface_rejects_changed_owner_receipt[float_integer]",
    "test_fair_interface_rejects_changed_owner_receipt[bool_integer]",
    "test_fair_interface_rejects_changed_owner_receipt[missing_number]",
    "test_fair_interface_rejects_changed_owner_receipt[missing_failure]",
    "test_fair_interface_rejects_changed_owner_receipt[stale_value]",
    "test_fair_interface_rejects_duplicate_receipt_keys",
    "test_nested_fair_replays_and_sources_are_registered_for_model_verify",
)
EXPECTED_IDS = frozenset(TEST + "::" + name for name in TEST_NAMES)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_guard():
    path = Path(__file__).with_name("immutable_custody_replay.py").resolve()
    if sha(path) != GUARD_SHA256:
        raise ValueError("guarded replay implementation requires renewed review")
    name = "_oph_immutable_custody_pytest_guard"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError("cannot load the reviewed custody implementation")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module, path


def authenticate_test(root):
    target = (root / TEST).resolve()
    if not target.is_relative_to(root) or sha(target) != TEST_SHA256:
        raise ValueError("original test module requires renewed review")
    return target


def authenticate_harness(root):
    if sha(root / "pytest.ini") != PYTEST_INI_SHA256:
        raise ValueError("repository pytest harness requires renewed review")
    for name in ABSENT_CONFTESTS:
        path = root / name
        if path.exists() or path.is_symlink():
            raise ValueError("isolated module gained an unreviewed conftest: " + name)


class ExactTests:
    """Reject deselection, duplicates, skips, xfails and incomplete execution."""
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.collected = []
        self.reports = {}
        self.errors = []

    def pytest_collection_finish(self, session):
        self.collected = [item.nodeid for item in session.items]
        if (len(self.collected) != len(EXPECTED_IDS)
                or set(self.collected) != EXPECTED_IDS
                or any(Path(item.path).resolve() != self.root / TEST for item in session.items)):
            self.errors.append("collection differs from the complete original 12 tests")
            raise ValueError(self.errors[-1])

    def pytest_collectreport(self, report):
        if report.outcome != "passed":
            self.errors.append("collection was not fully successful")

    def pytest_runtest_logreport(self, report):
        key = (report.nodeid, report.when)
        if key in self.reports:
            self.errors.append("duplicate test-phase result")
        self.reports[key] = report.outcome
        if (report.nodeid not in EXPECTED_IDS or report.when not in ("setup", "call", "teardown")
                or report.outcome != "passed" or hasattr(report, "wasxfail")):
            self.errors.append("unexpected, failed, skipped or xfailed test phase")

    def validate(self, exit_code):
        expected = {(node, phase) for node in EXPECTED_IDS for phase in ("setup", "call", "teardown")}
        if (int(exit_code) != 0 or self.errors or set(self.collected) != EXPECTED_IDS
                or len(self.collected) != 12 or set(self.reports) != expected
                or any(outcome != "passed" for outcome in self.reports.values())):
            raise ValueError("complete original 12-test execution was not accepted")

    def evidence(self):
        return {"collected_ids": self.collected,
                "phase_reports": [{"nodeid": key[0], "phase": key[1], "outcome": outcome}
                                  for key, outcome in sorted(self.reports.items())],
                "errors": self.errors}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args(argv)
    root, audit, adapter = args.root.resolve(), args.audit.resolve(), Path(__file__).resolve()
    guard, guard_path = load_guard()
    if (audit.is_relative_to(root) or guard.same_file(audit, adapter)
            or guard.same_file(audit, guard_path) or args.audit.is_symlink()
            or (audit.exists() and (not audit.is_file() or audit.stat().st_nlink != 1))):
        parser.error("audit must be an external unaliased report path")
    if any(os.environ.get(name) for name in ("PYTEST_ADDOPTS", "PYTEST_PLUGINS")):
        parser.error("PYTEST_ADDOPTS and PYTEST_PLUGINS must be empty for exact test routing")
    authenticate_test(root)
    authenticate_harness(root)
    old_cwd, old_path, old_argv = Path.cwd(), list(sys.path), list(sys.argv)
    old_autoload = os.environ.get("PYTEST_DISABLE_PLUGIN_AUTOLOAD")
    started = time.monotonic()
    snapshot = guard.FrozenSnapshot(root)
    slots, plugin = [], ExactTests(root)
    expected = root / guard.OUTPUT
    report = {"schema": "oph.adapted_metadata_pytest_replay.v1", "mode": "audited_path_inventory_cache_pytest",
              "root": str(root), "original_test": TEST, "original_test_sha256": TEST_SHA256,
              "adapter_sha256": sha(adapter), "guard_sha256": sha(guard_path),
              "input_sha256": snapshot.bindings, "expected_receipt_sha256_before": sha(expected),
              "python": sys.version, "platform": platform.platform(),
              "unmodified_cli_run": False, "adapted_pytest_run": True,
              "scientific_result_cache": False, "original_tests_modified": False,
              "external_plugin_autoload_disabled": True, "repository_conftest_disabled": True,
              "reviewed_repository_pytest_ini_sha256": PYTEST_INI_SHA256,
              "reviewed_absent_ancestor_conftests": list(ABSENT_CONFTESTS),
              "original_fixture_contract": "Only built-in tmp_path and monkeypatch plus in-module parametrization; no repository plugin fixtures."}
    try:
        os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
        os.chdir(root)
        sys.path.insert(0, str(root))
        if any(guard.same_file(audit, root / name) for name in snapshot.names):
            raise ValueError("audit aliases an immutable source")
        guard.authenticate_owner(root)
        guard.authenticate_loaded_modules(snapshot)
        guard.install(snapshot, frozenset({guard.OUTPUT}), slots)
        import pytest
        report["pytest_version"] = pytest.__version__
        with tempfile.TemporaryDirectory(prefix="oph-completion-pytest-") as temporary:
            temp = Path(temporary)
            config = temp / "pytest.ini"
            config.write_text("[pytest]\n")
            pytest_args = [TEST, "-q", "-W", "error", "-p", "no:cacheprovider", "--noconftest",
                           "--import-mode=importlib", "-c", str(config), "--rootdir", str(root),
                           "--basetemp", str(temp / "cases")]
            report["pytest_arguments"] = pytest_args
            code = pytest.main(pytest_args, plugins=[plugin])
            report["pytest_exit_code"] = int(code)
            plugin.validate(code)
        guard.check_expected(expected, report["expected_receipt_sha256_before"])
        for slot in slots:
            slot.verify()
        authenticate_test(root)
        authenticate_harness(root)
        guard.authenticate_loaded_modules(snapshot)
        snapshot.check()
        if sha(adapter) != report["adapter_sha256"] or sha(guard_path) != report["guard_sha256"]:
            raise ValueError("replay implementation changed during test execution")
        report.update(verdict="PASS_ADAPTED_COMPLETE_ORIGINAL_12_TEST_MODULE",
                      expected_receipt_sha256_after=sha(expected), tests_passed=12,
                      original_scientific_report_and_strict_comparison_executed=True,
                      immutable_source_bytes_and_filenames_preserved=True)
    except BaseException as error:
        report.update(verdict="FAIL", error=repr(error))
        raise
    finally:
        for slot in slots:
            slot.module.source_bindings = slot.original
        report["test_execution"] = plugin.evidence()
        report["metadata_caches"] = [slot.evidence() for slot in slots]
        report["elapsed_seconds"] = time.monotonic() - started
        if old_autoload is None:
            os.environ.pop("PYTEST_DISABLE_PLUGIN_AUTOLOAD", None)
        else:
            os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = old_autoload
        os.chdir(old_cwd)
        sys.path[:], sys.argv = old_path, old_argv
        if not any(guard.same_file(audit, root / name) for name in snapshot.names):
            audit.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(json.dumps({key: report.get(key) for key in ("verdict", "mode", "tests_passed", "elapsed_seconds")}), flush=True)


if __name__ == "__main__":
    main()
