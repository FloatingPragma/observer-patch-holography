"""Partition the actual pytest collection without changing its acceptance gates."""
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from .pytest_shards import recorded_costs, selected_modules


@pytest.mark.parametrize("count", range(1, 8))
def test_modules_are_exhaustive_disjoint_with_all_parametrizations(count):
    modules = [f"model/part{i}/test_example.py" for i in range(7)]
    nodes = [(module, f"case{case}") for module in modules for case in range(1, 5)]
    parts = [selected_modules(modules, index, count) for index in range(count)]
    assert Counter(m for part in parts for m in part) == Counter(modules)
    assert Counter(node for part in parts for node in nodes if node[0] in part) == Counter(nodes)
    assert all(parts)


def test_recorded_costs_balance_shards_without_changing_membership():
    modules = [f"model/part{i}/test_example.py" for i in range(6)]
    costs = dict(zip(modules, (50.0, 1.0, 1.0, 1.0, 1.0, 46.0)))
    parts = [selected_modules(modules, index, 2, costs) for index in range(2)]
    assert Counter(m for part in parts for m in part) == Counter(modules)
    assert sorted(sum(costs[m] for m in part) for part in parts) == [50.0, 50.0]
    assert recorded_costs() and all(path.startswith("model/") for path in recorded_costs())


def test_new_modules_are_included_without_a_second_allowlist():
    modules = ["model/a/test_one.py", "model/b/test_two.py", "model/new/test_three.py"]
    assert set.union(*(set(selected_modules(modules, i, 3)) for i in range(3))) == set(modules)
    assert selected_modules(reversed(modules), 0, 3) == {"model/a/test_one.py"}


@pytest.mark.parametrize("index,count", [(-1, 3), (3, 3), (0, 0), (0, 4),
                                        (True, 3), (0, True), (0.0, 3), (0, 3.0)])
def test_invalid_and_empty_partitions_fail(index, count):
    with pytest.raises(ValueError):
        selected_modules(["a", "b", "c"], index, count)
    with pytest.raises(ValueError):
        selected_modules([], 0, 1)


@pytest.fixture
def suite(tmp_path):
    (tmp_path / "pytest.ini").write_text("[pytest]\npython_files = test_*.py\n", encoding="utf-8")
    (tmp_path / "conftest.py").write_text('''
import json
from pathlib import Path
seen = []
deselected = []
def pytest_collection_modifyitems(items):
    seen.extend(item.nodeid for item in items)
def pytest_deselected(items):
    deselected.extend(item.nodeid for item in items)
def pytest_sessionfinish(session):
    Path("collection.json").write_text(json.dumps({
        "all": seen, "deselected": deselected,
        "selected": [item.nodeid for item in session.items]}), encoding="utf-8")
''', encoding="utf-8")
    model = tmp_path / "model"
    model.mkdir()
    for name in ("a", "b", "c"):
        (model / f"test_{name}.py").write_text('''
import pytest
@pytest.fixture(scope="module")
def value(): return 7
@pytest.mark.parametrize("n", [1, 2, 3])
def test_parameter(value, n): assert value == 7 and n > 0
@pytest.mark.skip(reason="declared test scope")
def test_declared_skip(): pass
''', encoding="utf-8")
    return tmp_path


def run_suite(suite, *options):
    env = os.environ.copy()
    root = str(Path(__file__).resolve().parents[1])
    env["PYTHONPATH"] = root + os.pathsep + env.get("PYTHONPATH", "")
    env.pop("PYTEST_ADDOPTS", None)
    env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
    result = subprocess.run([sys.executable, "-m", "pytest", "-q", "-W", "error",
                             "-p", "model.pytest_shards", *options, "model"],
                            cwd=suite, env=env, capture_output=True, text=True, timeout=45)
    record = suite / "collection.json"
    data = json.loads(record.read_text(encoding="utf-8")) if record.exists() else None
    record.unlink(missing_ok=True)
    return result, data


def test_actual_full_collection_and_three_shards_preserve_every_test_and_skip(suite):
    baseline, whole = run_suite(suite)
    assert baseline.returncode == 0, baseline.stdout + baseline.stderr
    assert "9 passed, 3 skipped" in baseline.stdout
    selected = []
    for index in range(3):
        result, part = run_suite(suite, "--model-shard-index", str(index), "--model-shard-count", "3")
        assert result.returncode == 0, result.stdout + result.stderr
        assert Counter(part["all"]) == Counter(whole["all"])
        assert Counter(part["selected"] + part["deselected"]) == Counter(whole["all"])
        assert "3 passed, 1 skipped" in result.stdout
        selected.extend(part["selected"])
    assert Counter(selected) == Counter(whole["selected"])


@pytest.mark.parametrize("options", [
    ["--model-shard-index", "0"], ["--model-shard-count", "3"],
    ["--model-shard-index", "3", "--model-shard-count", "3"],
    ["--model-shard-index", "0", "--model-shard-count", "4"],
    ["--model-shard-index", "0", "--model-shard-count", "3", "--continue-on-collection-errors"],
])
def test_invalid_cli_cannot_run_a_partial_acceptance_suite(suite, options):
    result, _ = run_suite(suite, *options)
    assert result.returncode == pytest.ExitCode.USAGE_ERROR


@pytest.mark.parametrize("index", range(3))
def test_collection_failure_in_unselected_module_fails_every_shard(suite, index):
    (suite / "model/test_broken.py").write_text("raise RuntimeError('broken collection')\n", encoding="utf-8")
    result, _ = run_suite(suite, "--model-shard-index", str(index), "--model-shard-count", "3")
    assert result.returncode == pytest.ExitCode.INTERRUPTED
    assert "broken collection" in result.stdout


def test_selected_failure_is_not_lost_or_converted_to_a_pass(suite):
    (suite / "model/test_a.py").write_text("def test_failure(): assert False\n", encoding="utf-8")
    result, data = run_suite(suite, "--model-shard-index", "0", "--model-shard-count", "3")
    assert result.returncode == pytest.ExitCode.TESTS_FAILED
    assert data["selected"] == ["model/test_a.py::test_failure"]
