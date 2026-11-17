"""Partition the complete pytest collection by module, preserving every node.

Load explicitly with ``-p model.pytest_shards``. Without both shard options,
the ordinary full-suite command remains unchanged. All jobs collect the same
suite before selecting modules, so a collection error fails every job.
Recorded module costs only balance the shards; they never select or skip a
module.
"""
import json
from pathlib import Path
from statistics import median

import pytest

COSTS = Path(__file__).with_name("pytest_shard_costs.json")


def recorded_costs():
    try:
        seconds = json.loads(COSTS.read_text(encoding="utf-8"))["seconds"]
    except (OSError, ValueError, KeyError):
        return {}
    return {path: float(cost) for path, cost in seconds.items() if float(cost) > 0}


def selected_modules(modules, index, count, costs=None):
    """A disjoint, exhaustive partition of sorted repository-relative paths.

    Modules are placed longest first on the least-loaded shard. A module
    without a recorded cost takes the median recorded cost; equal costs keep
    sorted path order, so without costs shard i receives the i-th module first.
    """
    if (type(index) is not int or type(count) is not int
            or count < 1 or not 0 <= index < count):
        raise ValueError("require 0 <= shard index < positive shard count")
    ordered = sorted(set(modules))
    if count > len(ordered):
        raise ValueError("shard count exceeds collected test modules")
    costs = recorded_costs() if costs is None else costs
    default = median(costs.values()) if costs else 1.0
    cost = {path: costs.get(path, default) for path in ordered}
    loads, owner = [0.0] * count, {}
    for path in sorted(ordered, key=lambda path: -cost[path]):
        shard = min(range(count), key=lambda shard: loads[shard])
        loads[shard] += cost[path]
        owner[path] = shard
    return frozenset(path for path in ordered if owner[path] == index)


def pytest_addoption(parser):
    group = parser.getgroup("model regression partitions")
    group.addoption("--model-shard-index", type=int, default=None)
    group.addoption("--model-shard-count", type=int, default=None)


def pytest_configure(config):
    index = config.getoption("model_shard_index")
    count = config.getoption("model_shard_count")
    if index is None and count is None:
        return
    if (index is None or count is None or count < 1 or not 0 <= index < count):
        raise pytest.UsageError("supply both model shard options with 0 <= index < count")
    if config.getoption("continue_on_collection_errors", default=False):
        raise pytest.UsageError("model shards require collection errors to fail closed")


@pytest.hookimpl(trylast=True)
def pytest_collection_modifyitems(config, items):
    index = config.getoption("model_shard_index")
    count = config.getoption("model_shard_count")
    if index is None:
        return
    root = Path(config.rootpath)
    keys = [Path(item.path).relative_to(root).as_posix() for item in items]
    try:
        selected = selected_modules(keys, index, count)
    except ValueError as error:
        raise pytest.UsageError(str(error)) from error
    kept, removed = [], []
    for key, item in zip(keys, items, strict=True):
        (kept if key in selected else removed).append(item)
    config.hook.pytest_deselected(items=removed)
    items[:] = kept
    reporter = config.pluginmanager.getplugin("terminalreporter")
    if reporter is not None:
        reporter.write_line(
            f"model shard {index}/{count}: {len(kept)}/{len(keys)} collected tests, "
            f"{len(selected)}/{len(set(keys))} complete modules")
