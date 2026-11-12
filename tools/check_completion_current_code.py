#!/usr/bin/env python3
"""Test current completion code with its two authenticated documentary inputs.

The finite campaign identifies documentary inputs by their frozen bytes. The
current hypothesis is neither substituted publicly nor claimed equivalent.
Every executable, test, receipt, and configuration is copied from the current
worktree; only the two named documentary files use their exact Git blobs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
TESTS = ("tools/test_completion_program.py", "model/completion", "model/clock_defect_energy")
DOCUMENTS = {
    "docs/OPH_EXECUTABLE_HYPOTHESIS.md": {
        "commit": "9de91f9f0e01cb8fb185886a4e79579bb8b5a7a8",
        "sha256": "c503247a2b333c2e27b3eba865d47f2d15c966a5052eceb4aadb14fc9002f35a",
        "consumer": "model/completion/campaign.json",
        "role": "frozen finite campaign documentary premise; no current hypothesis equivalence",
    },
    "docs/research/OPH_FIELD_RECONSTRUCTION.md": {
        "commit": "6d39bfd75dda28a7affb199f7cb5f70d39feedf5",
        "sha256": "c2aeed6e6b716785b96c88b06628232b2c5d8f06068925f12a2c7968ff33f1df",
        "consumer": "model/common_source/native_defect_decision_results.json",
        "role": "retained native-control documentary input; current producer checked separately",
    },
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def current_files(root):
    root = Path(root).resolve()
    names = subprocess.check_output(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"], cwd=root
    ).decode().split("\0")
    result = {}
    for name in sorted(set(names) - {""}):
        path = root / name
        require(path.resolve().is_relative_to(root), "escaping current input: " + name)
        if path.is_file():
            result[name] = sha(path.read_bytes())
    return result


def authenticated_documents(root, *, documentary_paths=tuple(DOCUMENTS)):
    root = Path(root).resolve()
    require(set(documentary_paths) <= set(DOCUMENTS), "unknown documentary input requested")
    result = {}
    for name in documentary_paths:
        spec = DOCUMENTS[name]
        commit = subprocess.check_output(
            ["git", "rev-parse", spec["commit"] + "^{commit}"], cwd=root, text=True).strip()
        require(commit == spec["commit"], "documentary commit identity changed")
        ancestry = subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"], cwd=root)
        require(ancestry.returncode == 0, "documentary commit is not an ancestor of current HEAD")
        data = subprocess.check_output(["git", "show", commit + ":" + name], cwd=root)
        require(sha(data) == spec["sha256"], "documentary Git blob changed: " + name)
        consumer = json.loads((root / spec["consumer"]).read_text())
        if name == "docs/OPH_EXECUTABLE_HYPOTHESIS.md":
            refs = [row for row in consumer["retained_controls"] if row["path"] == name]
            require(len(refs) == 1 and refs[0]["sha256"] == spec["sha256"],
                    "finite campaign documentary input changed")
        else:
            require(consumer["source_sha256"][name] == spec["sha256"],
                    "native-control documentary input changed")
        result[name] = data
    return result


def install_current_inputs(root, destination, files, documents, *, documentary_paths=tuple(DOCUMENTS)):
    """Copy actual current files; never replace code or scientific receipts."""
    root, destination = Path(root).resolve(), Path(destination).resolve()
    require(set(documentary_paths) <= set(DOCUMENTS) and set(documents) == set(documentary_paths),
            "documentary replacement inventory changed")
    for name, data in documents.items():
        require(sha(data) == DOCUMENTS[name]["sha256"], "unauthenticated documentary bytes")
    # A shared Git clone retains immutable objects needed by unchanged custody
    # tests. All worktree files below are then overwritten from current inputs.
    subprocess.run(["git", "clone", "--quiet", "--shared", "--no-checkout", str(root), str(destination)],
                   check=True, capture_output=True, text=True)
    subprocess.run(["git", "read-tree", "HEAD"], cwd=destination, check=True, capture_output=True)
    for name, expected in files.items():
        source, target = root / name, destination / name
        require(sha(source.read_bytes()) == expected, "current input changed during copy: " + name)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    for name, data in documents.items():
        require(name in files, "documentary file absent from current tree")
        (destination / name).write_bytes(data)
    verify_view(destination, files, documentary_paths=documentary_paths)


def verify_view(view, files, *, documentary_paths=tuple(DOCUMENTS)):
    view = Path(view)
    require(set(documentary_paths) <= set(DOCUMENTS), "unknown documentary input requested")
    require(set(current_files(view)) == set(files), "private current input inventory changed")
    for name, expected in files.items():
        want = DOCUMENTS[name]["sha256"] if name in documentary_paths else expected
        require((view / name).is_file() and sha((view / name).read_bytes()) == want,
                "current executable, input, or documentary view changed: " + name)


def environment(view):
    return dict(os.environ, PYTHONPATH=str(view), PYTHONDONTWRITEBYTECODE="1",
                PYTEST_DISABLE_PLUGIN_AUTOLOAD="1", PYTEST_ADDOPTS="", PYTEST_PLUGINS="", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")


def collect(view, tests=TESTS):
    command = [sys.executable, "-m", "pytest", "--collect-only", "-q", *tests]
    result = subprocess.run(command, cwd=view, env=environment(view), capture_output=True, text=True)
    require(result.returncode == 0, "current test collection failed:\n" + result.stdout + result.stderr)
    ids = tuple(line.strip() for line in result.stdout.splitlines() if "::" in line)
    require(ids and len(ids) == len(set(ids)), "empty or duplicate current test inventory")
    return ids



# This optional reporting plugin records actual execution independently of the
# process exit status. It has no collection, selection, or runtest override.
_EXECUTION_REPORTS = {}


def pytest_addoption(parser):
    parser.addoption("--oph-execution-report", default=None)


def pytest_runtest_logreport(report):
    _EXECUTION_REPORTS.setdefault(report.nodeid, []).append(
        {"phase": report.when, "outcome": report.outcome,
         "xfail": bool(getattr(report, "wasxfail", False)), "pytest_report": report._to_json()})


def pytest_sessionfinish(session, exitstatus):
    output = session.config.getoption("oph_execution_report", default=None)
    if output:
        Path(output).write_text(json.dumps({"collected": [item.nodeid for item in session.items],
            "reports": _EXECUTION_REPORTS, "exitstatus": int(exitstatus)}, sort_keys=True) + "\n")


def verify_execution_report(path, expected):
    value = json.loads(Path(path).read_text())
    require(value["collected"] == list(expected) and set(value["reports"]) == set(expected),
            "executed test inventory differs from complete collection")
    require(value["exitstatus"] == 0, "test execution failed")
    for name in expected:
        phases = value["reports"][name]
        require([row["phase"] for row in phases] == ["setup", "call", "teardown"]
                and all(row["outcome"] == "passed" and row["xfail"] is False for row in phases),
                "test did not execute and pass all phases: " + name)


def execute(*, root=ROOT):
    root = Path(root).resolve()
    before = current_files(root)
    documents = authenticated_documents(root)
    original_ids = collect(root)
    # Current prose remains current. Its scientific producer runs against that
    # prose and its entire non-provenance result must match retained evidence.
    live = subprocess.run([sys.executable, "tools/renew_native_defect_live.py", "--check"],
                          cwd=root, env=environment(root), capture_output=True, text=True)
    require(live.returncode == 0, "current native producer replay failed:\n" + live.stdout + live.stderr)
    with tempfile.TemporaryDirectory(prefix="oph-completion-current-code-") as folder:
        view = Path(folder) / "view"
        install_current_inputs(root, view, before, documents)
        copied_ids = collect(view)
        require(copied_ids == original_ids, "current test inventory changed in documentary view")
        report_path = Path(folder) / "execution.json"
        result = subprocess.run([sys.executable, "-m", "pytest", "-q", "--durations=10", "-W", "error",
                                 "-p", "tools.check_completion_current_code",
                                 "--oph-execution-report", str(report_path), *TESTS],
                                cwd=view, env=environment(view), capture_output=True, text=True)
        print(result.stdout, end="", flush=True)
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="", flush=True)
        require(result.returncode == 0, "current completion executable regression failed")
        verify_execution_report(report_path, original_ids)
        verify_view(view, before)
    require(current_files(root) == before, "current worktree changed during complete regression")
    return {
        "schema": "oph.current_completion_code_frozen_documentary_inputs.v1",
        "scope": "current executable code and complete tests on the retained finite documentary contract",
        "current_input_file_count": len(before), "test_count": len(original_ids),
        "test_inventory_sha256": sha("\n".join(original_ids).encode()),
        "current_worktree_manifest_sha256": sha(json.dumps(before, sort_keys=True).encode()),
        "documentary_inputs": DOCUMENTS, "current_hypothesis_prose_validated": False,
        "current_and_frozen_hypotheses_claimed_equivalent": False,
        "current_executable_test_and_scientific_input_bytes_preserved": True,
        "scientific_assertions_removed": False,
        "current_native_producer_replayed_and_complete_scientific_payload_compared": True,
        "historical_sources_or_public_documents_rewritten": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = execute()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
