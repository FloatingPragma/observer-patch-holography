"""Synthetic receipts exercise the maintenance boundary, never live evidence."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import threading
import time
import types

import pytest
import yaml

import replay_model_receipts as replay


def test_completion_maintenance_artifact_manifest():
    """Named artifact integrity is distinct from scientific replay admission."""
    root = Path(replay.__file__).resolve().parents[1]
    path = root / "model/completion/evidence_maintenance_manifest.json"
    manifest = json.loads(path.read_text())
    assert manifest["schema"] == "oph.evidence_maintenance_artifacts.v1"
    assert set(manifest["files"]) == {
        "model/completion/evidence_maintenance_review.json",
        "model/completion/evidence_maintenance_transitions.json",
        "model/completion/ordered_spin_energy_readout.json",
        "model/completion/ORDERED_SPIN_ENERGY_READOUT.md",
        "tools/read_ordered_spin_energies.py",
        "tools/test_read_ordered_spin_energies.py",
    }
    for name, expected in manifest["files"].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected


def test_completion_ci_cold_checks_every_new_package_and_propagates_failure(monkeypatch):
    import subprocess
    workflow = yaml.safe_load((Path(replay.__file__).resolve().parents[1] /
                              ".github/workflows/completion-program.yml").read_text())
    packages = ("spin_goldstone_sector", "native_spin_observer", "lattice_weyl_response",
                "cubic_particle_channels", "lattice_modular_limit", "spin_particle_scaling",
                "native_ordered_ground", "spin_shell_diagnostic")
    events = workflow.get("on", workflow.get(True))
    steps = workflow["jobs"]["native-field-thermodynamics"]["steps"]
    cold = next(s for s in steps if s.get("name") == "Cold replay continuum and ordered-spin receipts")
    tests = next(s for s in steps if s.get("name") == "Check continuum and ordered-spin observer evidence")
    for package in packages:
        assert "model/" + package in tests["run"]
        for event in ("push", "pull_request"):
            assert f"model/{package}/**" in events[event]["paths"]
    calls = []

    def checked(argv, *, check):
        assert check is True
        calls.append(argv)

    monkeypatch.setattr(subprocess, "run", checked)
    exec(compile(cold["run"], "<CI cold checks>", "exec"), {})
    assert [argv[-2:] for argv in calls] == [["model." + p, "--check"] for p in packages]

    def failed(argv, *, check):
        raise subprocess.CalledProcessError(1, argv)

    monkeypatch.setattr(subprocess, "run", failed)
    with pytest.raises(subprocess.CalledProcessError):
        exec(compile(cold["run"], "<CI fail-closed cold checks>", "exec"), {})


def test_dependency_cycle_is_not_silently_rehashed():
    with pytest.raises(ValueError, match="dependency cycle"):
        replay.topological({"first": {"second"}, "second": {"first"}})


def test_repeated_invocation_preserves_original_audit_bytes(tmp_path, monkeypatch):
    output = tmp_path / "audit"
    output.mkdir()
    original = output / "report.json"
    original.write_bytes(b"original audited report")
    monkeypatch.setattr(sys, "argv", ["replay", "--root", str(tmp_path),
                                     "--output", str(output), "--execute"])
    with pytest.raises(SystemExit):
        replay.main()
    assert original.read_bytes() == b"original audited report"


def test_payload_audit_preserves_types_scope_and_source_path_catalog(tmp_path):
    (tmp_path / "input").write_text("source")
    old = {"sources": {"input": "old"}, "theorem": {"proved": True, "count": 3}}
    refreshed = {"sources": {"input": "new"}, "theorem": {"proved": True, "count": 3}}
    assert replay.differences(replay.science_payload(tmp_path, "receipt.json", old),
                              replay.science_payload(tmp_path, "receipt.json", refreshed)) == []
    refreshed["theorem"]["proved"] = 1
    assert replay.differences(replay.science_payload(tmp_path, "receipt.json", old),
                              replay.science_payload(tmp_path, "receipt.json", refreshed))
    refreshed = {"sources": {}, "theorem": old["theorem"]}
    assert replay.differences(replay.science_payload(tmp_path, "receipt.json", old),
                              replay.science_payload(tmp_path, "receipt.json", refreshed))


@pytest.mark.parametrize("bad_path", ["../outside", "/tmp/outside"])
def test_local_source_map_rejects_escape(tmp_path, bad_path):
    with pytest.raises(ValueError, match="escaping|Absolute"):
        list(replay.bindings(tmp_path, "receipt.json", {"source_bindings": {bad_path: "hash"}}))


def test_archive_ownership_includes_chunks_but_not_imported_outputs(tmp_path):
    for name in ("a.json", "a.npz", "a.state-0.npz", "imported.npz"):
        (tmp_path / name).write_bytes(b"bytes")
    value = {"complete_record_archive": {"files": [{"filename": "a.npz"},
                                                    {"filename": "a.state-0.npz"}]},
             "prior": {"artifact": {"filename": "imported.npz"}}}
    owners, _ = replay.output_owners(tmp_path, {"a.json": value})
    assert set(owners) == {tmp_path / p for p in ("a.json", "a.npz", "a.state-0.npz")}
    with pytest.raises(ValueError, match="Multiple producers"):
        replay.output_owners(tmp_path, {"a.json": value,
                                      "b.json": {"artifact": {"filename": "a.npz"}}})
    with pytest.raises(ValueError, match="Escaping archive"):
        replay.owned_outputs(tmp_path, "a.json", {"artifact": {"filename": "../a.npz"}})


@pytest.mark.parametrize("failure", [None, "producer", "checker", "science"])
def test_parallel_replay_never_admits_failed_or_changed_premise(tmp_path, monkeypatch, failure):
    root = tmp_path.resolve()
    (root / "model").mkdir()
    (root / "model/verify.py").write_text("synthetic verifier")
    (root / "input").write_text("new source")
    sha = lambda path: hashlib.sha256((root / path).read_bytes()).hexdigest()
    paths = {name: f"model/{name}.json" for name in ("a", "b", "c")}
    for name in ("a", "b"):
        (root / paths[name]).write_text(json.dumps({"sources": {"input": "stale"}, "science": 13}))
    archive = root / "model/a.npz"
    archive.write_bytes(b"original array bytes")
    # c consumes a's array directly, without a reference to a's JSON receipt.
    (root / paths["c"]).write_text(json.dumps({"sources": {"model/a.npz": sha("model/a.npz"),
                                                              paths["b"]: sha(paths["b"])}, "science": 17}))
    frozen = root / "frozen.json"
    frozen.write_text('{"historical": "not a registered live receipt"}')
    frozen_bytes = frozen.read_bytes()

    def final_check():
        for path in paths.values():
            value = json.loads((root / path).read_text())
            assert all(expected == sha(source) for source, expected in value["sources"].items())

    verifier = types.SimpleNamespace(ROOT=root, REPLAYS=[(f"fake.{n}", None, p) for n, p in paths.items()],
                                     verify_receipt_source_bindings=final_check)
    monkeypatch.setattr(replay, "SUPPORTED_WRITERS", {
        (f"fake.{name}", None): ("model/verify.py", sha("model/verify.py"))
        for name in paths
    })
    fake_model = types.ModuleType("model")
    fake_model.verify = verifier
    fake_tools = types.ModuleType("tools")
    fake_tools.strict_json = types.SimpleNamespace(loads=json.loads)
    monkeypatch.setitem(sys.modules, "model", fake_model)
    monkeypatch.setitem(sys.modules, "tools", fake_tools)
    lock = threading.Lock()
    entered, checked, events = set(), set(), []
    both_roots = threading.Event()

    def fake_run(argv, **kwargs):
        operation = argv[-2]
        name = Path(argv[-1]).stem
        with lock:
            events.append((name, operation))
            if name == "c":
                assert checked == {"a", "b"}
            if operation == "--write" and name in ("a", "b"):
                entered.add(name)
                if len(entered) == 2:
                    both_roots.set()
        if operation == "--write" and name in ("a", "b"):
            assert both_roots.wait(3), "Independent producers should overlap"
            if name == "b":
                time.sleep(.05)
        if name == "a" and ((failure == "producer" and operation == "--write") or
                             (failure == "checker" and operation == "--check")):
            return types.SimpleNamespace(returncode=23)
        if operation == "--write":
            destination = Path(argv[-1])
            if name == "a":
                destination.with_suffix(".npz").write_bytes(b"recomputed array bytes")
            path = root / paths[name]
            value = json.loads(path.read_text())
            value["sources"] = {source: sha(source) for source in value["sources"]}
            if name == "a" and failure == "science":
                value["science"] = 14
            # Match exclusive-create producers, which refuse in-place writes.
            with destination.open("x") as stream:
                json.dump(value, stream)
        else:
            with lock:
                checked.add(name)
        return types.SimpleNamespace(returncode=0)

    monkeypatch.setattr(replay, "subprocess", types.SimpleNamespace(run=fake_run, STDOUT=-2))
    monkeypatch.setattr(sys, "argv", ["replay", "--root", str(root), "--output", str(root / "out"),
                                     "--execute", "--jobs", "2"])
    old_path = list(sys.path)
    try:
        if failure:
            with pytest.raises(SystemExit, match="Stopped"):
                replay.main()
        else:
            replay.main()
    finally:
        sys.path[:] = old_path
    report = json.loads((root / "out/report.json").read_text())
    assert frozen.read_bytes() == frozen_bytes
    assert (root / "out/before/model/a.json").is_file()
    assert (root / "out/before/model/b.json").is_file()
    assert (root / "out/before/model/a.npz").read_bytes() == b"original array bytes"
    if failure:
        assert not any(name == "c" for name, _ in events)
        assert checked == {"b"}
        assert report["running"] == []
        assert report["not_scheduled"] == [paths["c"]]
        assert report["failures"][0]["receipt"] == paths["a"]
        assert archive.read_bytes() == b"original array bytes"
    else:
        assert checked == {"a", "b", "c"}
        assert report["authoritative_final_source_check"] == "PASS"
        assert all(row["non_provenance_differences"] == [] for row in report["results"])


def test_pathless_store_true_writer_is_rejected_before_any_producer(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    (root / "model").mkdir()
    (root / "model/verify.py").write_text("synthetic verifier")
    (root / "input").write_text("changed source")
    marker = root / "CANONICAL_WAS_OVERWRITTEN"
    source = root / "pathless.py"
    source.write_text(
        "import argparse\nfrom pathlib import Path\n"
        "p = argparse.ArgumentParser()\np.add_argument('--write', action='store_true')\n"
        "args, extra = p.parse_known_args()\n"
        f"if args.write: Path({str(marker)!r}).write_text('unsafe canonical write')\n")
    receipt = root / "model/result.json"
    receipt.write_text(json.dumps({"sources": {"input": "stale"}, "science": True}))
    original = receipt.read_bytes()
    model = types.ModuleType("model")
    model.verify = types.SimpleNamespace(ROOT=root, REPLAYS=[("pathless", None, "model/result.json")])
    tools = types.ModuleType("tools")
    tools.strict_json = types.SimpleNamespace(loads=json.loads)
    monkeypatch.setitem(sys.modules, "model", model)
    monkeypatch.setitem(sys.modules, "tools", tools)
    monkeypatch.setattr(sys, "argv", ["replay", "--root", str(root), "--output", str(root / "audit"), "--execute"])
    def forbidden(*args, **kwargs):
        pytest.fail("unsupported writer must never be invoked, even with --help")
    monkeypatch.setattr(replay, "subprocess", types.SimpleNamespace(run=forbidden))
    with pytest.raises(SystemExit, match="No producers started"):
        replay.main()
    assert receipt.read_bytes() == original
    assert not marker.exists()
    assert not (root / "audit/before").exists()
    report = json.loads((root / "audit/report.json").read_text())
    assert report["unsupported_affected_writers"] == {"model/result.json": "no reviewed path-output writer contract"}
    assert report["results"] == []


def test_reviewed_writer_parser_changes_require_new_review(tmp_path, monkeypatch):
    source = tmp_path / "writer.py"
    source.write_text("reviewed path-output writer")
    monkeypatch.setattr(replay, "SUPPORTED_WRITERS", {("writer", None): ("writer.py", replay.digest(source))})
    row = ("writer", None, "result.json")
    assert replay.writer_problem(tmp_path, row) is None
    source.write_text("changed parser or side effects")
    assert replay.writer_problem(tmp_path, row) == "reviewed writer source changed: writer.py"
