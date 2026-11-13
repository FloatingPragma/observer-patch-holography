"""Mutation controls for metadata reuse across complete scientific replays."""

import hashlib
import subprocess
import sys
import types
from pathlib import Path

import pytest

from tools import immutable_custody_replay as replay


def fixture(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    (root / "model").mkdir()
    (root / "tools").mkdir()
    (root / "model/a.py").write_text("x = 1\n")
    (root / ".gitignore").write_text("model/ignored.py\n")
    calls = []

    def original():
        calls.append(True)
        return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted((root / "model").glob("*.py"))}

    module = types.SimpleNamespace(ROOT=root, __name__="fixture", source_bindings=original)
    return root, module, replay.FrozenSnapshot(root), calls


def slot(module, snapshot):
    return replay.PathCache(module, snapshot, frozenset({replay.OUTPUT}))


def test_original_discovery_and_caller_mutation_isolation(tmp_path):
    root, module, snapshot, calls = fixture(tmp_path)
    cached = slot(module, snapshot)
    first = cached()
    expected = module.source_bindings()
    assert first == expected
    first["model/a.py"] = "caller mutation"
    assert cached() == expected and len(calls) == 2
    cached.verify()
    assert cached.final_discovery_verified and len(calls) == 3
    snapshot.check()


def test_every_cache_hit_rehashes_and_mutation_invalidates_acceptance(tmp_path):
    root, module, snapshot, calls = fixture(tmp_path)
    cached = slot(module, snapshot)
    before = cached()
    (root / "model/a.py").write_text("x = 2\n")
    assert cached() != before and len(calls) == 1
    with pytest.raises(ValueError, match="source bytes changed"):
        snapshot.check()


def test_ignored_import_candidate_addition_fails_both_guards(tmp_path):
    root, module, snapshot, _ = fixture(tmp_path)
    cached = slot(module, snapshot)
    cached()
    (root / "model/ignored.py").write_text("x = 3\n")
    with pytest.raises(ValueError, match="filename inventory changed"):
        snapshot.check()
    with pytest.raises(ValueError, match="uncaptured input|path keys changed"):
        cached.verify()


def test_cache_cannot_read_its_mutable_output(tmp_path):
    root, module, snapshot, _ = fixture(tmp_path)
    output = root / replay.OUTPUT
    output.parent.mkdir()
    output.write_text("{}\n")
    module.source_bindings = lambda: {replay.OUTPUT: replay.sha(output)}
    with pytest.raises(ValueError, match="mutable or uncaptured"):
        slot(module, replay.FrozenSnapshot(root))()


def test_changed_root_and_escaping_symlink_fail(tmp_path):
    root, module, snapshot, _ = fixture(tmp_path)
    cached = slot(module, snapshot)
    cached()
    module.ROOT = tmp_path
    with pytest.raises(ValueError, match="root changed"):
        cached()
    module.ROOT = root
    external = tmp_path / "outside.py"
    external.write_bytes((root / "model/a.py").read_bytes())
    (root / "model/a.py").unlink()
    (root / "model/a.py").symlink_to(external)
    with pytest.raises(ValueError, match="escapes"):
        cached()


def test_final_original_discovery_rejects_different_keys(tmp_path):
    root, module, snapshot, _ = fixture(tmp_path)
    cached = slot(module, snapshot)
    cached()
    # The file exists in the immutable snapshot but was not a discovered input.
    cached.original = lambda: {".gitignore": replay.sha(root / ".gitignore")}
    with pytest.raises(ValueError, match="path keys changed"):
        cached.verify()


def test_original_inventory_contract_requires_explicit_review(tmp_path):
    root, _, snapshot, _ = fixture(tmp_path)
    name = next(iter(replay.TARGET_SOURCE_SHA256))
    path = root / (name.replace(".", "/") + ".py")
    path.parent.mkdir()
    path.write_text("def source_bindings(): return {}\n")
    with pytest.raises(ValueError, match="contract requires review"):
        replay.install(snapshot, frozenset(), [])


def test_expected_receipt_is_bound_across_comparison(tmp_path):
    expected = tmp_path / "expected.json"
    expected.write_text('{"value": 1}\n')
    before = replay.sha(expected)
    replay.check_expected(expected, before)
    expected.write_text('{"value": 2}\n')
    with pytest.raises(ValueError, match="expected receipt changed"):
        replay.check_expected(expected, before)


@pytest.mark.parametrize("kind", ["repository", "output", "adapter", "hardlink"])
def test_cli_rejects_audit_aliases_without_overwriting(tmp_path, kind):
    root, _, _, _ = fixture(tmp_path)
    adapter = Path(replay.__file__).resolve()
    before = adapter.read_bytes()
    output = tmp_path / "candidate.json"
    output.write_text("preserve output\n")
    hardlink = tmp_path / "alias.json"
    hardlink.hardlink_to(output)
    audit = {"repository": root / "unsafe.json", "output": output,
             "adapter": adapter, "hardlink": hardlink}[kind]
    result = subprocess.run(
        [sys.executable, str(adapter), "--root", str(root), "--audit", str(audit),
         "--", "--write", str(output)], capture_output=True, text=True,
    )
    assert result.returncode == 2 and "must be distinct" in result.stderr
    assert adapter.read_bytes() == before and output.read_text() == "preserve output\n"
    assert not (root / "unsafe.json").exists()


def test_cli_cannot_overwrite_arbitrary_repository_file(tmp_path):
    root, _, _, _ = fixture(tmp_path)
    target = root / "model/a.py"
    before = target.read_bytes()
    result = subprocess.run(
        [sys.executable, replay.__file__, "--root", str(root), "--audit",
         str(tmp_path / "audit.json"), "--", "--write", str(target)],
        capture_output=True, text=True,
    )
    assert result.returncode == 2 and "declared receipt" in result.stderr
    assert target.read_bytes() == before


def test_owner_origin_must_belong_to_selected_checkout(tmp_path, monkeypatch):
    root, _, _, _ = fixture(tmp_path)
    foreign = tmp_path / "foreign_owner.py"
    foreign.write_text("# Different checkout\n")
    monkeypatch.setattr(
        replay.importlib.util, "find_spec",
        lambda name: types.SimpleNamespace(origin=str(foreign)),
    )
    with pytest.raises(ValueError):
        replay.authenticate_owner(root)


def test_late_guard_failure_preserves_existing_canonical_receipt(tmp_path, monkeypatch):
    import json

    root, _, _, _ = fixture(tmp_path)
    output = root / replay.OUTPUT
    output.parent.mkdir(parents=True)
    original = b'{"accepted": "old receipt"}\n'
    output.write_bytes(original)
    audit = tmp_path / "late-failure-audit.json"

    # Keep real staging, snapshot and audit logic. The scientific runner is a
    # small controlled writer so the test targets a failure after it returns.
    monkeypatch.setattr(replay, "install", lambda *args: None)
    monkeypatch.setattr(replay, "authenticate_owner", lambda selected: selected)
    monkeypatch.setattr(replay, "authenticate_loaded_modules", lambda snapshot: None)
    monkeypatch.setattr(sys, "argv", list(sys.argv))
    monkeypatch.setattr(sys, "path", list(sys.path))
    monkeypatch.chdir(root)

    def write_then_mutate(*args, **kwargs):
        Path(sys.argv[-1]).write_text('{"candidate": "new receipt"}\n')
        (root / "model/a.py").write_text("x = 2\n")

    monkeypatch.setattr(replay.runpy, "run_module", write_then_mutate)
    with pytest.raises(ValueError, match="source bytes changed"):
        replay.main(["--root", str(root), "--audit", str(audit),
                     "--", "--write", str(output)])
    assert output.read_bytes() == original
    assert json.loads(audit.read_text())["verdict"] == "FAIL"


def loaded_module_inventory(monkeypatch, name, module):
    """Isolate only the model/tools entries; preserve Python's other imports."""
    inventory = {key: value for key, value in sys.modules.items()
                 if not (key == "model" or key.startswith("model.")
                         or key == "tools" or key.startswith("tools."))}
    inventory[replay.__name__] = replay
    inventory[name] = module
    monkeypatch.setattr(sys, "modules", inventory)


def test_foreign_preloaded_model_module_is_rejected(tmp_path, monkeypatch):
    root, _, snapshot, _ = fixture(tmp_path)
    foreign = tmp_path / "foreign.py"
    foreign.write_text("x = 1\n")
    module = types.ModuleType("model.foreign_control")
    module.__file__ = str(foreign)
    loaded_module_inventory(monkeypatch, module.__name__, module)
    with pytest.raises(ValueError):
        replay.authenticate_loaded_modules(snapshot)


def test_loaded_module_bytes_must_match_the_frozen_snapshot(tmp_path, monkeypatch):
    root, _, snapshot, _ = fixture(tmp_path)
    module = types.ModuleType("model.a")
    module.__file__ = str(root / "model/a.py")
    loaded_module_inventory(monkeypatch, module.__name__, module)
    replay.authenticate_loaded_modules(snapshot)
    (root / "model/a.py").write_text("x = 2\n")
    with pytest.raises(ValueError):
        replay.authenticate_loaded_modules(snapshot)
