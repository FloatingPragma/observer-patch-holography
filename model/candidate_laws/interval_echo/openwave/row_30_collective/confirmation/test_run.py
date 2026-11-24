"""Check launcher custody and failure statuses using disposable official-history fixtures."""

import json
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pytest

from model.candidate_laws.interval_echo.openwave.row_30_collective.confirmation import (
    run,
)


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args]).decode().strip()


def commit(root):
    git(root, "add", "-A")
    git(
        root,
        "-c",
        "user.name=fixture",
        "-c",
        "user.email=fixture@example.invalid",
        "-c",
        "commit.gpgsign=false",
        "commit",
        "-qm",
        "fixture",
    )
    return git(root, "rev-parse", "HEAD")


def write(root, rel, value):
    file = root / rel
    file.parent.mkdir(parents=True, exist_ok=True)
    file.write_text(value)
    return file


@pytest.fixture
def official(tmp_path, monkeypatch):
    p = tmp_path / "platform"
    p.mkdir()
    git(p, "init", "-qb", "main")
    write(p, "README", "baseline")
    monkeypatch.setattr(run, "BASE", commit(p))
    monkeypatch.setattr(
        run, "official_head", lambda root: git(root, "rev-parse", "main")
    )
    return p


def test_never_landed_port_exits_before_vendor_import(official, monkeypatch):
    p = official
    port = p / run.PORT
    port.mkdir(parents=True)
    monkeypatch.setattr(run, "HERE", port)

    def surprise(*args):
        pytest.fail("unaccepted port imported code or derived selectors")

    monkeypatch.setattr(run.importlib, "import_module", surprise)
    monkeypatch.setattr(run, "public_settings", surprise)
    args = SimpleNamespace(
        command="acquire", mode="public", openwave=None, source=p, release=p / "unused"
    )
    with pytest.raises(run.NoPublicTask):
        run.run(args)
    assert not args.release.exists()


def test_official_port_tracks_launcher_and_every_vendor_blob(official):
    p = official
    write(p, run.TASK, "candidate terms")
    commit(p)
    launcher = write(p, run.PORT / "run.py", "fixture launcher")
    vendor = write(p, run.PORT / "vendor/module.py", "fixture numerical source")
    with pytest.raises(ValueError, match="official|Official"):
        run.authenticate_port(p)
    approved = commit(p)
    assert run.authenticate_port(p) == approved
    vendor.write_text("working edit")
    with pytest.raises(ValueError, match="official blob"):
        run.authenticate_port(p)
    git(p, "restore", str(run.PORT))
    launcher.write_text("working launcher edit")
    with pytest.raises(ValueError, match="official blob"):
        run.authenticate_port(p)
    git(p, "restore", str(run.PORT))
    write(p, run.PORT / "vendor/untracked.py", "undeclared")
    with pytest.raises(ValueError, match="inventory"):
        run.authenticate_port(p)


def test_withdrawn_task_history_cannot_return_pending_for_absent_port(official):
    p = official
    write(p, run.TASK, "landed")
    commit(p)
    (p / run.TASK).unlink()
    commit(p)
    with pytest.raises(ValueError, match="Official OpenWave launcher"):
        run.authenticate_port(p)


def test_official_symlink_is_never_a_runtime_blob(official):
    p = official
    write(p, run.TASK, "landed")
    write(p, run.PORT / "run.py", "fixture")
    target = write(p, run.PORT / "vendor/a.py", "data")
    (target.parent / "b.py").symlink_to("a.py")
    commit(p)
    with pytest.raises(ValueError, match="nonregular"):
        run.authenticate_port(p)


def test_unlisted_import_is_blocked_before_module_body(tmp_path):
    body = "raise AssertionError('unfrozen module was executed')"
    write(tmp_path, "model/phantom_probe.py", body)
    guard = run.FrozenModelImports(tmp_path, set())
    with pytest.raises(ValueError, match="Unbound model import"):
        guard.find_spec("model.phantom_probe", [str(tmp_path / "model")])


def test_missing_science_is_a_distinct_incomplete_failure(monkeypatch):
    def absent(name):
        raise ModuleNotFoundError(f"No module named {name}")

    monkeypatch.setattr(run.importlib, "import_module", absent)
    with pytest.raises(run.ScientificUnavailable, match="numpy"):
        run.scientific_engine()


def test_malformed_source_manifest_cli_returns_invalid(tmp_path):
    git(tmp_path, "init", "-qb", "main")
    place = tmp_path / run.LOCAL
    place.mkdir(parents=True)
    shutil.copyfile(run.__file__, place / "run.py")
    (place / "files.json").write_text("[]\n")
    (place / "TASK.md").write_text("fixture")
    commit(tmp_path)
    result = subprocess.run(
        [sys.executable, "-I", "-B", str(place / "run.py"), "check-sources"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 4 and json.loads(result.stderr)["status"] == "invalid"
