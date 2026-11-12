"""The documentary adapter cannot admit changed science or executable inputs."""
from copy import deepcopy
from contextlib import contextmanager
from pathlib import Path
import json
from types import SimpleNamespace

import pytest

from model import native_control_custody as custody
from model import verify, verify_current


@pytest.fixture
def source_view(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[1]
    historical = (root / custody.CONTROL).read_bytes()
    control = json.loads(historical)
    live = json.loads((root / custody.LIVE).read_text())
    documents = custody.views.authenticated_documents(root, documentary_paths=(custody.DOCUMENT,))
    for name in (*control["source_sha256"], custody.CONTROL, custody.LIVE):
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((root / name).read_bytes())
    monkeypatch.setattr(custody.subprocess, "check_output", lambda *args, **kwargs: historical)
    monkeypatch.setattr(custody.views, "authenticated_documents",
                        lambda *args, **kwargs: deepcopy(documents))
    return tmp_path, live


def test_authenticates_all_six_current_executables_and_full_science(source_view):
    root, live = source_view
    assert set(custody.authenticate(root, live)) == {custody.DOCUMENT}


@pytest.mark.parametrize("damage", ["control", "science", "type", "inventory", "executable", "live_pin"])
def test_every_non_documentary_change_fails(source_view, damage):
    root, live = source_view
    if damage == "control":
        (root / custody.CONTROL).write_bytes((root / custody.CONTROL).read_bytes() + b" ")
    elif damage == "executable":
        name = next(p for p in live["source_sha256"] if p != custody.DOCUMENT)
        (root / name).write_bytes((root / name).read_bytes() + b"\n")
    else:
        if damage == "science":
            live["scope"]["physical_energy_or_mass_derived"] = True
        elif damage == "type":
            live["carrier"]["states"] = float(live["carrier"]["states"])
        elif damage == "inventory":
            live["source_sha256"].pop(next(p for p in live["source_sha256"] if p != custody.DOCUMENT))
        elif damage == "live_pin":
            live["source_sha256"][custody.DOCUMENT] = "0" * 64
        (root / custody.LIVE).write_text(json.dumps(live))
    with pytest.raises(ValueError):
        custody.authenticate(root, live)


def test_only_exact_native_document_error_enters_adapter(monkeypatch):
    calls = []
    def original():
        calls.append("original")
        raise ValueError(custody.SOURCE_ERROR)
    monkeypatch.setattr(verify, "verify_artifacts", original)
    monkeypatch.setattr(custody, "verify_in_documentary_view", lambda retained: calls.append("adapter") or 17)
    monkeypatch.setattr(verify_current, "verify_successor_source_bindings", lambda: calls.append("successor"))
    assert verify_current.verify_all_hashes() == 17
    assert calls == ["original", "adapter", "successor"]


def test_live_replay_must_match_actual_complete_producer(source_view):
    root, live = source_view
    altered = deepcopy(live)
    altered["carrier"]["states"] += 1
    with pytest.raises(ValueError, match="live native producer replay"):
        custody.authenticate(root, altered)


@pytest.mark.parametrize("real_root,sharded", [(True, True), (True, False), (False, True)])
def test_documentary_plugin_is_scoped_to_real_repository_shards(tmp_path, real_root, sharded):
    from model import pytest_shards
    calls = []
    options = {"model_shard_index": 0 if sharded else None,
               "model_shard_count": 16 if sharded else None}
    config = SimpleNamespace(rootpath=Path(__file__).resolve().parents[1] if real_root else tmp_path,
                             getoption=lambda name, default=None: options.get(name, default),
                             pluginmanager=SimpleNamespace(import_plugin=calls.append))
    pytest_shards.pytest_configure(config)
    assert calls == (["tools.completion_contract_inputs"] if real_root and sharded else [])


def test_nightly_dispatch_preserves_cold_modules_arguments_and_atomic_worker(tmp_path, monkeypatch):
    from tools import completion_contract_inputs
    calls = []
    ordinary = ("model.example", "subcommand", "model/example/results.json")
    groups = (verify.EMPIRICAL_REPLAYS, (custody.REPLAY,), (ordinary,))
    monkeypatch.setattr(verify_current, "selected_groups", lambda *args: groups)
    monkeypatch.setattr(verify_current, "verify_all_hashes", lambda: 17)
    monkeypatch.setattr(verify_current.sys, "argv", ["verify_current"])
    monkeypatch.setattr(verify_current.subprocess, "run", lambda argv, **kwargs: calls.append(argv))
    monkeypatch.setattr(custody, "replay_live_control", lambda root: calls.append("native_live_replay"))
    @contextmanager
    def inputs(root):
        calls.append("view_enter")
        yield tmp_path
        calls.append("view_verified")
    monkeypatch.setattr(completion_contract_inputs, "finite_inputs", inputs)
    verify_current.main()
    assert calls == [
        "view_enter",
        [verify_current.sys.executable, "-m", "model.replay_empirical"],
        "native_live_replay",
        [verify_current.sys.executable, "-m", "tools.finite_contract_command", "--view", str(tmp_path),
         "--module", ordinary[0], "--", ordinary[1], "--check", ordinary[2]],
        "view_verified",
    ]
