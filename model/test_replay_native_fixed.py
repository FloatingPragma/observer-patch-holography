"""Fixed-output dispatch rejects mismatched receipts and preserves owner checks."""
import ast
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from model import replay_native_fixed as replay, verify


@pytest.mark.parametrize('owner,module', tuple(replay.OWNERS.items()))
def test_actual_owner_output_and_flag_parser(owner, module):
    receipt = f'model/{owner}/results.json'
    command = replay.command_for(owner, receipt)
    assert command == [sys.executable, '-m', module, '--check']
    source = replay.ROOT.joinpath(*module.split('.')).with_suffix('.py')
    tree = ast.parse(source.read_text())
    flags = [node for node in ast.walk(tree)
             if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
             and node.func.attr == 'add_argument'
             and any(isinstance(a, ast.Constant) and a.value == '--check' for a in node.args)]
    assert len(flags) == 1
    assert any(k.arg == 'action' and isinstance(k.value, ast.Constant)
               and k.value.value == 'store_true' for k in flags[0].keywords)


def test_registered_receipts_and_adapter_owners_are_bijective():
    actual = [row for row in verify.REPLAYS if row[0] == 'model.replay_native_fixed']
    expected = {('model.replay_native_fixed', owner, f'model/{owner}/results.json')
                for owner in replay.OWNERS}
    assert len(actual) == len(expected) == 14
    assert set(actual) == expected
    assert not any(row[0] in replay.OWNERS.values() for row in verify.REPLAYS)


@pytest.mark.parametrize('receipt', [
    'model/native_source_chsh/results.json',
    './model/native_local_energy_clock/results.json',
    'model/native_local_energy_clock/../native_local_energy_clock/results.json',
    '/tmp/results.json', '', None,
])
def test_noncanonical_or_wrong_receipt_rejected(receipt):
    with pytest.raises(ValueError):
        replay.command_for('native_local_energy_clock', receipt)


def test_unknown_owner_rejected():
    with pytest.raises(ValueError):
        replay.command_for('os', 'model/os/results.json')


def test_missing_receipt_fails_before_import(tmp_path, monkeypatch):
    monkeypatch.setattr(replay, 'ROOT', tmp_path)
    monkeypatch.setattr(replay.importlib, 'import_module', lambda _: pytest.fail('must not import missing owner'))
    with pytest.raises(FileNotFoundError):
        replay.command_for('native_source_chsh', 'model/native_source_chsh/results.json')


def test_symlink_receipt_rejected(tmp_path, monkeypatch):
    target = tmp_path/'actual.json'
    target.write_text('{}')
    receipt = tmp_path/'model/native_source_chsh/results.json'
    receipt.parent.mkdir(parents=True)
    receipt.symlink_to(target)
    monkeypatch.setattr(replay, 'ROOT', tmp_path)
    with pytest.raises(ValueError, match='symlink'):
        replay.command_for('native_source_chsh', 'model/native_source_chsh/results.json')


def test_changed_producer_output_rejected(monkeypatch):
    monkeypatch.setattr(replay.importlib, 'import_module',
                        lambda _: SimpleNamespace(OUTPUT=Path('/tmp/another.json')))
    with pytest.raises(ValueError, match='OUTPUT'):
        replay.command_for('native_source_chsh', 'model/native_source_chsh/results.json')


def test_actual_argv_and_cwd_preserve_failure(monkeypatch):
    owner = 'native_source_chsh'
    monkeypatch.setattr(sys, 'argv', ['adapter', owner, '--check', f'model/{owner}/results.json'])
    called = []
    def fail(command, **kwargs):
        called.append((command, kwargs))
        raise subprocess.CalledProcessError(7, command)
    monkeypatch.setattr(replay.subprocess, 'run', fail)
    with pytest.raises(subprocess.CalledProcessError) as error:
        replay.main()
    assert error.value.returncode == 7
    assert called == [([sys.executable, '-m', replay.OWNERS[owner], '--check'],
                       {'cwd': replay.ROOT, 'check': True})]
