"""Retained scientific scope and exact provider replay are binding."""
import copy
from pathlib import Path
import subprocess
import sys

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def result():
    value = packet.build()
    assert packet.canonical(strict_json.load(packet.OUTPUT)) == packet.canonical(value)
    return value


def test_complete_four_provider_replay_preserves_boundaries(result):
    assert result['verification_checks_pass'] is True
    assert tuple(result['providers']) == packet.PROVIDERS
    scaling = result['providers']['spin_particle_scaling']
    assert scaling['provider_status_flags'] == {
        'new_checks_pass': True, 'dependency_custody_pass': False, 'all_pass': False}
    assert result['scope']['all_ancestral_producers_replayed'] is False
    assert result['scope']['physical_particle_species_or_masses_identified'] is False


def test_actual_retained_provider_tampering_rejected(monkeypatch):
    load = packet.strict_json.load
    def tampered(path):
        value = load(path)
        if Path(path).parent.name == 'spin_canonical_field':
            value['theorem']['sharp_mass_shell_proved'] = True
        return value
    monkeypatch.setattr(packet.strict_json, 'load', tampered)
    with pytest.raises(ValueError, match='scope changed'):
        packet.build()


@pytest.mark.parametrize('mutation', ['scope', 'source', 'provider', 'type'])
def test_aggregate_tampering_rejected(result, mutation, monkeypatch):
    value = copy.deepcopy(result)
    if mutation == 'scope':
        value['scope']['all_ancestral_producers_replayed'] = True
    elif mutation == 'source':
        value['source_bindings'][next(iter(value['source_bindings']))] = '0' * 64
    elif mutation == 'provider':
        value['providers']['spin_particle_scaling']['provider_status_flags']['all_pass'] = True
    else:
        value['verification_checks_pass'] = 1
    monkeypatch.setattr(packet, 'build', lambda: result)
    with pytest.raises(ValueError, match='scope changed'):
        packet.validate(value)


def test_source_changes_during_replay_rejected(result, monkeypatch):
    pins = dict(result['source_bindings'])
    changed = dict(pins, injected='0' * 64)
    sequence = iter((pins, changed))
    monkeypatch.setattr(packet, 'source_bindings', lambda: next(sequence))
    with pytest.raises(RuntimeError, match='source changed'):
        packet.build()


def test_explicit_check_path_loaded_before_replay(tmp_path):
    path = tmp_path / 'missing.json'
    run = subprocess.run([sys.executable, '-m', 'model.spin_source_chain', '--check', str(path)],
                         cwd=packet.ROOT, capture_output=True, text=True)
    assert run.returncode != 0
    assert str(path) in run.stderr


def test_duplicate_keys_rejected_before_replay(tmp_path):
    path = tmp_path / 'duplicate.json'
    path.write_text('{"schema": "oph.spin_source_chain.v1", "schema": "forged"}')
    run = subprocess.run([sys.executable, '-m', 'model.spin_source_chain', '--check', str(path)],
                         cwd=packet.ROOT, capture_output=True, text=True)
    assert run.returncode != 0
    assert 'duplicate JSON key' in run.stderr
