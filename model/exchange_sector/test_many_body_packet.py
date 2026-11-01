"""Public fresh replay, complete receipt mutation and source-custody controls."""
from copy import deepcopy
import json
from types import SimpleNamespace

import pytest

from . import many_body_packet as packet
from . import many_body
from tools import strict_json


@pytest.fixture(scope='module')
def replay():
    # A warmed producer is deliberately wrong. Public replay must execute the
    # source files in its own interpreter, not attach new hashes to this value.
    original = many_body.certificate
    many_body.certificate = lambda: {'fabricated_warm_result': True}
    try:
        result = packet.report()
    finally:
        many_body.certificate = original
    assert 'fabricated_warm_result' not in result['complete_instruments']
    return result


def test_retained_receipt_matches_complete_fresh_execution(replay):
    packet.compare(strict_json.load(packet.OUTPUT), replay)


@pytest.mark.parametrize('key', sorted(packet.RESULT_KEYS))
def test_removing_any_root_receipt_field_fails(replay, key):
    altered = deepcopy(replay)
    del altered[key]
    with pytest.raises(ValueError):
        packet.compare(altered, replay)


@pytest.mark.parametrize('key', ['topology', 'compatibility', 'complete_instruments', 'source_observer'])
def test_trash_proof_or_execution_is_not_admitted(replay, key):
    altered = deepcopy(replay)
    altered[key] = {}
    with pytest.raises(ValueError):
        packet.compare(altered, replay)


def test_probability_resource_and_failure_budget_mutations_fail(replay):
    original = replay['source_observer']['fixed_number_phase_readouts'][0]
    for mutate in (
        lambda r: r.update(complete_output_amplitude_error_ceiling='0'),
        lambda r: r['resources'].update(first_failure_Kraus_slots=0),
        lambda r: r['resources'].update(native_pulses=0),
        lambda r: r['ideal_nonzero_terminal_probabilities'].update({'0': '1/100000000000'}),
    ):
        altered = deepcopy(replay)
        row = altered['source_observer']['fixed_number_phase_readouts'][0]
        mutate(row)
        assert row != original
        with pytest.raises(ValueError):
            packet.compare(altered, replay)


@pytest.mark.parametrize('stdout', ['', '{}', '[]', '{"x":1,"x":2}', '{}\nextra output'])
def test_fresh_worker_protocol_rejects_missing_duplicate_or_extra_output(monkeypatch, stdout):
    monkeypatch.setattr(packet.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout=stdout))
    with pytest.raises((ValueError, RuntimeError)):
        packet.report()


def test_fresh_worker_uses_isolation_and_empty_bytecode_prefix(monkeypatch, replay):
    def owned(command, **kwargs):
        assert command[:5] == [packet.sys.executable, '-I', '-B', '-W', 'error']
        prefix = next(x.split('=', 1)[1] for x in command if x.startswith('pycache_prefix='))
        assert packet.Path(prefix).is_dir()
        assert list(packet.Path(prefix).iterdir()) == []
        assert kwargs['check'] is True
        assert kwargs['env']['OPENBLAS_NUM_THREADS'] == '1'
        return SimpleNamespace(stdout=json.dumps(replay))
    monkeypatch.setattr(packet.subprocess, 'run', owned)
    packet.compare(packet.report(), replay)


def test_source_change_during_child_execution_rejects(monkeypatch, replay):
    before = replay['source_sha256']
    calls = 0
    def bindings():
        nonlocal calls
        calls += 1
        return before if calls == 1 else {**before, 'changed.py': '0'*64}
    monkeypatch.setattr(packet, 'source_bindings', bindings)
    monkeypatch.setattr(packet.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout=json.dumps(replay)))
    with pytest.raises(RuntimeError):
        packet.report()


@pytest.mark.parametrize('key', ['physical_particle_identification', 'integrated_completion_gate_accepted'])
def test_worker_cannot_promote_finite_scope(monkeypatch, replay, key):
    altered = deepcopy(replay)
    altered[key] = True
    monkeypatch.setattr(packet.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout=json.dumps(altered)))
    with pytest.raises(ValueError, match='contract'):
        packet.report()


def test_nonfinite_json_is_not_an_exact_receipt():
    with pytest.raises(ValueError):
        packet.compare({'probability': float('nan')}, {'probability': float('nan')})
