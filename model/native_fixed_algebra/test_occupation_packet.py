"""The public catalogue verifier must execute sources and reject changed evidence."""
from copy import deepcopy
import importlib.util
import json
import os
import py_compile
import subprocess
import sys
from types import SimpleNamespace

import pytest

from tools import strict_json
from . import occupation, occupation_packet as packet


@pytest.fixture(scope='module')
def replay():
    original = occupation.certificate
    occupation.certificate = lambda: {'invented_cached_answer': True}
    try:
        result = packet.report()
    finally:
        occupation.certificate = original
    assert 'invented_cached_answer' not in result['algebra_and_protection']
    return result


def test_retained_complete_receipt_matches_owned_execution(replay):
    packet.compare(strict_json.load(packet.OUTPUT), replay)


@pytest.mark.parametrize('key', sorted(packet.RESULT_KEYS))
def test_every_root_field_is_required(replay, key):
    altered = deepcopy(replay)
    del altered[key]
    with pytest.raises(ValueError):
        packet.compare(altered, replay)


@pytest.mark.parametrize('key', ['algebra_and_protection', 'source_observer', 'source_sha256'])
@pytest.mark.parametrize('junk', [{}, [], None, True, 'verified'])
def test_trash_cannot_replace_proof_execution_or_custody(replay, key, junk):
    altered = deepcopy(replay)
    altered[key] = junk
    with pytest.raises(ValueError):
        packet.compare(altered, replay)


@pytest.mark.parametrize('key', sorted(packet.FALSE_KEYS))
def test_no_physical_or_program_admission_promotion(monkeypatch, replay, key):
    altered = deepcopy(replay)
    altered[key] = True
    monkeypatch.setattr(packet.subprocess, 'run', lambda *a, **kw:
                        SimpleNamespace(stdout=json.dumps(altered)))
    with pytest.raises(ValueError, match='contract'):
        packet.report()


@pytest.mark.parametrize('stdout', ['', '{}', '[]', '{"x":1,"x":2}', '{}\ntrailing',
                                    '{"x":NaN}', '{"x":Infinity}'])
def test_malformed_owned_worker_output_is_rejected(monkeypatch, stdout):
    monkeypatch.setattr(packet.subprocess, 'run', lambda *a, **kw: SimpleNamespace(stdout=stdout))
    with pytest.raises((ValueError, RuntimeError)):
        packet.report()


def test_owned_replay_cannot_read_parent_modules_or_existing_bytecode(monkeypatch, replay):
    def worker(command, **kwargs):
        assert command[:5] == [packet.sys.executable, '-I', '-B', '-W', 'error']
        prefix = next(part.split('=', 1)[1] for part in command
                      if part.startswith('pycache_prefix='))
        assert packet.Path(prefix).is_dir() and not list(packet.Path(prefix).iterdir())
        assert kwargs['check'] and kwargs['cwd'] == packet.ROOT
        assert kwargs['env']['OPENBLAS_NUM_THREADS'] == '1'
        return SimpleNamespace(stdout=json.dumps(replay))
    monkeypatch.setattr(packet.subprocess, 'run', worker)
    packet.compare(packet.report(), replay)


def test_owned_worker_flags_bypass_valid_poisoned_timestamp_bytecode(tmp_path, monkeypatch):
    """Exercise real imports while keeping every poisoned byte outside the repo."""
    source = packet.ROOT/'model/native_fixed_algebra/occupation.py'
    stat = source.stat()
    poison_prefix, fresh_prefix = tmp_path/'old-cache', tmp_path/'fresh-cache'
    poison_prefix.mkdir()
    fresh_prefix.mkdir()
    with monkeypatch.context() as context:
        context.setattr(sys, 'pycache_prefix', str(poison_prefix))
        destination = packet.Path(importlib.util.cache_from_source(str(source)))
    assert destination.resolve().is_relative_to(poison_prefix.resolve())
    destination.parent.mkdir(parents=True)
    # Matching the original file size and timestamp makes this an accepted
    # cache candidate, rather than a malformed .pyc Python would reject anyway.
    payload = b"SCHEMA='POISONED_OLD_BYTECODE'\n"
    assert len(payload) < stat.st_size
    poison_source = tmp_path/'poison.py'
    poison_source.write_bytes(payload+b'#'*(stat.st_size-len(payload)))
    os.utime(poison_source, (stat.st_atime, stat.st_mtime))
    py_compile.compile(str(poison_source), cfile=str(destination), dfile=str(source),
                       doraise=True, invalidation_mode=py_compile.PycInvalidationMode.TIMESTAMP)
    body = ('import sys;sys.path.insert(0,'+repr(str(packet.ROOT))+');'
            'from model.native_fixed_algebra import occupation;print(occupation.SCHEMA)')
    environment = dict(os.environ, PYTHONPYCACHEPREFIX=str(poison_prefix),
                       OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')

    def imported_schema(flags):
        return subprocess.run([sys.executable, *flags, '-c', body],
                              cwd=packet.ROOT, env=environment, check=True,
                              capture_output=True, text=True, encoding='utf-8',
                              timeout=60).stdout.strip()

    assert imported_schema(['-B']) == 'POISONED_OLD_BYTECODE'
    assert imported_schema(['-I', '-B', '-W', 'error', '-X', 'utf8', '-X',
                            'pycache_prefix='+str(fresh_prefix)]) == occupation.SCHEMA
    assert occupation.SCHEMA != 'POISONED_OLD_BYTECODE'


def test_source_change_during_execution_is_rejected(monkeypatch, replay):
    before, calls = replay['source_sha256'], 0
    def bindings():
        nonlocal calls
        calls += 1
        return before if calls == 1 else {**before, 'unbound.py': '0'*64}
    monkeypatch.setattr(packet, 'source_bindings', bindings)
    monkeypatch.setattr(packet.subprocess, 'run', lambda *a, **kw:
                        SimpleNamespace(stdout=json.dumps(replay)))
    with pytest.raises(RuntimeError):
        packet.report()


def test_custody_covers_the_fixed_number_contract_and_new_sources(replay):
    for path in ('model/completion/program_class_v2.py',
                 'model/native_fixed_algebra/occupation.py',
                 'model/native_fixed_algebra/occupation_service.py',
                 'model/native_fixed_algebra/occupation_packet.py',
                 'model/native_fixed_algebra/occupation_contract.json'):
        assert path in replay['source_sha256']
    assert 'model/completion/program_class_v2.json' not in replay['source_sha256']
    assert replay['program_class_contract_sha256'] == packet.CLASS_CONTRACT_SHA256
    assert packet.class_contract_digest() == packet.CLASS_CONTRACT_SHA256
    assert all('\\' not in path for path in replay['source_sha256'])


@pytest.mark.parametrize('key', sorted(set(strict_json.load(
    packet.ROOT/'model/completion/program_class_v2.json'))-{'adoption'}))
def test_every_non_adoption_class_clause_remains_pinned(monkeypatch, key):
    value = strict_json.load(packet.ROOT/'model/completion/program_class_v2.json')
    value[key] = {'changed_class_clause': True}
    monkeypatch.setattr(packet.strict_json, 'load', lambda _: value)
    with pytest.raises(ValueError, match='fixed-number V2 class contract'):
        packet.class_contract_digest()


def test_new_class_clause_cannot_evade_the_semantic_pin(monkeypatch):
    value = strict_json.load(packet.ROOT/'model/completion/program_class_v2.json')
    value['additional_class_clause'] = 'arbitrary total number now allowed'
    monkeypatch.setattr(packet.strict_json, 'load', lambda _: value)
    with pytest.raises(ValueError, match='fixed-number V2 class contract'):
        packet.class_contract_digest()


@pytest.mark.parametrize('receipt_sha256', [
    '54101b232466509924b18f3a830ae1524dbb2d6cb0b2e516e92b2035d037fc11',
    '0'*64,
])
def test_changed_interpreter_member_receipt_needs_review(monkeypatch, receipt_sha256):
    value = strict_json.load(packet.ROOT/'model/completion/program_class_v2.json')
    value['interpreter_binding']['source_hashes'][
        'model/native_program/results.json'] = receipt_sha256
    monkeypatch.setattr(packet.strict_json, 'load', lambda _: value)
    with pytest.raises(ValueError, match='fixed-number V2 class contract'):
        packet.class_contract_digest()


def test_independent_adoption_refresh_does_not_make_receipt_hash_recursive(monkeypatch):
    value = strict_json.load(packet.ROOT/'model/completion/program_class_v2.json')
    value['adoption']['independent_review']['sha256'] = 'f'*64
    monkeypatch.setattr(packet.strict_json, 'load', lambda _: value)
    # This is only the class-clause pin: the outer adoption checker must still
    # authenticate any replacement independent review before accepting it.
    assert packet.class_contract_digest() == packet.CLASS_CONTRACT_SHA256


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -float('inf')])
def test_nonfinite_receipts_do_not_pass_by_self_comparison(value):
    with pytest.raises(ValueError):
        packet.compare({'residual': value}, {'residual': value})


def test_exact_json_scalar_types_are_part_of_the_contract():
    for altered in (True, 1.0, '1'):
        with pytest.raises(ValueError):
            packet.compare({'count': altered}, {'count': 1})


@pytest.mark.parametrize('path,junk', [
    (('algebra_and_protection',), None),
    (('source_observer',), []),
    (('algebra_and_protection', 'schema'), 'unowned.schema'),
    (('algebra_and_protection', 'cases'), {}),
    (('algebra_and_protection', 'cases', 'C:+1', 'algebra'), 'verified'),
    (('algebra_and_protection', 'cases', 'E:-1', 'algebra', 'full_witness_replayed'), False),
    (('algebra_and_protection', 'cases', 'C:-1', 'algebra', 'commutant_dimension'), True),
    (('algebra_and_protection', 'cases', 'C:+1', 'algebra', 'witness_sha256'), '0'*63),
    (('algebra_and_protection', 'cases', 'E:+1', 'references'), []),
    (('algebra_and_protection', 'cases', 'E:-1', 'peripheral'), None),
    (('algebra_and_protection', 'scope', 'class_V2_carrier_revised'), True),
    (('algebra_and_protection', 'scope', 'active_recovery_codes_classified'), True),
    (('algebra_and_protection', 'locality', 'moving_particles_with_global_charge_excluded'), True),
    (('algebra_and_protection', 'extremal', 'class_V2_fixed_global_number_admission_satisfied'), True),
    (('source_observer', 'all_actual_complex_columns_and_failure_prefixes_freshly_admitted'), False),
    (('source_observer', 'all_failure_outputs'), 'success only'),
    (('source_observer', 'postselection'), True),
    (('source_observer', 'program_class_v2_execution_claimed'), True),
    (('source_observer', 'physical_particle_or_clock_identified'), True),
    (('source_observer', 'numerical_admissions_are_interval_or_hardware_certificates'), True),
    (('source_observer', 'extreme_cross_number_phase_apparatus_supplied'), False),
    (('source_observer', 'encoded_preparation_preserves_declared_total', 'extreme'), True),
    (('source_observer', 'general_input_repairs'), []),
    (('source_observer', 'experiments'), []),
    (('source_observer', 'experiments', 0, 'work_bits'), True),
    (('source_observer', 'experiments', 0, 'ideal_terminal_laws'), [{}, {}]),
    (('source_observer', 'experiments', 0, 'ideal_terminal_laws', 0), {'0': '1/0'}),
    (('source_observer', 'experiments', 0, 'ideal_terminal_laws', 0), {'0': '1/2'}),
    (('source_observer', 'experiments', 0, 'ideal_terminal_laws', 0), {'00': '1'}),
    (('source_observer', 'experiments', 0, 'ideal_terminal_laws', 0), {'131072': '1'}),
])
def test_malformed_nested_worker_evidence_is_rejected(monkeypatch, replay, path, junk):
    altered = deepcopy(replay)
    target = altered
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = junk
    # The worker-output parser is under test here; custody has separate tests.
    monkeypatch.setattr(packet, 'source_bindings', lambda: replay['source_sha256'])
    monkeypatch.setattr(packet.subprocess, 'run', lambda *a, **kw:
                        SimpleNamespace(stdout=json.dumps(altered)))
    with pytest.raises(ValueError):
        packet.report()


@pytest.mark.parametrize('field', ('general_input_repairs', 'experiments'))
def test_duplicate_worker_cases_cannot_replace_an_unexecuted_case(monkeypatch, replay, field):
    altered = deepcopy(replay)
    altered['source_observer'][field][1] = deepcopy(altered['source_observer'][field][0])
    monkeypatch.setattr(packet, 'source_bindings', lambda: replay['source_sha256'])
    monkeypatch.setattr(packet.subprocess, 'run', lambda *a, **kw:
                        SimpleNamespace(stdout=json.dumps(altered)))
    with pytest.raises(ValueError, match='omitted or duplicated'):
        packet.report()
