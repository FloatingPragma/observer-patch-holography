"""Adversarial replay checks using isolated files and a synthetic experiment."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from . import experiment as e


def write_json(path, value):
    path.write_text(json.dumps(value, allow_nan=False), encoding='utf-8')


@pytest.fixture
def replay_fixture(tmp_path, monkeypatch):
    original_root = e.ROOT
    root = tmp_path / 'isolated'
    for name in e.source_pins():
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original_root / name, target)
    here = root / 'model/np_scalar'
    monkeypatch.setattr(e, 'ROOT', root)
    monkeypatch.setattr(e, 'HERE', here)
    freeze = {'schema': 'oph.native_base_freeze.v1', 'source_sha256': e.source_pins()}
    freeze_path = here / 'FREEZE.json'
    write_json(freeze_path, freeze)
    result = {'schema': 'oph.native_base_replay.v1', 'extended': False,
              'scientific': {'pointer_probabilities': [.75, 0.0, 0.0, .25],
                             'branch_count': 4, 'accepted': True, 'decode_residual': 1e-15},
              'physical_rows_validated': []}

    def compute(extended=False):
        value = copy.deepcopy(result)
        value['extended'] = extended
        return value

    monkeypatch.setattr(e, 'scientific_result', compute)
    saved = dict(copy.deepcopy(result), source_sha256=e.source_pins())
    receipt = here / 'results.json'
    write_json(receipt, saved)
    return {'root': root, 'here': here, 'freeze': freeze, 'freeze_path': freeze_path,
            'result': result, 'saved': saved, 'receipt': receipt, 'original_root': original_root}


def test_verification_reads_without_rewriting(replay_fixture):
    f = replay_fixture
    before = {p: p.read_bytes() for p in f['root'].rglob('*') if p.is_file()}
    assert e.main([]) == 0
    assert {p: p.read_bytes() for p in f['root'].rglob('*') if p.is_file()} == before


def test_recording_refuses_existing_file_before_scientific_work(replay_fixture, monkeypatch):
    f = replay_fixture
    original = f['receipt'].read_bytes()
    monkeypatch.setattr(e, 'scientific_result', lambda *_: pytest.fail('experiment ran before overwrite guard'))
    with pytest.raises(ValueError, match='Refusing to replace'):
        e.main(['--record', str(f['receipt'])])
    assert f['receipt'].read_bytes() == original


def test_invalid_generated_probability_law_is_not_recorded(replay_fixture, monkeypatch):
    f = replay_fixture
    invalid = copy.deepcopy(f['result'])
    invalid['scientific']['pointer_probabilities'] = [-1e-16, 0.0, 0.0, 1.0]
    monkeypatch.setattr(e, 'scientific_result', lambda *_: invalid)
    output = f['here'] / 'synthetic-invalid.json'
    with pytest.raises(ValueError):
        e.main(['--record', str(output)])
    assert not output.exists()


def test_extended_mode_has_a_distinct_complete_receipt(replay_fixture):
    f = replay_fixture
    saved = copy.deepcopy(f['saved'])
    saved['extended'] = True
    write_json(f['here'] / 'extended.json', saved)
    assert e.main(['--extended']) == 0
    with pytest.raises(ValueError):
        e.main(['--extended', '--verify', str(f['receipt'])])


@pytest.mark.parametrize('name', [
    'model/np_scalar/carrier.py', 'model/np_scalar/__init__.py',
    'model/np_scalar/protocol.json', 'model/np_repair.py',
    'model/np_export/inputs/CELL.npz', 'tools/pin_scope.py', 'code/pin_scope.py',
    'model/polynomial_sv/carrier.py', 'model/polynomial_sv/law.py', 'model/polynomial_sv/sectors.py',
])
def test_production_code_and_input_changes_fail_before_compute(replay_fixture, monkeypatch, name):
    f = replay_fixture
    path = f['root'] / name
    assert name in f['freeze']['source_sha256']
    path.write_bytes(path.read_bytes() + b'\nchanged scientific source\n')
    monkeypatch.setattr(e, 'scientific_result', lambda *_: pytest.fail('unfrozen source executed'))
    with pytest.raises(ValueError, match='source freeze differs'):
        e.main([])


def test_source_inventory_rejects_an_added_production_module(replay_fixture):
    (replay_fixture['here'] / 'unexpected.py').write_text('VALUE = 1\n')
    with pytest.raises(ValueError, match='source freeze differs'):
        e.main([])


def test_scientific_pins_exclude_prose_tests_and_receipts(replay_fixture):
    f = replay_fixture
    pins = e.source_pins()
    for relative in ('README.md', 'PROOF.md', 'test_extra.py', 'FREEZE.json', 'results.json',
                     'extended.json', 'receipts/v1/results.json'):
        path = f['here'] / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('editorial or replay-only bytes\n')
    (f['root'] / 'model/__init__.py').write_text('# unrelated root initializer\n')
    assert e.source_pins() == pins
    assert 'model/np_scalar/__init__.py' in pins
    assert not any('/test_' in name or name.endswith('.md') for name in pins)


@pytest.mark.parametrize('mutation', ['missing_source', 'extra_source', 'bad_schema', 'nonmap'])
def test_freeze_integrity_is_required(replay_fixture, mutation):
    f = replay_fixture
    freeze = copy.deepcopy(f['freeze'])
    if mutation == 'missing_source':
        freeze['source_sha256'].pop('model/np_scalar/carrier.py')
    elif mutation == 'extra_source':
        freeze['source_sha256']['model/np_scalar/foreign.py'] = '0' * 64
    elif mutation == 'bad_schema':
        freeze['schema'] = True
    else:
        freeze['source_sha256'] = False
    write_json(f['freeze_path'], freeze)
    with pytest.raises(ValueError):
        e.main([])


@pytest.mark.parametrize('key,value', [('mass_coefficient', 2), ('mass_coefficient', True),
                                     ('mass_coefficient', 1.0), ('penalty', 1025),
                                     ('loop_sign_seed', 0)])
def test_refreezing_cannot_relabel_hardcoded_scientific_parameters(replay_fixture, monkeypatch, key, value):
    f = replay_fixture
    protocol_path = f['here'] / 'protocol.json'
    protocol = e.load_json(protocol_path)
    protocol['scalar'][key] = value
    write_json(protocol_path, protocol)
    freeze = copy.deepcopy(f['freeze'])
    freeze['source_sha256'] = e.source_pins()
    write_json(f['freeze_path'], freeze)
    monkeypatch.setattr(e, 'scientific_result', lambda *_: pytest.fail('mislabeled scientific protocol executed'))
    with pytest.raises(ValueError):
        e.main([])


@pytest.mark.parametrize('mutation', ['missing_branch', 'missing_field', 'unknown_field', 'bool_count',
                                     'int_readiness', 'wrong_mode', 'awarded_row', 'altered_source'])
def test_receipt_structure_types_and_source_identity_are_exact(replay_fixture, mutation):
    f = replay_fixture
    saved = copy.deepcopy(f['saved'])
    if mutation == 'missing_branch':
        saved['scientific']['pointer_probabilities'].pop()
    elif mutation == 'missing_field':
        saved['scientific'].pop('accepted')
    elif mutation == 'unknown_field':
        saved['scientific']['unreported_intervention'] = 0
    elif mutation == 'bool_count':
        saved['scientific']['branch_count'] = True
    elif mutation == 'int_readiness':
        saved['scientific']['accepted'] = 1
    elif mutation == 'wrong_mode':
        saved['extended'] = True
    elif mutation == 'awarded_row':
        saved['physical_rows_validated'] = ['gauge_response']
    else:
        saved['source_sha256']['model/np_scalar/carrier.py'] = '0' * 64
    write_json(f['receipt'], saved)
    with pytest.raises(ValueError):
        e.main([])


@pytest.mark.parametrize('text', ['{"key":0,"key":1}', '{"nested":{"key":0,"key":1}}',
                                '{"value":NaN}', '{"value":Infinity}', '{"value":-Infinity}',
                                '{"value":1e999}'])
def test_json_loader_rejects_ambiguous_or_nonfinite_input(tmp_path, text):
    path = tmp_path / 'bad.json'
    path.write_text(text)
    with pytest.raises(ValueError):
        e.load_json(path)


@pytest.mark.parametrize('key', ['pointer_probabilities', 'probabilities', 'complete_basis_law'])
@pytest.mark.parametrize('bad', [[True, 0.0], [float('nan'), 0.0], [float('inf'), 0.0],
                               [-1e-16, 1.0], [0.0, 1.0 + 1e-15], [], [0.0, 0.0]])
def test_probability_vectors_reject_invalid_domains_on_both_sides(key, bad):
    actual = {key: [1e-16, 1.0 - 1e-16]}
    with pytest.raises(ValueError):
        e.compare({key: bad}, actual)
    with pytest.raises(ValueError):
        e.compare(actual, {key: bad})


def test_probability_mass_drift_cannot_hide_in_element_tolerances():
    actual = {'pointer_probabilities': [1 / 64] * 64}
    saved = {'pointer_probabilities': [1 / 64 + 1e-12] * 64}
    with pytest.raises(ValueError):
        e.compare(saved, actual)


def test_forbidden_pointer_support_stays_exact():
    actual = {'pointer_probabilities': [1.0, 0.0]}
    saved = {'pointer_probabilities': [1.0 - 1e-16, 1e-16]}
    with pytest.raises(ValueError):
        e.compare(saved, actual)


def test_energy_difference_roundoff_is_not_a_forbidden_outcome():
    e.compare({'source_energy_change': 3e-14}, {'source_energy_change': 0.0})
    e.compare({'source_energy_change': 0.0}, {'source_energy_change': 3e-14})
    with pytest.raises(ValueError):
        e.compare({'source_energy_change': 1e-5}, {'source_energy_change': 0.0})


def test_saved_residual_cannot_exceed_its_scientific_gate():
    with pytest.raises(ValueError, match='exceeds its gate'):
        e.compare({'decode_residual': 2e-12}, {'decode_residual': 1e-15})


def test_optimized_cli_keeps_replay_guards(replay_fixture):
    f = replay_fixture
    saved = copy.deepcopy(f['saved'])
    saved['scientific']['pointer_probabilities'] = [0.75, 1e-16, 0.0, 0.25]
    write_json(f['receipt'], saved)
    script = (
        'import json\nfrom pathlib import Path\n'
        'from model.np_scalar import experiment as e\n'
        f'e.ROOT=Path({str(f["root"])!r})\n'
        f'e.HERE=Path({str(f["here"])!r})\n'
        f'e.scientific_result=lambda extended=False: json.loads({json.dumps(f["result"])!r})\n'
        'try:\n    e.main([])\n'
        'except ValueError:\n    raise SystemExit(7)\n'
        'raise SystemExit(0)\n'
    )
    process = subprocess.run([sys.executable, '-B', '-O', '-c', script], cwd=f['original_root'],
                             text=True, capture_output=True, timeout=30)
    assert process.returncode == 7, process.stdout + process.stderr
