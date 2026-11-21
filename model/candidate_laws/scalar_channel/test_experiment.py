"""Adversarial replay checks using isolated files and complete saved receipts."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import re
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
    here = root / 'model/candidate_laws/scalar_channel'
    monkeypatch.setattr(e, 'ROOT', root)
    monkeypatch.setattr(e, 'HERE', here)
    freeze = {'schema': 'oph.native_base_freeze.v1', 'source_sha256': e.source_pins()}
    freeze_path = here / 'FREEZE.json'
    write_json(freeze_path, freeze)
    result = {'schema': 'oph.native_base_replay.v1', 'extended': False,
              'scientific': {'pointer_probabilities': [.75, 0.0, 0.0, .25],
                             'branch_count': 4, 'accepted': True, 'decode_residual': 1e-15}}

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
    'model/candidate_laws/scalar_channel/carrier.py', 'model/candidate_laws/scalar_channel/__init__.py',
    'model/candidate_laws/scalar_channel/protocol.json', 'model/candidate_laws/contrast_memory_spherical_v1/replay.py',
    'model/executable_model/carriers/CELL.npz', 'tools/pin_scope.py', 'code/pin_scope.py',
    'model/executable_model/carrier.py', 'model/candidate_laws/contrast_memory_spherical_v1/law.py', 'model/candidate_laws/contrast_memory_spherical_v1/sectors.py',
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
    assert 'model/candidate_laws/scalar_channel/__init__.py' in pins
    assert not any('/test_' in name or name.endswith('.md') for name in pins)


@pytest.mark.parametrize('mutation', ['missing_source', 'extra_source', 'bad_schema', 'nonmap'])
def test_freeze_integrity_is_required(replay_fixture, mutation):
    f = replay_fixture
    freeze = copy.deepcopy(f['freeze'])
    if mutation == 'missing_source':
        freeze['source_sha256'].pop('model/candidate_laws/scalar_channel/carrier.py')
    elif mutation == 'extra_source':
        freeze['source_sha256']['model/candidate_laws/scalar_channel/foreign.py'] = '0' * 64
    elif mutation == 'bad_schema':
        freeze['schema'] = True
    else:
        freeze['source_sha256'] = False
    write_json(f['freeze_path'], freeze)
    with pytest.raises(ValueError):
        e.main([])


@pytest.mark.parametrize('key,value', [('rest_level', 2), ('rest_level', True),
                                     ('rest_level', 1.0), ('penalty', 1025),
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
                                     'int_readiness', 'wrong_mode', 'added_top_level_field', 'altered_source'])
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
    elif mutation == 'added_top_level_field':
        saved['unreported_claim'] = ['gauge_response']
    else:
        saved['source_sha256']['model/candidate_laws/scalar_channel/carrier.py'] = '0' * 64
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


def test_eigensolver_zero_in_the_contrast_law_is_backend_roundoff():
    # An eigendecomposition gives an exact zero on one LAPACK backend and a
    # value near 1e-34 on another. Neither is a forbidden outcome.
    path = 'result.scientific.contrast_alternative.cell[1].records[0]'
    actual = {'probabilities': [1.0, 0.0]}
    saved = {'probabilities': [1.0, 1.8e-34]}
    e.compare(saved, actual, path)
    e.compare(actual, saved, path)
    with pytest.raises(ValueError):
        e.compare({'probabilities': [1.0 - 1e-6, 1e-6]}, actual, path)


def test_energy_difference_roundoff_is_not_a_forbidden_outcome():
    e.compare({'source_energy_change': 3e-14}, {'source_energy_change': 0.0})
    e.compare({'source_energy_change': 0.0}, {'source_energy_change': 3e-14})
    with pytest.raises(ValueError):
        e.compare({'source_energy_change': 1e-5}, {'source_energy_change': 0.0})


def test_allowed_deterministic_pointer_weights_have_numerical_roundoff():
    actual = {'pointer_probabilities': [0.0, 0.0, 0.0, 1.0]}
    e.compare({'pointer_probabilities': [1e-15, 0.0, 0.0, 1.0-1e-15]}, actual)
    for index in (1, 2):
        invalid = [0.0, 0.0, 0.0, 1.0-1e-15]
        invalid[index] = 1e-15
        with pytest.raises(ValueError):
            e.compare({'pointer_probabilities': invalid}, actual)


def test_saved_residual_cannot_exceed_its_scientific_gate():
    with pytest.raises(ValueError, match='exceeds its gate'):
        e.compare({'decode_residual': 2.1e-11}, {'decode_residual': 1.5e-11})


@pytest.fixture(scope='module', params=['results.json', 'extended.json'])
def complete_receipt(request):
    return e.load_json(Path(__file__).with_name(request.param))


def test_complete_saved_receipt_satisfies_its_scientific_domains(complete_receipt):
    e.compare(complete_receipt, complete_receipt)


RECEIPT_MUTATIONS = [
    ('blank_memory_leakage', 'Scientific field exceeds its gate: .*blank_memory_leakage'),
    ('complex_hopping', 'Scientific field exceeds its gate: .*complex_hopping_formula_dense_error'),
    ('stage_balance', 'Source energy/work balance exceeds its gate'),
    ('stage_continuity', 'Inconsistent saved relation: .*stage_continuity'),
    ('switch_work', 'Inconsistent saved relation: .*external_driving_work'),
    ('total_work', 'Inconsistent saved relation: .*total_external_driving_work'),
    ('maximum_work_residual', 'Inconsistent saved relation: .*largest_source_work_balance_residual'),
    ('raw_fuel_pointer_link', 'Inconsistent saved relation: .*pointer_probabilities'),
    ('raw_fuel_domain', 'Raw fuel probability exceeds its gate'),
    ('clipping_roundoff', 'Inconsistent saved relation: .*largest_probability_roundoff'),
    ('contrast_weight_norm', 'Contrast normalization exceeds its gate'),
    ('contrast_normalization', 'Inconsistent saved relation: .*numerical_normalization'),
    ('contrast_norm_error', 'Inconsistent saved relation: .*state_norm_error'),
    ('contrast_weight_law', 'Inconsistent saved relation: .*probabilities'),
    ('contrast_coarse_probability', 'Inconsistent saved relation: .*no_adjacent_pair_probability'),
]


def mutate_scientific_receipt(receipt, mutation):
    """Break a producer gate or a linked equality within replay tolerances."""
    changed = copy.deepcopy(receipt)
    scientific = changed['scientific']
    cases = scientific['instruments']['cases']
    ledger = cases[0]['records'][0]['source_energy_ledger']
    contrast = scientific['contrast_alternative']['cell'][0]['records'][1]
    if mutation == 'blank_memory_leakage':
        scientific['contrast_alternative']['seam']['blank_memory_leakage'] = 1e-11
    elif mutation == 'complex_hopping':
        scientific['compiler']['compiled_law']['complex_hopping_formula_dense_error'] = 5e-13
    elif mutation == 'stage_balance':
        # The unchanged reported residual cannot certify these raw energies.
        ledger['stages'][0]['source_energy_before'] += 1e-9
    elif mutation == 'stage_continuity':
        # Each stage's own energy difference and work remain consistent.
        ledger['stages'][1]['source_energy_before'] += 1e-9
        ledger['stages'][1]['source_energy_after'] += 1e-9
    elif mutation == 'switch_work':
        ledger['stages'][0]['switch_off_work'] += 5e-12
    elif mutation == 'total_work':
        ledger['total_external_driving_work'] += 5e-12
    elif mutation == 'maximum_work_residual':
        cases[0]['checks']['largest_source_work_balance_residual'] += 5e-12
    elif mutation == 'raw_fuel_pointer_link':
        cases[0]['records'][0]['raw_final_fuel_probability'] += 5e-12
    elif mutation == 'raw_fuel_domain':
        circulation = next(case for case in cases if case['target'] == 'circulation')
        circulation['records'][1]['raw_final_fuel_probability'] = 1 + 3e-11
    elif mutation == 'clipping_roundoff':
        circulation = next(case for case in cases if case['target'] == 'circulation')
        # The clipped pointer law remains valid; its clipping allowance does not.
        circulation['records'][1]['raw_final_fuel_probability'] = 1 + 1e-11
    elif mutation == 'contrast_weight_norm':
        contrast['unnormalized_weights'][1] += 5e-12
    elif mutation == 'contrast_normalization':
        contrast['numerical_normalization'] += 5e-13
    elif mutation == 'contrast_norm_error':
        contrast['state_norm_error'] += 5e-13
    elif mutation == 'contrast_weight_law':
        law = contrast['probabilities']
        first, second = sorted(range(len(law)), key=law.__getitem__, reverse=True)[:2]
        law[first] -= 5e-12
        law[second] += 5e-12
    elif mutation == 'contrast_coarse_probability':
        contrast['no_adjacent_pair_probability'] += 5e-12
    else:
        raise ValueError('Unknown test mutation: ' + mutation)
    return changed


@pytest.mark.parametrize('mutation,message', RECEIPT_MUTATIONS)
def test_saved_scientific_relations_cannot_hide_within_replay_tolerance(
        complete_receipt, mutation, message):
    changed = mutate_scientific_receipt(complete_receipt, mutation)
    # A generated invalid result also fails self-comparison before recording.
    for saved, actual in ((changed, complete_receipt), (complete_receipt, changed),
                          (changed, changed)):
        with pytest.raises(ValueError, match=message):
            e.compare(saved, actual)


@pytest.mark.parametrize('mutation,message', [RECEIPT_MUTATIONS[index] for index in (0, 1, 2, 9, 14)])
def test_optimized_cli_checks_saved_scientific_relations(
        replay_fixture, complete_receipt, mutation, message):
    f = replay_fixture
    actual = copy.deepcopy(complete_receipt)
    actual['source_sha256'] = e.source_pins()
    expected_path = f['root'] / 'expected.json'
    write_json(expected_path, actual)
    write_json(f['receipt'], mutate_scientific_receipt(actual, mutation))
    arguments = ['--verify', str(f['receipt'])]
    if actual['extended']:
        arguments.append('--extended')
    script = (
        'from pathlib import Path\n'
        'from model.candidate_laws.scalar_channel import experiment as e\n'
        f'e.ROOT=Path({str(f["root"])!r})\n'
        f'e.HERE=Path({str(f["here"])!r})\n'
        f'e.scientific_result=lambda extended=False: e.load_json({str(expected_path)!r})\n'
        'try:\n'
        f'    e.main({arguments!r})\n'
        'except ValueError as error:\n'
        '    print(error)\n'
        '    raise SystemExit(7)\n'
        'raise SystemExit(0)\n'
    )
    process = subprocess.run([sys.executable, '-B', '-O', '-c', script], cwd=f['original_root'],
                             text=True, capture_output=True, timeout=30)
    assert process.returncode == 7, process.stdout + process.stderr
    assert re.search(message, process.stdout), process.stdout + process.stderr


def test_optimized_cli_keeps_replay_guards(replay_fixture):
    f = replay_fixture
    saved = copy.deepcopy(f['saved'])
    saved['scientific']['pointer_probabilities'] = [0.75, 1e-16, 0.0, 0.25]
    write_json(f['receipt'], saved)
    script = (
        'import json\nfrom pathlib import Path\n'
        'from model.candidate_laws.scalar_channel import experiment as e\n'
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
