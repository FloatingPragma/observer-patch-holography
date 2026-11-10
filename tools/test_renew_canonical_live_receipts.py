"""The enumerated replay interface retains numerical data and scientific scope."""
from copy import deepcopy
import hashlib
import pytest
from tools import renew_canonical_live_receipts as renewal


def test_unknown_owner_and_modified_catalog_cannot_select_a_producer(monkeypatch):
    with pytest.raises(ValueError, match='unknown canonical live owner'):
        renewal.execute('unreviewed_owner')
    monkeypatch.setattr(renewal, 'CATALOG_SHA256', 'f' * 64)
    with pytest.raises(ValueError, match='authority changed'):
        renewal.load_catalog()


def test_metadata_projection_preserves_science_and_source_scope():
    path = 'tools/strict_json.py'
    old = {'source_sha256': {path: 'old'}, 'diagnostic': 1e-15,
           'scope': {'physical_mass_derived': False}}
    fresh = deepcopy(old)
    fresh['source_sha256'][path] = hashlib.sha256((renewal.ROOT/path).read_bytes()).hexdigest()
    fresh['diagnostic'] = 2e-15
    actual = renewal.carry_metadata('model/native_program/results.json', old, fresh)
    assert actual['diagnostic'] == old['diagnostic']
    assert actual['scope'] == old['scope']
    assert actual['source_sha256'] == fresh['source_sha256']
    fresh['source_sha256']['model/completion/family.py'] = 'extra'
    with pytest.raises(ValueError, match='scope changed'):
        renewal.carry_metadata('model/native_program/results.json', old, fresh)


def test_every_registered_producer_and_retained_archive_has_its_exact_pin():
    for spec in renewal.load_catalog().values():
        assert renewal.digest(renewal.ROOT / spec['producer_source']) == spec['producer_source_sha256']
        assert renewal.digest(renewal.ROOT / spec['retained_path']) == spec['retained_sha256']


def test_transitive_scientific_helper_cannot_change_under_metadata_maintenance(tmp_path):
    helper = tmp_path / 'helper.py'
    helper.write_text('value = 1\n')
    pins = {'helper.py': renewal.digest(helper)}
    renewal.verify_executable_source_pins(pins, tmp_path)
    helper.write_text('value = 2\n')
    with pytest.raises(ValueError, match='transitive executable source changed'):
        renewal.verify_executable_source_pins(pins, tmp_path)
    with pytest.raises(ValueError):
        renewal.verify_executable_source_pins({'../unreviewed.py': 'f' * 64}, tmp_path)


def test_reviewed_four_input_append_is_exact_and_cannot_delete_or_expand():
    receipt = 'model/native_energy_measurement/results.json'
    spec = renewal.load_catalog()['native_energy_measurement']
    old = renewal.strict_json.load(renewal.ROOT/spec['retained_path'])
    additions = spec['reviewed_source_scope_append']
    fresh = deepcopy(old)
    fresh['source_bindings'].update(additions)
    assert renewal.carry_metadata(receipt, old, fresh, additions) == fresh
    for damage in ('delete', 'extra', 'wrong_hash', 'partial_append'):
        bad = deepcopy(fresh)
        if damage == 'delete':
            del bad['source_bindings'][next(iter(old['source_bindings']))]
        elif damage == 'extra':
            bad['source_bindings']['tools/strict_json.py/extra'] = 'f' * 64
        elif damage == 'wrong_hash':
            bad['source_bindings'][next(iter(additions))] = 'f' * 64
        else:
            del bad['source_bindings'][next(iter(additions))]
        with pytest.raises(ValueError):
            renewal.carry_metadata(receipt, old, bad, additions)


def test_only_named_receipt_pointer_can_follow_its_actual_owner():
    spec = renewal.load_catalog()['particle_comparison']
    old = renewal.strict_json.load(renewal.ROOT / spec['retained_path'])
    pointers = spec['reviewed_receipt_pointers']
    fresh = deepcopy(old)
    renewal.set_path(fresh, pointers[0]['json_path'],
                     renewal.digest(renewal.ROOT / pointers[0]['target']))
    assert renewal.carry_metadata(spec['receipt'], old, fresh, pointers=pointers) == fresh
    for damage in ('forged_pointer', 'different_field', 'retained_authority'):
        bad = deepcopy(fresh)
        policy = deepcopy(pointers)
        if damage == 'forged_pointer':
            renewal.set_path(bad, policy[0]['json_path'], 'f' * 64)
        elif damage == 'different_field':
            bad['scope'] = {'unreviewed_acceptance': True}
        else:
            policy[0]['retained_sha256'] = 'f' * 64
        if damage == 'different_field':
            projected = renewal.carry_metadata(spec['receipt'], old, bad, pointers=policy)
            assert renewal.metadata_view(spec['receipt'], projected, policy) != renewal.metadata_view(spec['receipt'], bad, policy)
        else:
            with pytest.raises(ValueError):
                renewal.carry_metadata(spec['receipt'], old, bad, pointers=policy)


def test_array_replay_preserves_json_object_and_original_npz(tmp_path, monkeypatch):
    import json
    from types import SimpleNamespace
    import numpy as np
    (tmp_path / 'producer.py').write_text('pinned producer\n')
    arrays = {'complete': np.array([1.0, 2.0])}
    np.savez(tmp_path / 'results.npz', **arrays)
    original_arrays = (tmp_path / 'results.npz').read_bytes()
    old = {'source_bindings': {'producer.py': renewal.digest(tmp_path / 'producer.py')},
           'scope': {'physical_mass_derived': False},
           'artifact': {'filename': 'results.npz', 'sha256': renewal.digest(tmp_path / 'results.npz')}}
    for name in ('results.json', 'retained.json'):
        (tmp_path / name).write_text(json.dumps(old))
    spec = {'receipt': 'results.json', 'retained_path': 'retained.json',
            'retained_sha256': renewal.digest(tmp_path / 'retained.json'),
            'producer_source': 'producer.py', 'producer_source_sha256': renewal.digest(tmp_path / 'producer.py'),
            'module': 'typed_test_owner', 'entrypoint': 'produce', 'comparison': 'exact'}
    monkeypatch.setattr(renewal, 'ROOT', tmp_path)
    monkeypatch.setattr(renewal, 'load_catalog', lambda: {'test': spec})
    monkeypatch.setattr(renewal.importlib, 'import_module', lambda name: SimpleNamespace(produce=lambda: (deepcopy(old), arrays)))
    assert renewal.execute('test', renew=True)['canonical_scientific_payload_exactly_preserved']
    assert json.loads((tmp_path / 'results.json').read_text()) == old
    assert (tmp_path / 'results.npz').read_bytes() == original_arrays
