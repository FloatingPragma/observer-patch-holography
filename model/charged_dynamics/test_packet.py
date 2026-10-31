"""Challenge the complete receipt, not just its headline success flags."""
from copy import deepcopy
import hashlib
import json
import shutil

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def leaves(value, path=()):
    if type(value) is dict:
        for key, item in value.items():
            yield from leaves(item, (*path, key))
    elif type(value) is list:
        for index, item in enumerate(value):
            yield from leaves(item, (*path, index))
    else:
        yield path, value


def test_every_receipt_leaf_is_checked(replay):
    count = 0
    for path, value in leaves(replay):
        mutant = deepcopy(replay)
        parent = mutant
        for key in path[:-1]:
            parent = parent[key]
        parent[path[-1]] = (not value if type(value) is bool else
                            value+1 if type(value) in (int, float) else
                            value+' altered' if type(value) is str else False)
        with pytest.raises(ValueError):
            packet.compare(mutant, replay)
        count += 1
    assert count > 500


def test_receipt_requires_all_fields_exact_types_and_finite_probabilities(replay):
    packet.compare(deepcopy(replay), replay)
    for mutant in ({}, {**replay, 'extra': 1}, {**replay, 'integrated_gate_accepted': 0},
                   {**replay, 'physical_comparison_count_added': False}):
        with pytest.raises(ValueError):
            packet.compare(mutant, replay)
    for value in (-1e-300, float('nan'), float('inf'), -float('inf')):
        actual = {'output_probability': value}
        with pytest.raises(ValueError):
            packet.compare(actual, {'output_probability': 0.})
    mutant = deepcopy(replay)
    bins = mutant['source_history']['rows'][0]['block_number_bins']
    index = min(range(len(bins)), key=bins.__getitem__)
    assert bins[index] < 1e-20
    bins[index] = -1e-300
    with pytest.raises(ValueError, match='invalid charged dynamics probability'):
        packet.compare(mutant, replay)
    with pytest.raises(ValueError):
        packet.compare({'output_probability': 1.1}, {'output_probability': 1.1})


def test_stored_receipt_replays_and_duplicate_keys_fail(replay):
    packet.compare(strict_json.load(packet.OUTPUT), replay)
    with pytest.raises(strict_json.DuplicateKeyError):
        strict_json.loads('{"physical_comparison_count_added":1,"physical_comparison_count_added":0}')
    # Unsupported JSON NaN is rejected by the consuming verifier as well.
    with pytest.raises(ValueError):
        packet.compare(strict_json.loads('{"value":NaN}'), {'value': 0.})


def test_contract_change_cannot_relabel_the_scientific_result(tmp_path, monkeypatch):
    data = packet.contract()
    data['separate_obligations'] = []
    (tmp_path/'contract.json').write_text(json.dumps(data), encoding='utf-8')
    monkeypatch.setattr(packet, 'HERE', tmp_path)
    with pytest.raises(ValueError, match='contract changed'):
        packet.contract()


def test_recursive_source_custody_covers_primitives_and_all_owning_tests():
    bindings = packet.source_bindings()
    required = {p.relative_to(packet.ROOT).as_posix() for p in packet.HERE.glob('*.py')}
    required |= {'model/charged_dynamics/contract.json', 'model/charged_dynamics/README.md',
                 'model/__init__.py', 'model/observer/source_scalar/compiler.py',
                 'model/full_fock_field/field.py', 'model/full_fock_field/local_simulation.py',
                 'model/particle_response/dimer.py', 'model/particle_response/spectrum.py',
                 'tools/strict_json.py', 'model/requirements.txt', 'requirements.txt',
                 'code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json'}
    assert required <= bindings.keys()
    assert len(bindings) == len(set(bindings))
    for relative, digest in bindings.items():
        assert hashlib.sha256((packet.ROOT/relative).read_bytes()).hexdigest() == digest


def test_warm_packet_refuses_new_file_hashes_for_stale_imported_code(tmp_path, monkeypatch):
    original = packet.source_bindings()
    for relative in original:
        destination = tmp_path/relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(packet.ROOT/relative, destination)
    changed = tmp_path/'model/charged_dynamics/chain.py'
    text = changed.read_text(encoding='utf-8')
    old = 'return FieldGeometry(tuple(sides), speed=1/8, interaction=1.)'
    assert old in text
    changed.write_text(text.replace(old, 'return FieldGeometry(tuple(sides), speed=1/7, interaction=1.)'),
                       encoding='utf-8')
    monkeypatch.setattr(packet, 'ROOT', tmp_path)
    monkeypatch.setattr(packet, 'HERE', tmp_path/'model/charged_dynamics')
    assert packet.source_bindings()['model/charged_dynamics/chain.py'] != original['model/charged_dynamics/chain.py']
    # The old imported producer must not be relabelled with the changed file's
    # new hash. Refusal precedes every scientific producer and source service.
    def should_not_execute():
        pytest.fail('stale imported producer was executed')
    monkeypatch.setattr(packet.certification, 'certificate', should_not_execute)
    with pytest.raises(ArithmeticError, match='changed since packet import'):
        packet.report()


def test_workflow_and_replay_registry_cover_the_current_receipt():
    from model import verify
    assert ('model.charged_dynamics.packet', None, 'model/charged_dynamics/results.json') in verify.REPLAYS
    workflow = (packet.ROOT/'.github/workflows/completion-program.yml').read_text()
    assert 'model/charged_dynamics/**' in workflow
    assert 'model/charged_dynamics' in workflow.split('particle-services:', 1)[1]
    assert 'model.charged_dynamics.packet --check model/charged_dynamics/results.json' in workflow
