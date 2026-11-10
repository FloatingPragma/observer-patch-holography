"""Guard the one proof insertion and unchanged exact scientific owner payloads."""
import copy
import json
import pytest
from tools import renew_field_current_documentary_owners as current


def fixture():
    policy = current.authority()
    spec = policy['owners']['spin_dynamical_memory']
    old = json.loads((current.ROOT/spec['retained_path']).read_text())
    fresh = copy.deepcopy(old)
    for _, pins, _ in current.live.source_maps(fresh):
        path = policy['only_documentary_delta']['path']
        if path in pins:
            pins[path] = policy['only_documentary_delta']['current_sha256']
    return policy, spec, old, fresh


def test_exact_documentary_delta_keeps_every_scientific_value():
    policy, spec, old, fresh = fixture()
    original = copy.deepcopy(old)
    assert current.project(spec, old, fresh, policy) == fresh
    assert old == original
    assert len(policy['owners']) == 8
    assert not policy['physical_statuses_changed']


@pytest.mark.parametrize('mutation', ['science', 'added_source', 'lost_source', 'wrong_proof', 'other_source'])
def test_unreviewed_science_or_source_cannot_be_renewed(mutation):
    policy, spec, old, fresh = fixture()
    maps = next(current.live.source_maps(fresh))[1]
    path = policy['only_documentary_delta']['path']
    if mutation == 'science':
        fresh['schema'] += '.changed'
    elif mutation == 'added_source':
        maps['unreviewed.py'] = 'a'*64
    elif mutation == 'lost_source':
        maps.pop(path)
    elif mutation == 'wrong_proof':
        maps[path] = 'a'*64
    else:
        maps[next(name for name in maps if name != path)] = 'b'*64
    with pytest.raises(ValueError):
        current.project(spec, old, fresh, policy)


def test_edited_proof_is_rejected_before_any_producer(tmp_path):
    policy = current.authority()
    for name in (current.AUTHORITY, policy['only_documentary_delta']['retained_path'],
                 policy['only_documentary_delta']['path']):
        dest = tmp_path/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes((current.ROOT/name).read_bytes())
    proof = tmp_path/policy['only_documentary_delta']['path']
    proof.write_bytes(proof.read_bytes()+b'\n')
    with pytest.raises(ValueError, match='reviewed proof changed'):
        current.authority(tmp_path)


def test_only_one_exact_eleven_path_scope_append_is_admitted():
    policy = current.authority()
    appended = {name: row['reviewed_source_scope_append'] for name, row in policy['owners'].items()
                if row.get('reviewed_source_scope_append')}
    assert set(appended) == {'spin_charge_dispersion'}
    assert len(appended['spin_charge_dispersion']) == 11
    assert all(policy['unchanged_source_sha256'][name] == sha
               for name, sha in appended['spin_charge_dispersion'].items())
    numerical = {name for name, row in policy['owners'].items()
                 if row.get('comparison') == 'unchanged_complete_owner_validator'}
    assert numerical == {'source_energy_effort', 'spin_boost_stress'}
