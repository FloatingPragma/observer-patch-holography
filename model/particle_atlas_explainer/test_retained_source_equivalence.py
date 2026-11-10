"""Typed custody equivalence cannot alter an operator, record law or physical ID."""
import copy
import hashlib
import json

import pytest

from . import retained_source_equivalence as admission


def fixture(tmp_path):
    def write(name, value):
        p = tmp_path / name
        p.write_text(json.dumps(value, sort_keys=True) + '\n')
        return hashlib.sha256(p.read_bytes()).hexdigest()
    old_input = write('input-old.json', {'value': 7})
    new_input = write('input.json', {'value': 7, 'source_bindings': {'code.py': 'metadata'}})
    (tmp_path / 'code.py').write_text('H = 7\n')
    code_hash = admission.sha(tmp_path / 'code.py')
    old = {'source_bindings': {'input.json': old_input, 'code.py': code_hash},
           'record': {'energy': '7/3', 'source_id': 'retained-id', 'physical_mass': None}}
    old_hash = write('old.json', old)
    live = copy.deepcopy(old)
    live['source_bindings']['input.json'] = new_input
    write('live.json', live)
    spec = {'retained_path': 'old.json', 'retained_sha256': old_hash,
            'reviewed_nested_receipt_pointers': []}
    return old, live, spec


def test_only_source_custody_differs(tmp_path):
    old, live, spec = fixture(tmp_path)
    result = admission.validate_transition('live.json', spec, root=tmp_path)
    assert result['non_provenance_payload_equal']
    assert old['record'] == live['record']


@pytest.mark.parametrize('mutation', ['energy', 'physical_mass', 'source_id', 'added_source',
                                     'removed_source', 'archive'])
def test_transition_rejects_science_scope_and_archive_changes(tmp_path, mutation):
    old, live, spec = fixture(tmp_path)
    if mutation == 'energy':
        live['record']['energy'] = '8/3'
    elif mutation == 'physical_mass':
        live['record']['physical_mass'] = '0.511 MeV'
    elif mutation == 'source_id':
        live['record']['source_id'] = 'rewritten-id'
    elif mutation == 'added_source':
        (tmp_path / 'unreviewed.py').write_text('H = 8\n')
        live['source_bindings']['unreviewed.py'] = admission.sha(tmp_path / 'unreviewed.py')
    elif mutation == 'removed_source':
        del live['source_bindings']['code.py']
    elif mutation == 'archive':
        (tmp_path / 'old.json').write_text('{}')
    (tmp_path / 'live.json').write_text(json.dumps(live))
    with pytest.raises(ValueError):
        admission.validate_transition('live.json', spec, root=tmp_path)


def test_historical_view_checks_code_and_original_receipt_bytes(tmp_path):
    _, _, spec = fixture(tmp_path)
    approved = {'allowed_live_receipt_transitions': {'live.json': spec},
                'retained_input_sha256': {'live.json': spec['retained_sha256'],
                                         'code.py': admission.sha(tmp_path / 'code.py')}}
    paths = admission.historical_input_paths(root=tmp_path, approved=approved)
    assert paths['live.json'] == tmp_path / 'old.json'
    (tmp_path / 'code.py').write_text('H = 8\n')
    with pytest.raises(ValueError, match='historical source change'):
        admission.historical_input_paths(root=tmp_path, approved=approved)


def member():
    return {'source': {'implementation_sha256': {'code.py': 'old'}, 'H': '7*N'},
            'state': {'source_id': 's', 'density': 'P/Tr(P)'},
            'reference': {'source_id': 's', 'state_id': 'v', 'idle_H': 0},
            'clock': {'source_id': 's', 'formula': 'exp(-itH)'},
            'ids': {'source_id': 's'}, 'member_id': 'm', 'preparation_plan': {'failure': '1/4'}}


def test_only_content_address_fields_are_removed_from_physical_comparison():
    old = member()
    other = copy.deepcopy(old)
    other['source']['implementation_sha256']['code.py'] = 'new'
    other['ids']['source_id'] = 'other'
    other['member_id'] = 'other-member'
    assert admission.physical_specification(old) == admission.physical_specification(other)
    assert old['member_id'] == 'm'


@pytest.mark.parametrize('section,key,value', [
    ('source', 'H', '8*N'), ('state', 'density', 'I/d'),
    ('reference', 'idle_H', 1), ('clock', 'formula', 'exp(-2itH)'),
    ('preparation_plan', 'failure', '0')])
def test_physical_change_survives_projection(section, key, value):
    old, changed = member(), member()
    changed[section][key] = value
    assert admission.physical_specification(old) != admission.physical_specification(changed)
