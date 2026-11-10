"""Hostile conservation tests for the two named nested atlas receipt pointers."""
import copy
import hashlib
import json

import pytest

from tools.renew_atlas_parent_custody import conserve


def fixture(tmp_path):
    (tmp_path / 'old.json').write_text('{"old":true}\n')
    (tmp_path / 'parent.json').write_text('{"current":true}\n')
    (tmp_path / 'source.py').write_text('value = 1\n')
    sha = lambda name: hashlib.sha256((tmp_path / name).read_bytes()).hexdigest()
    retained = {'source_bindings': {'parent.json': sha('old.json'), 'source.py': sha('source.py')},
                'certificate': {'authenticated_source_bindings': {'parent.json': sha('old.json')},
                                'mass': '1/3', 'physical_particle_identified': False}}
    fresh = copy.deepcopy(retained)
    fresh['source_bindings']['parent.json'] = sha('parent.json')
    fresh['certificate']['authenticated_source_bindings']['parent.json'] = sha('parent.json')
    spec = {'receipt': 'output.json', 'reviewed_nested_receipt_pointers': [
        {'json_path': ['certificate', 'authenticated_source_bindings', 'parent.json'],
         'target': 'parent.json', 'retained_path': 'old.json', 'retained_sha256': sha('old.json')}]}
    return spec, retained, fresh


def test_exact_typed_pointer_transition(tmp_path):
    spec, retained, fresh = fixture(tmp_path)
    original = json.dumps(retained)
    assert conserve(spec, retained, fresh, root=tmp_path) == fresh
    assert json.dumps(retained) == original
    assert conserve(spec, retained, retained, allow_retained=True, root=tmp_path) == retained


@pytest.mark.parametrize('mutation', ['mass', 'status', 'extra_nested_pointer', 'extra_source',
                                     'deleted_source', 'wrong_live_pointer', 'altered_archive'])
def test_reject_scientific_or_unreviewed_source_change(tmp_path, mutation):
    spec, retained, fresh = fixture(tmp_path)
    if mutation == 'mass':
        fresh['certificate']['mass'] = '1/2'
    elif mutation == 'status':
        fresh['certificate']['physical_particle_identified'] = True
    elif mutation == 'extra_nested_pointer':
        fresh['certificate']['authenticated_source_bindings']['unreviewed'] = '0' * 64
    elif mutation == 'extra_source':
        (tmp_path / 'new.py').write_text('value = 2\n')
        fresh['source_bindings']['new.py'] = hashlib.sha256((tmp_path / 'new.py').read_bytes()).hexdigest()
    elif mutation == 'deleted_source':
        del fresh['source_bindings']['source.py']
    elif mutation == 'wrong_live_pointer':
        fresh['certificate']['authenticated_source_bindings']['parent.json'] = '0' * 64
    elif mutation == 'altered_archive':
        (tmp_path / 'old.json').write_text('{"old":false}\n')
    with pytest.raises(ValueError):
        conserve(spec, retained, fresh, root=tmp_path)
