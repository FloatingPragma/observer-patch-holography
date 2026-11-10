"""Readout maintenance cannot alter source scope or scientific payload."""
import copy,hashlib
import pytest
from tools import renew_source_record_readouts as current


def fixture(tmp_path):
    (tmp_path/'parent.json').write_text('{"current":true}\n');(tmp_path/'source.py').write_text('value=1\n')
    sha=lambda p:hashlib.sha256((tmp_path/p).read_bytes()).hexdigest()
    old={'source_bindings':{'parent.json':'0'*64,'source.py':sha('source.py')},
         'energy':'7/3','physical_species_identified':False,'complete_record_archive':{'sha256':'a'*64}}
    spec={'source_field':'source_bindings','only_renewable_source_paths':['parent.json']}
    return spec,old,{p:sha(p) for p in old['source_bindings']}


def test_one_typed_metadata_parent_preserves_numbers_and_record_archive(tmp_path):
    spec,old,pins=fixture(tmp_path);before=copy.deepcopy(old)
    value=current.project(spec,old,pins,root=tmp_path)
    assert old==before
    assert {k:v for k,v in value.items() if k!='source_bindings'}=={k:v for k,v in old.items() if k!='source_bindings'}


@pytest.mark.parametrize('mutation',['append','delete','code','incorrect_parent'])
def test_unlisted_source_edits_are_rejected(tmp_path,mutation):
    spec,old,pins=fixture(tmp_path)
    if mutation=='append':pins['extra.py']='f'*64
    elif mutation=='delete':pins.pop('source.py')
    elif mutation=='code':pins['source.py']='f'*64
    else:pins['parent.json']='f'*64
    with pytest.raises(ValueError):current.project(spec,old,pins,root=tmp_path)


def test_only_the_four_declared_live_owners_are_admitted():
    value=current.authority()
    assert set(value['owners'])==set(current.OWNERS)
    assert all(not row['only_renewable_source_paths'][0].endswith('.py') for row in value['owners'].values())
