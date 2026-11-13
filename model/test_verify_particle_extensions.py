import json
from types import SimpleNamespace
import pytest
from . import verify_particle_extensions as v

def test_new_shards_partition_without_changing_retained_order():
    for n in (1,2,6,16):
        groups=[v.selected(i,n)for i in range(n)]
        flat=[x for g in groups for x in g]
        assert len(flat)==len(set(flat))==len(v.PACKAGES)
        assert set(flat)==set(v.PACKAGES)
    with pytest.raises(ValueError):v.selected(1,None)
    with pytest.raises(ValueError):v.selected(3,2)

def test_full_inventory_must_match_not_just_recorded_subset(tmp_path,monkeypatch):
    monkeypatch.setattr(v,'ROOT',tmp_path);monkeypatch.setattr(v,'PACKAGES',('test',))
    target=tmp_path/'model/test/results.json';target.parent.mkdir(parents=True)
    monkeypatch.setattr(v.importlib,'import_module',lambda _:SimpleNamespace(pins=lambda:{'a':'one','b':'two'}))
    target.write_text(json.dumps({'source_sha256':{'a':'one'},'verdict':'PASS'}))
    with pytest.raises(ValueError,match='inventory'):v.source_bindings()
    target.write_text(json.dumps({'source_sha256':{'a':'one','b':'two'},'verdict':'FAIL'}))
    with pytest.raises(ValueError,match='producer failed'):v.source_bindings()
    target.write_text(json.dumps({'source_sha256':{'a':'one','b':'two'},'verdict':'PASS'}))
    assert v.source_bindings()==1
