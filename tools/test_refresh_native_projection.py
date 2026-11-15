from copy import deepcopy
import pytest
from tools.refresh_native_projection import (
 classify_changes,retained_collector,completion_projection,NAMES,PAPER)

@pytest.mark.parametrize('name',[
 'model/openwave_native_coherent_clock/source.py',
 'model/openwave_native_coherent_clock/results.json',
 'Lean/OPH/Invariant.lean','claims/selection_ledger.json',
 'tracking/premise_register.json','claims/physical_identification_registry.json',
 'code/particles/scripts/build_postdiction_ledger.py'])
def test_existing_operative_drift_cannot_be_blessed_by_review_manifest(name):
    with pytest.raises(ValueError,match='operative'):classify_changes({name:'a'},{name:'b'})

def test_deleted_and_unknown_inputs_are_rejected():
    with pytest.raises(ValueError,match='deleted'):classify_changes({'old':'a'}, {})
    with pytest.raises(ValueError,match='added'):classify_changes({}, {'Lean/Extra.lean':'a'})
    with pytest.raises(ValueError,match='added'):classify_changes({}, {'model/native_gauss_walk/override.py':'a'})
    assert classify_changes({'README.md':'a'},{'README.md':'b'})=={'README.md':'b'}

def test_retained_collector_exact_values_not_just_labels():
    old={'receipts_sha256':{'old':'a'},'science':{'value':1,'failed':False}}
    fresh=deepcopy(old);pins={f'model/{n}/results.json':n for n in NAMES}
    fresh['receipts_sha256'].update(pins)
    fresh['joint_observer_certificates']={'source_packages':{n:{}for n in NAMES}}
    assert retained_collector(old,fresh,pins)
    for path,value in [(('science','value'),2),(('science','failed'),True),(('receipts_sha256','old'),'b')]:
        mutated=deepcopy(fresh);mutated[path[0]][path[1]]=value
        with pytest.raises(ValueError,match='retained'):retained_collector(old,mutated,pins)
    bad=deepcopy(fresh);bad['receipts_sha256'][next(iter(pins))]='wrong'
    with pytest.raises(ValueError,match='binding'):retained_collector(old,bad,pins)

def test_only_authenticated_same_paper_hash_can_renew():
    old={'decision':'blocked','retained_foundations':[{'source_sha256':{}}for _ in range(5)]}
    old['retained_foundations'][4]['source_sha256'][PAPER]='old'
    fresh=deepcopy(old);fresh['retained_foundations'][4]['source_sha256'][PAPER]='current'
    report=completion_projection(old,fresh,'old','current')
    assert report['path']==PAPER
    for mutate in ('decision','stale','extra'):
        bad=deepcopy(fresh)
        if mutate=='decision':bad['decision']='complete'
        elif mutate=='stale':bad['retained_foundations'][4]['source_sha256'][PAPER]='stale'
        else:bad['retained_foundations'][3]['source_sha256'][PAPER]='current'
        with pytest.raises(ValueError):completion_projection(old,bad,'old','current')

def test_file_type_and_executable_modes_are_bound(tmp_path):
    from tools.refresh_native_projection import read_path
    p=tmp_path/'source';p.write_text('target');p.chmod(0o644)
    assert read_path(tmp_path,'source','100644')==b'target'
    with pytest.raises(ValueError,match='executable'):read_path(tmp_path,'source','100755')
    with pytest.raises(ValueError,match='symlink'):read_path(tmp_path,'source','120000')
    p.unlink();p.symlink_to('target')
    assert read_path(tmp_path,'source','120000')==b'target'
    with pytest.raises(ValueError,match='symlink'):read_path(tmp_path,'source','100644')
    with pytest.raises(ValueError,match='symlink'):read_path(tmp_path,'source')
    p.unlink();p.mkdir()
    with pytest.raises(ValueError,match='missing'):read_path(tmp_path,'source','100644')

def test_foreign_root_cannot_reuse_module_global_owners(tmp_path):
    from tools.refresh_native_projection import refresh
    with pytest.raises(ValueError,match='producer module root'):
        refresh(tmp_path,'irrelevant',{},write=False)

def test_pure_parser_cache_keeps_actual_closure_and_invalidates_changed_text(tmp_path,monkeypatch):
    import ast
    from model.particle_common_contract import provenance
    from tools.refresh_native_projection import source_ast_cache
    root=tmp_path/'source';(root/'model').mkdir(parents=True)
    (root/'model/__init__.py').write_text('')
    (root/'model/a.py').write_text('from . import b\n')
    (root/'model/b.py').write_text('x=1\n')
    (root/'model/c.py').write_text('y=2\n')
    monkeypatch.setattr(provenance,'ROOT',root)
    original=provenance.ast
    baseline=provenance.import_closure({'model/a.py'})
    with source_ast_cache()as stats:
        assert ast.parse is original.parse
        assert provenance.import_closure({'model/a.py'})==baseline
        assert provenance.import_closure({'model/a.py'})==baseline
        (root/'model/b.py').write_text('from . import c\n')
        changed=provenance.import_closure({'model/a.py'})
        assert 'model/c.py'in changed and changed>baseline
        # Full parser arguments, including mode, form the cache key.
        assert isinstance(provenance.ast.parse('x','f','eval'),ast.Expression)
        assert isinstance(provenance.ast.parse('x','f','exec'),ast.Module)
    assert provenance.ast is original and ast.parse is original.parse
    assert stats['hits']>0 and stats['misses']>=5
    assert provenance.import_closure({'model/a.py'})==changed
    with pytest.raises(RuntimeError):
        with source_ast_cache():
            raise RuntimeError('controlled failure')
    assert provenance.ast is original
