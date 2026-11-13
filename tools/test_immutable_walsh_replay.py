"""Fail-closed arithmetic-adapter routing; full original science runs separately."""
from pathlib import Path
import hashlib
import json
import os
import sys
from types import SimpleNamespace

import pytest

from tools import immutable_walsh_replay as r
from tools import immutable_custody_pytest as harness


def test_every_reviewed_helper_is_authenticated():
    modules,paths=r.load_siblings()
    assert set(modules)==set(r.SIBLING_SHA256)
    assert all(r.sha(paths[name])==digest for name,digest in r.SIBLING_SHA256.items())


def test_changed_helper_rejected_before_import(monkeypatch):
    monkeypatch.setattr(r,'sha',lambda _: '0'*64)
    with pytest.raises(ValueError,match='renewed review'):r.load_siblings()


def test_final_helper_mutation_rejected(tmp_path):
    path=tmp_path/'exact_walsh.py';path.write_text('changed')
    with pytest.raises(ValueError,match='changed during execution'):
        r.verify_helpers({'exact_walsh.py':path})


@pytest.mark.parametrize('name',['PYTEST_ADDOPTS','PYTEST_PLUGINS'])
def test_environment_cannot_change_selected_tests(tmp_path,monkeypatch,name):
    monkeypatch.setenv(name,'-k omitted')
    with pytest.raises(SystemExit):
        r.main(['--root',str(tmp_path),'--audit',str(tmp_path.parent/'audit.json')])


def test_audit_cannot_overwrite_repository_source(tmp_path):
    with pytest.raises(SystemExit):
        r.main(['--root',str(tmp_path),'--audit',str(tmp_path/'owned.json')])


@pytest.mark.parametrize('failure',['none','math','metadata','expected','snapshot','helper'])
def test_final_guards_restore_runtime_and_label_algorithm_change(tmp_path,monkeypatch,failure):
    root=tmp_path/'repo';root.mkdir();(root/'expected.json').write_text('{}')
    output=tmp_path/'audit.json';state={'installed':False,'restored':False,'math_verified':False}
    digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    class Snapshot:
        names=('expected.json',)
        bindings={'expected.json':digest(root/'expected.json')}
        def __init__(self,root):self.root=root
        def check(self):
            if failure=='snapshot':raise ValueError('snapshot failure')
    class Math:
        def __init__(self,*args):pass
        def install(self):state['installed']=True
        def verify(self):
            if failure=='math':raise ValueError('math failure')
            state['math_verified']=True
        def restore(self):state['restored']=True;state['installed']=False
        def evidence(self):return {'algorithm_changed':True}
    class Metadata:
        module=SimpleNamespace(source_bindings='adapted')
        original='original'
        def verify(self):
            if failure=='metadata':raise ValueError('metadata failure')
        def evidence(self):return {'fresh_hashes':True}
    slot=Metadata()
    def check_expected(path,expected):
        if failure=='expected':raise ValueError('expected failure')
        assert digest(path)==expected
    guard=SimpleNamespace(OUTPUT='expected.json',FrozenSnapshot=Snapshot,
        same_file=lambda a,b:Path(a)==Path(b),authenticate_owner=lambda _:None,
        authenticate_loaded_modules=lambda _:None,
        install=lambda snapshot,forbidden,slots:slots.append(slot),check_expected=check_expected)
    reviewed=SimpleNamespace(TEST=harness.TEST,TEST_SHA256=harness.TEST_SHA256,
        PYTEST_INI_SHA256=harness.PYTEST_INI_SHA256,ABSENT_CONFTESTS=harness.ABSENT_CONFTESTS,
        ExactTests=harness.ExactTests,authenticate_test=lambda _:None,authenticate_harness=lambda _:None)
    modules={'immutable_custody_replay.py':guard,'immutable_custody_pytest.py':reviewed,
             'exact_walsh.py':SimpleNamespace(SynthesisAdapter=Math)}
    monkeypatch.setattr(r,'load_siblings',lambda:(modules,{}))
    def helpers(_):
        if failure=='helper':raise ValueError('helper failure')
    monkeypatch.setattr(r,'verify_helpers',helpers)
    import_original=r.importlib.import_module
    monkeypatch.setattr(r.importlib,'import_module',lambda name,*a,**kw:SimpleNamespace() if name=='model.native_effect_synthesis.synthesis' else import_original(name,*a,**kw))
    def pretend_pytest(args,plugins):
        assert args[0]==harness.TEST and '--noconftest' in args
        plugin=plugins[0]
        plugin.pytest_collection_finish(SimpleNamespace(items=[SimpleNamespace(nodeid=node,path=root/harness.TEST) for node in sorted(harness.EXPECTED_IDS)]))
        for node in harness.EXPECTED_IDS:
            for phase in ('setup','call','teardown'):
                plugin.pytest_runtest_logreport(SimpleNamespace(nodeid=node,when=phase,outcome='passed'))
        return 0
    monkeypatch.setattr(pytest,'main',pretend_pytest)
    monkeypatch.delenv('PYTEST_ADDOPTS',raising=False);monkeypatch.delenv('PYTEST_PLUGINS',raising=False)
    monkeypatch.setenv('PYTEST_DISABLE_PLUGIN_AUTOLOAD','previous')
    before=(Path.cwd(),list(sys.path),list(sys.argv),os.environ['PYTEST_DISABLE_PLUGIN_AUTOLOAD'])
    arguments=['--root',str(root),'--audit',str(output)]
    if failure=='none':r.main(arguments)
    else:
        with pytest.raises(ValueError,match=failure+' failure'):r.main(arguments)
    assert state['restored'] and not state['installed']
    assert slot.module.source_bindings=='original'
    assert (Path.cwd(),list(sys.path),list(sys.argv),os.environ['PYTEST_DISABLE_PLUGIN_AUTOLOAD'])==before
    audit=json.loads(output.read_text())
    assert audit['arithmetic_algorithm_changed'] and not audit['metadata_only_replay']
    assert not audit['native_observer_clock_or_resource_counts_changed']
    assert not audit['scientific_result_cache'] and not audit['comparison_function_changed']
    assert audit['verdict'].startswith('PASS_') if failure=='none' else audit['verdict']=='FAIL'
