"""Receipt custody and typed mutations."""
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import pytest
from . import packet


@pytest.fixture(scope='module')
def receipt():return packet.produce()


def test_retained_receipt_replays(receipt):
    assert packet.canonical(receipt)==packet.canonical(packet.strict_json.load(packet.OUTPUT))


@pytest.mark.parametrize('change',('minorization','archive','type','source'))
def test_receipt_mutations_fail(receipt,tmp_path,change):
    bad=deepcopy(receipt)
    if change=='minorization':bad['exact_all_sector_minorization'][6]['minorization_numerator']+=1
    elif change=='archive':bad['scope']['mixing_archive_traced_only_for_comparison']=False
    elif change=='type':bad['chosen_budget']['attempts']=1931.0
    else:bad['source_bindings'][next(iter(bad['source_bindings']))]='0'*64
    name=tmp_path/'changed.json';name.write_text(packet.canonical(bad))
    with pytest.raises(ValueError):packet.verify(name)


def test_duplicate_keys_fail(tmp_path):
    name=tmp_path/'duplicate.json';name.write_text('{"schema":1,"schema":2}')
    with pytest.raises(ValueError):packet.verify(name)


def test_current_source_mutation_rejected(receipt,monkeypatch,tmp_path):
    name=tmp_path/'receipt.json';name.write_text(packet.canonical(receipt))
    changed=dict(receipt['source_bindings']);changed[next(iter(changed))]='f'*64
    monkeypatch.setattr(packet,'source_bindings',lambda:changed)
    with pytest.raises(ValueError):packet.verify(name)


def test_no_floating_diagnostics_retained(receipt):
    def visit(value):
        assert type(value) is not float
        if isinstance(value,dict):
            for x in value.values():visit(x)
        if isinstance(value,list):
            for x in value:visit(x)
    visit(receipt)


def test_loaded_local_import_closure_is_bound():
    code='''
from pathlib import Path
import sys
from model.native_fair_local_preparation import packet
paths=packet.source_bindings();missing=[]
for module in tuple(sys.modules.values()):
 name=getattr(module,'__file__',None)
 if not name:continue
 try:rel=Path(name).resolve().relative_to(packet.ROOT).as_posix()
 except ValueError:continue
 if rel.endswith('.py') and rel.startswith(('model/','tools/')) and rel not in paths:missing.append(rel)
assert not missing,missing
'''
    subprocess.run([sys.executable,'-c',code],cwd=packet.ROOT,check=True)
