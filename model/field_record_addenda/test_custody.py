"""Every retained test is collected; archives, scientific boundaries and routing fail closed."""
from copy import deepcopy
import ast
import importlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import os
import pytest
from . import certificate as c
from .registry import PACKAGES
from .provenance import ROOT


@pytest.mark.parametrize('name',tuple(PACKAGES))
def test_immutable_research_receipt_and_original_owned_sources(name):
    here=c.HERE/name;value=c.authenticate_archive(name,here)
    assert c.digest(here/'research_reference.json')==PACKAGES[name]['research_sha256']
    assert value['owned_source_sha256']
    c.authenticate_adaptations(name,here,value)
    assert set(PACKAGES[name]['external_receipts'])<=c.authenticate_external(name)


@pytest.mark.parametrize('name',tuple(PACKAGES))
def test_every_original_test_function_survives_public_collection(name):
    here=c.HERE/name;value=c.authenticate_archive(name,here)
    def functions(text):
        return sorted(n.name for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))and n.name.startswith('test_'))
    count=0
    for filename in value['owned_source_sha256']:
        if filename.startswith('test_')and filename.endswith('.py'):
            old=functions((here/'retained'/(filename+'.raw')).read_text())
            assert old==functions((here/filename).read_text())
            count+=len(old)
    assert count>=4


@pytest.mark.parametrize('name',tuple(PACKAGES))
def test_any_scientific_change_extra_field_or_deleted_field_is_rejected(name):
    old=c.authenticate_archive(name,c.HERE/name)
    key=next(k for k in old if k not in c.METADATA and k!='schema')
    changed=deepcopy(old);changed[key]={'unreviewed':True}
    with pytest.raises(ValueError,match='scientific'):c.compare_science(old,changed)
    changed=deepcopy(old);del changed[key]
    with pytest.raises(ValueError,match='scientific'):c.compare_science(old,changed)
    changed=deepcopy(old);changed['unlisted_physical_status_promotion']=True
    with pytest.raises(ValueError,match='scientific'):c.compare_science(old,changed)
    changed=deepcopy(old);changed['input_sha256']={};changed['owned_source_sha256']={}
    c.compare_science(old,changed)  # Custody is checked separately against actual files.


@pytest.mark.parametrize('mutation',('receipt','original_code','original_proof','extra_archive','manifest'))
def test_archive_or_reviewed_adaptation_mutation_is_rejected(mutation,tmp_path):
    name='native_orthogonal_wilson_bridge';here=tmp_path/name
    shutil.copytree(c.HERE/name,here)
    retained=c.authenticate_archive(name,here)
    if mutation=='receipt':(here/'research_reference.json').write_text('{}\n')
    elif mutation=='extra_archive':(here/'retained'/'extra.py.raw').write_text('pass\n')
    elif mutation=='manifest':
        data=json.loads((here/'port_adaptations.json').read_text());data['scientific_payload_changes']=['refit']
        (here/'port_adaptations.json').write_text(json.dumps(data))
    else:
        filename='bridge.py.raw'if mutation=='original_code'else'PROOF.md.raw'
        with (here/'retained'/filename).open('a')as f:f.write('\n')
    with pytest.raises(ValueError):
        record=c.authenticate_archive(name,here);c.authenticate_adaptations(name,here,record)


def test_no_private_runtime_resolver_or_generic_module_side_loading():
    for path in c.HERE.rglob('*.py'):
        if path==Path(__file__).resolve():continue
        source=path.read_text();assert 'sys.path.insert'not in source
        assert "spec_from_file_location"not in source
        assert 'reverse-engineering-reality'not in source
        assert 'CODEX_JOIN_REFINEMENT'not in source
    assert set(c.METADATA)=={'input_sha256','owned_source_sha256'}
    assert sum(bool(v['dependencies'])for v in PACKAGES.values())==6


def test_all_eleven_cold_public_producers_and_closed_manifests():
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    result=subprocess.run([sys.executable,'-B','-m','model.field_record_addenda','--check'],
                          cwd=ROOT,env=env,capture_output=True,text=True,check=False)
    assert result.returncode==0,result.stdout+result.stderr
    for name in PACKAGES:
        here=c.HERE/name;value=json.loads((here/'results.json').read_text())
        assert value['public_port']['exact_scientific_payload_equal']
        assert not value['public_port']['physical_status_promotion']
        for path,sha in value['input_sha256'].items():assert c.digest(ROOT/path)==sha,path
        assert str((here/'research_reference.json').relative_to(ROOT))in value['input_sha256']
        assert str((here/'port_adaptations.json').relative_to(ROOT))in value['input_sha256']
        assert set(PACKAGES[name]['external_receipts'])<=set(value['input_sha256'])


def test_changed_declared_external_receipt_is_rejected(monkeypatch):
    name='native_wilson_cone_audit'
    path=next(iter(PACKAGES[name]['external_receipts']))
    monkeypatch.setitem(PACKAGES[name]['external_receipts'],path,'0'*64)
    with pytest.raises(ValueError,match='fixed external source receipt'):
        c.authenticate_external(name)
