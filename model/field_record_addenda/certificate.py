"""Strict archive comparison and source closure for additive public field records."""
from __future__ import annotations
import argparse
from copy import deepcopy
import hashlib
import importlib
import json
from pathlib import Path
from .registry import PACKAGES
from .provenance import ROOT,import_closure
from tools import strict_json

METADATA=frozenset(('input_sha256','owned_source_sha256'))
HERE=Path(__file__).resolve().parent


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scientific(value):return {k:v for k,v in value.items()if k not in METADATA}


def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)


def compare_science(expected,actual):
    if canonical(scientific(expected))!=canonical(scientific(actual)):
        raise ValueError('public port changed a scientific field or field inventory')


def authenticate_archive(name,here):
    if name not in PACKAGES:raise ValueError('unregistered additive package')
    here=Path(here);path=here/'research_reference.json'
    authority=PACKAGES[name]
    if digest(path)!=authority['research_sha256']:
        raise ValueError('immutable research receipt changed')
    retained=strict_json.load(path)
    expected={filename+'.raw':sha for filename,sha in retained['owned_source_sha256'].items()}
    if {p.name for p in (here/'retained').iterdir()if p.is_file()}!=set(expected):
        raise ValueError('immutable research source archive inventory changed')
    for filename,sha in expected.items():
        if digest(here/'retained'/filename)!=sha:
            raise ValueError('immutable original producer or proof changed: '+filename)
    return retained


def authenticate_adaptations(name,here,retained):
    here=Path(here);value=strict_json.load(here/'port_adaptations.json')
    if (set(value)!={'research_receipt_sha256','files','scientific_payload_changes'}
            or value['research_receipt_sha256']!=PACKAGES[name]['research_sha256']
            or value['scientific_payload_changes']!=[]):
        raise ValueError('reviewed import/path adaptation scope changed')
    rows=value['files'];old=retained['owned_source_sha256']
    if len(rows)!=len(old)or {r['original']for r in rows}!=set(old):
        raise ValueError('reviewed original source mapping changed')
    for row in rows:
        original=row['original'];target='source_producer.py'if original=='producer.py'else original
        if (set(row)!={'original','public','retained_sha256','adaptations'}
                or row['public']!=target or row['retained_sha256']!=old[original]
                or not(here/target).is_file()):
            raise ValueError('reviewed original source path changed')


def authenticate_external(name,root=ROOT):
    seeds=set()
    for path,sha in PACKAGES[name]['external_receipts'].items():
        if Path(path).is_absolute()or '..'in Path(path).parts:
            raise ValueError('external source path escapes the public repository')
        if digest(root/path)!=sha:raise ValueError('fixed external source receipt changed: '+path)
        value=strict_json.load(root/path);seeds.add(path)
        for source,expected in value['input_sha256'].items():
            if Path(source).is_absolute()or '..'in Path(source).parts:
                raise ValueError('inherited source path escapes the public repository')
            if digest(root/source)!=expected:raise ValueError('fixed external source closure changed: '+source)
            seeds.add(source)
    return seeds


def public_certificate(name,calculate,here):
    here=Path(here).resolve();retained=authenticate_archive(name,here)
    authenticate_adaptations(name,here,retained)
    dependency_inputs=authenticate_external(name)
    for dependency in PACKAGES[name]['dependencies']:
        module=importlib.import_module('model.field_record_addenda.'+dependency+'.producer')
        fresh=module.certificate();saved=strict_json.load(module.HERE/'results.json')
        if canonical(fresh)!=canonical(saved):raise ValueError('sibling public producer differs: '+dependency)
        dependency_inputs.update(saved['input_sha256'])
        dependency_inputs.add(str((module.HERE/'results.json').relative_to(ROOT)))
    value=json.loads(json.dumps(calculate(),default=str,allow_nan=False))
    compare_science(retained,value)
    for path,sha in retained.get('input_sha256',{}).items():
        if digest(ROOT/path)!=sha:raise ValueError('immutable research source ancestry changed: '+path)
    seeds=set(retained.get('input_sha256',{}))|set(value.get('input_sha256',{}))|dependency_inputs
    seeds.update(str(p.relative_to(ROOT))for p in here.iterdir()if p.suffix in ('.py','.md'))
    seeds.update(str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in ('.py','.md'))
    seeds.update(str(p.relative_to(ROOT))for p in (here/'retained').iterdir()if p.is_file())
    seeds.update(str((here/p).relative_to(ROOT))for p in ('research_reference.json','port_adaptations.json'))
    value['input_sha256']={path:digest(ROOT/path)for path in sorted(import_closure(seeds))}
    value['owned_source_sha256']={p.name:digest(p)for p in sorted(here.iterdir())if p.suffix in ('.py','.md')}
    value['public_port']={'package':'model/field_record_addenda/'+name,
        'retained_research_receipt':'research_reference.json',
        'retained_research_sha256':PACKAGES[name]['research_sha256'],
        'retained_original_sources':'retained/*.raw; exact original owned-source hashes authenticated',
        'exact_scientific_payload_equal':True,
        'ignored_comparison_fields':sorted(METADATA),
        'comparison_stage':'calculated owner payload, before public_port metadata is appended',
        'adaptation_manifest':'port_adaptations.json',
        'repository_local_runtime':True,'physical_status_promotion':False}
    return value


def run_cli(certificate,here):
    p=argparse.ArgumentParser();m=p.add_mutually_exclusive_group(required=True)
    m.add_argument('--write',action='store_true');m.add_argument('--check',action='store_true')
    args=p.parse_args();here=Path(here);value=certificate();target=here/'results.json'
    if args.write:target.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')
    elif canonical(value)!=canonical(strict_json.load(target)):
        raise ValueError('public record, complete source inventory or live source pins differ')
    print('PASS: '+here.name+'; exact retained scientific result and complete public source custody')
