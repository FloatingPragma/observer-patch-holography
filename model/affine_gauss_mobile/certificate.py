"""Fresh source replay, immutable ancestry and independent scientific comparison."""
import copy
import hashlib
import json
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from .model import calculate

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
REFERENCE=HERE/'research_reference.json'
REFERENCE_SHA='79de12b140bf5e4ce8335ad1eef97bc851882ebd11d96f6a2152dc6d79bcdd17'
IGNORED_COMPARISON_FIELDS={'schema','input_sha256','owned_source_sha256','parent_receipt_sha256'}


def validate_pins(pins):
    for name,expected in pins.items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('source custody requires repository-relative paths')
        path=ROOT/relative
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('source custody mismatch: '+name)


def scientific_fields(value):
    return {key:item for key,item in value.items() if key not in IGNORED_COMPARISON_FIELDS}


def compare_research(value):
    if hashlib.sha256(REFERENCE.read_bytes()).hexdigest()!=REFERENCE_SHA:
        raise ValueError('immutable independent research receipt changed')
    reference=json.loads(REFERENCE.read_text())
    if scientific_fields(value)!=scientific_fields(reference):
        raise ValueError('fresh public scientific calculation differs from independent research receipt')
    return {'reference_sha256':REFERENCE_SHA,'all_scientific_fields_equal':True,
            'ignored_only_provenance_fields':sorted(IGNORED_COMPARISON_FIELDS),
            'fresh_calculation_was_replayed':True}


def input_hashes():
    seeds={str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py','.md')}
    seeds.update({str(REFERENCE.relative_to(ROOT)),'model/requirements.txt'})
    from .model import source,PARENT
    parent,parent_sha=source()
    validate_pins(parent['input_sha256'])
    seeds.update(parent['input_sha256'])
    seeds.add(PARENT)
    names=import_closure(seeds)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def certificate():
    pins=input_hashes()
    value=copy.deepcopy(calculate())
    comparison=compare_research(value)
    value.pop('owned_source_sha256',None)
    value['schema']='oph.affine_gauss_mobile.v1'
    value['input_sha256']=pins
    value['independent_research_comparison']=comparison
    return value


def write():
    value=certificate()
    (HERE/'results.json').write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    return value
