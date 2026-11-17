"""Portable replay of independent finite Taylor and rational inertia controls."""
import copy,hashlib,json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.affine_gauss_spectrum.model import certificate as parent_certificate
from model.affine_gauss_mobile.certificate import validate_pins
from .model import audit
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
REFERENCE_SHA='e9d1180aab228045451fb245c9422ba7c9dccfbcb9e702aeebaa71c642ba5bef'
PARENT='model/affine_gauss_spectrum/results.json'
PARENT_SHA='632a3b6ae85c7f1fb67ce768a78b684972831eb1b194517b1dc9c2cb1affb9ff'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def compare_research(value):
    path=HERE/'research_reference.json'
    if digest(path)!=REFERENCE_SHA:raise ValueError('independent research receipt changed')
    old=json.loads(path.read_text())
    ignored={'schema','owned_source_sha256','input_sha256'}
    if {k:v for k,v in value.items() if k not in ignored}!={k:v for k,v in old.items() if k not in ignored}:
        raise ValueError('fresh independent validation differs from research evidence')
    return {'reference_sha256':REFERENCE_SHA,'all_scientific_fields_equal':True,
            'fresh_public_calculation_replayed':True,'ignored_only_provenance_fields':sorted(ignored)}

def certificate():
    if digest(ROOT/PARENT)!=PARENT_SHA:raise ValueError('frozen charged spectrum changed')
    parent=json.loads((ROOT/PARENT).read_text());validate_pins(parent['input_sha256'])
    if parent!=parent_certificate():raise ValueError('fresh charged spectrum replay differs')
    value=copy.deepcopy(audit());comparison=compare_research(value)
    value['schema']='oph.affine_gauss_spectrum_validation.v1'
    value.pop('owned_source_sha256')
    pins=dict(parent['input_sha256']);pins[PARENT]=PARENT_SHA
    seeds={p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.add((HERE/'research_reference.json').relative_to(ROOT).as_posix())
    for relative in import_closure(seeds):pins[relative]=digest(ROOT/relative)
    value['input_sha256']=dict(sorted(pins.items()))
    value['independent_research_comparison']=comparison
    return value
