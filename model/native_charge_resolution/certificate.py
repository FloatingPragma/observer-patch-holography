"""Exact root, full-source residual and refined charge-selection receipt."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from . import bounds as b
from model.particle_common_contract.provenance import import_closure
PARENT='model/native_abelian_field_join/results.json'
PARENT_SHA='f6aadcd79eec889de0c3243925c950e482bd456ea5c22fcab6a8406978472327'

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def validate(pins):
    for name,sha in pins.items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('repository-relative source required')
        if digest(ROOT/relative)!=sha:raise ValueError('changed native source: '+name)

def raw_certificate():
    if digest(ROOT/PARENT)!=PARENT_SHA:raise ValueError('frozen native join changed')
    parent=json.loads((ROOT/PARENT).read_text());validate(parent['input_sha256'])
    from model.native_abelian_field_join.certificate import certificate as replay_join
    if replay_join()!=parent:raise ValueError('fresh native source/reference/G1 replay differs')
    value=b.calculate();pins=dict(parent['input_sha256']);pins[PARENT]=PARENT_SHA
    value['input_sha256']=dict(sorted(pins.items()))
    value['owned_source_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}
    value['fresh_replay']={'full_native_reference_Gauss_and_G1_source':True,
      'actual_local_CAR_target_annihilator_and_positive_weights':True,
      'neutral_and_all9_charged_roots_refined_exactly':True,
      'full_source_residual_weight_transfer_and_neutral_shift_exact':True,
      'full_native_eigensolver_used':False}
    return value

REFERENCE_SHA='9a88dc18a69a9d12eb7058d1bd62b42f015f7e0531147d3a885e82559ab6a086'

def compare_research(value):
    path=HERE/'research_reference.json'
    if digest(path)!=REFERENCE_SHA:raise ValueError('immutable resolution research receipt changed')
    old=json.loads(path.read_text());validate(old['input_sha256'])
    ignored={'schema','input_sha256','owned_source_sha256'}
    if {k:v for k,v in old.items() if k not in ignored}!={k:v for k,v in value.items() if k not in ignored}:
        raise ValueError('fresh resolution science differs from research result')
    return {'reference_sha256':REFERENCE_SHA,'all_scientific_fields_equal':True,
      'fresh_native_source_local_resolvent_roots_and_full_source_bounds_replayed':True,
      'ignored_only_provenance_fields':sorted(ignored)}

def certificate():
    value=raw_certificate();comparison=compare_research(value)
    value.pop('owned_source_sha256');value['schema']='oph.native_charge_resolution.v1'
    seeds={p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.add((HERE/'research_reference.json').relative_to(ROOT).as_posix())
    for relative in import_closure(seeds):value['input_sha256'][relative]=digest(ROOT/relative)
    value['input_sha256']=dict(sorted(value['input_sha256'].items()))
    value['independent_research_comparison']=comparison
    value['construction_scope']={'readout_source':'original full native h plus full30-mode static CAR reference with energy4N',
      'outer_affine_coupled_H_is_the_readout_source':False,
      'full_opposite_charge_cost_equality_requires_declared_dual_reference':True,
      'frozen_native_constituent_and_all_existing_outer_couplers_changed':False}
    return value
