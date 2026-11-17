"""Fresh structural source/reference/record proof with immutable prerequisites."""
from pathlib import Path
import argparse,hashlib,json,sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from .model import calculate
from model.particle_common_contract.provenance import import_closure
PARENTS=(
 'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json',
 'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json',
 'model/coherent_defect_wave/chiral_joint/interacting_charge_conjugation/results.json',
 'model/affine_gauss_mobile/results.json',
 'model/affine_gauss_spectrum/results.json',
 'model/affine_gauss_full_spectrum/results.json')

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def raw_certificate():
    pins={};receipts={}
    for relative in PARENTS:
        path=ROOT/relative;value=json.loads(path.read_text())
        for source,expected in value['input_sha256'].items():
            if digest(ROOT/source)!=expected:raise ValueError('changed frozen source: '+source)
            if source in pins and pins[source]!=expected:raise ValueError('inconsistent inherited source')
            pins[source]=expected
        pins[relative]=digest(path);receipts[relative]=digest(path)
    value=calculate()
    seeds={'model/coherent_defect_wave/physical_attachment/chiral.py',
      'model/coherent_defect_wave/chiral_source/model.py',
      'model/coherent_defect_wave/chiral_joint/oriented_links/rotations.py',
      'model/particle_common_contract/provenance.py'}
    for relative in import_closure(seeds):pins[relative]=digest(ROOT/relative)
    value['prerequisite_receipts']=receipts
    value['input_sha256']=dict(sorted(pins.items()))
    value['owned_source_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}
    value['scope_boundaries']={'declared_new_Hamiltonian_and_enlarged_reference':True,
      'original_uniform_native_lattice_H_preserved':False,
      'each_original_intrablock_H_preserved':True,
      'whole_bin_native_H_replaced_by_center':False,
      'entire_native_complement_retained':True,
      'selected_native_eigenvectors_or_full_native_source_evolution_numerically_executed':False,
      'finite_spectral_and_state_synthesis_services_declared':True,
      'efficient_native_gate_implementation_proved':False,
      'unpenalized_full_charged_ground_identified':False,
      'penalty_alternative_is_separate_Hamiltonian':True,
      'continuous_physical_spin_Dirac_or_SM_identification':False}
    return value

REFERENCE_SHA='c16970e826229c2c3489962ddddc5f2ff439f5d4e0f67d7802147f75e8cf0170'

def validate_pins(pins):
    for name,expected in pins.items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('repository-relative source pin required')
        if digest(ROOT/relative)!=expected:raise ValueError('source custody mismatch: '+name)

def compare_research(value):
    path=HERE/'research_reference.json'
    if digest(path)!=REFERENCE_SHA:raise ValueError('immutable native join research receipt changed')
    old=json.loads(path.read_text());validate_pins(old['input_sha256'])
    ignored={'schema','input_sha256','owned_source_sha256'}
    if {k:v for k,v in value.items() if k not in ignored}!={k:v for k,v in old.items() if k not in ignored}:
        raise ValueError('fresh native join differs from independent research result')
    return {'reference_sha256':REFERENCE_SHA,'all_scientific_fields_equal':True,
      'fresh_source_reference_Gauss_interval_and_record_controls_replayed':True,
      'ignored_only_provenance_fields':sorted(ignored)}

def certificate():
    value=raw_certificate();comparison=compare_research(value)
    value.pop('owned_source_sha256')
    value['schema']='oph.native_abelian_field_join.v1'
    seeds={p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.add((HERE/'research_reference.json').relative_to(ROOT).as_posix())
    for relative in import_closure(seeds):value['input_sha256'][relative]=digest(ROOT/relative)
    value['input_sha256']=dict(sorted(value['input_sha256'].items()))
    value['independent_research_comparison']=comparison
    return value
