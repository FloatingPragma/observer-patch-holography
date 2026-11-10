"""Fresh native preparation/resolution and finite-pointer centroid receipt."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from . import bounds as b
from model.particle_common_contract.provenance import import_closure
PARENTS={
 'model/native_charge_resolution/results.json':'f22c9a266a43de8aaae776b1293894081c3dc20b7d8e5678e30ef6af993dd19e',
 'model/affine_gauss_charged_clock/results.json':'7c8f0214f1412e38bdc410c607d0a8f8765df323c3a33a778da923f3892e8613'}
PREPARATION='model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json'

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def validate(pins):
    for name,sha in pins.items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('repository-relative source pin required')
        if digest(ROOT/relative)!=sha:raise ValueError('changed frozen source: '+name)

def raw_certificate():
    pins={}
    for name,sha in PARENTS.items():
        if digest(ROOT/name)!=sha:raise ValueError('frozen parent changed: '+name)
        value=json.loads((ROOT/name).read_text());validate(value['input_sha256'])
        for path,expected in value['input_sha256'].items():
            if path in pins and pins[path]!=expected:raise ValueError('inconsistent source closure')
            pins[path]=expected
        pins[name]=sha
    from model.native_charge_resolution.certificate import certificate as resolution
    from model.coherent_defect_wave.chiral_joint.interacting_ground_preparation.certificate import certificate as preparation
    for name,producer in (('model/native_charge_resolution/results.json',resolution),(PREPARATION,preparation)):
        if producer()!=json.loads((ROOT/name).read_text()):raise ValueError('fresh native source/preparation differs: '+name)
    value=b.calculate()
    value['prerequisite_receipts']={**PARENTS,PREPARATION:digest(ROOT/PREPARATION)}
    value['fresh_replay']={'actual_physical_two_copy_preparation_including_failure_and_trace_error':True,
      'native_complete_local_resolvent_and_refined_actual_H_bin':True,
      'all59_complete_PW119_blank_matrix_elements':True,
      'centroid_conditioning_and_finite_pointer_errors_exact':True,
      'clock_parent_only_source_independent_identities_used':True,
      'affine_state_energy_probability_or_H_transferred':False,
      'full_native_eigensolver_or_source_pointer_matrix_execution':False}
    value['input_sha256']=dict(sorted(pins.items()))
    value['owned_source_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}
    return value

REFERENCE_SHA='412146a491d580c14280b77739a4d997d04e24b021c1b8fd3d09982786c7faa4'

def compare_research(value):
    path=HERE/'research_reference.json'
    if digest(path)!=REFERENCE_SHA:raise ValueError('immutable centroid research receipt changed')
    expected=json.loads(path.read_text());validate(expected['input_sha256'])
    ignored={'schema','input_sha256','owned_source_sha256'}
    if {k:v for k,v in expected.items() if k not in ignored}!={k:v for k,v in value.items() if k not in ignored}:
        raise ValueError('fresh centroid science differs from research result')
    return {'reference_sha256':REFERENCE_SHA,'all_scientific_fields_equal':True,
      'fresh_actual_preparation_native_resolution_and_finite_mean_bounds_replayed':True,
      'ignored_only_provenance_fields':sorted(ignored)}

def certificate():
    value=raw_certificate();comparison=compare_research(value)
    value.pop('owned_source_sha256');value['schema']='oph.native_charge_centroid.v1'
    seeds={p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.add((HERE/'research_reference.json').relative_to(ROOT).as_posix())
    for relative in import_closure(seeds):value['input_sha256'][relative]=digest(ROOT/relative)
    value['input_sha256']=dict(sorted(value['input_sha256'].items()))
    value['independent_research_comparison']=comparison
    value['construction_scope']={
      'autonomous_source':'two complete original eight-cell h copies plus full30-mode static reference4N; preparation references retained',
      'outer_affine_coupled_source_used':False,
      'actual_prepared_source_state_asserted_exact_or_stationary':False,
      'global_ground_projection_added':False,
      'prior_native_bin_admission_is_declared_spectral_service':True,
      'efficient_native_gate_synthesis_or_controller_workspace_proved':False}
    return value
