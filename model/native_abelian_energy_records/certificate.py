"""Fresh source replay and exact finite arrival/energy certificate."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from . import model
from model.particle_common_contract.provenance import import_closure
PARENTS={
 'model/native_abelian_field_join/results.json':'f6aadcd79eec889de0c3243925c950e482bd456ea5c22fcab6a8406978472327',
 'model/affine_gauss_charged_clock/results.json':'7c8f0214f1412e38bdc410c607d0a8f8765df323c3a33a778da923f3892e8613'}

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def validate_pins(pins):
    for name,expected in pins.items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('repository-relative pin required')
        if digest(ROOT/relative)!=expected:raise ValueError('frozen source changed: '+name)

def raw_certificate():
    pins={}
    for path,expected in PARENTS.items():
        if digest(ROOT/path)!=expected:raise ValueError('frozen parent changed: '+path)
        value=json.loads((ROOT/path).read_text());validate_pins(value['input_sha256'])
        for name,sha in value['input_sha256'].items():
            if name in pins and pins[name]!=sha:raise ValueError('inconsistent parent pins')
            pins[name]=sha
        pins[path]=expected
    # Scientific parents are freshly regenerated, not admitted by hash alone.
    from model.native_abelian_field_join.certificate import certificate as join
    from model.affine_gauss_full_spectrum.certificate import certificate as spectrum
    for path,producer in (('model/native_abelian_field_join/results.json',join),
                          ('model/affine_gauss_full_spectrum/results.json',spectrum)):
        if producer()!=json.loads((ROOT/path).read_text()):raise ValueError('fresh parent science differs: '+path)
    result=model.calculate()
    result['prerequisite_receipts']=PARENTS
    result['fresh_replay']={'native_complete_reference_Gauss_and_whole_bin_controls':True,
      'full70928_actual_affine_H_residual':True,'all_new_bounds_exact_rational_or_outward_square_roots':True,
      'large_native_source_or_pointer_matrix_numerically_executed':False}
    result['capacity']={'retained_dimension_formula':'D^26 (2^31+1)^13 2^52 3^12 5456 2^416 Mprep^13 2^17 2^20',
      'native_D':'2^240 5^8 119^24','Mprep':1197721483502453176152084669,
      'energy_pointer_dimension':2**20,'arrival_pointer_reuses_parent_bit':True,
      'gate_synthesis_workspace_included':False}
    result['input_sha256']=dict(sorted(pins.items()))
    result['owned_source_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}
    return result

REFERENCE_SHA='e8a0a69f3162f8bfbdf8154ff84f6214ade453e801eeda1066f44df42a68805d'

def compare_research(value):
    path=HERE/'research_reference.json'
    if digest(path)!=REFERENCE_SHA:raise ValueError('immutable energy research receipt changed')
    expected=json.loads(path.read_text());validate_pins(expected['input_sha256'])
    ignored={'schema','input_sha256','owned_source_sha256'}
    if {k:v for k,v in value.items() if k not in ignored}!={k:v for k,v in expected.items() if k not in ignored}:
        raise ValueError('fresh energy science differs from immutable research result')
    return {'reference_sha256':REFERENCE_SHA,'all_scientific_fields_equal':True,
      'fresh_native_join_full_affine_residual_and_energy_bounds_replayed':True,
      'ignored_only_provenance_fields':sorted(ignored)}

def certificate():
    value=raw_certificate();comparison=compare_research(value)
    value.pop('owned_source_sha256');value['schema']='oph.native_abelian_energy_records.v1'
    seeds={p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.add((HERE/'research_reference.json').relative_to(ROOT).as_posix())
    for relative in import_closure(seeds):value['input_sha256'][relative]=digest(ROOT/relative)
    value['input_sha256']=dict(sorted(value['input_sha256'].items()))
    value['independent_research_comparison']=comparison
    value['arrival_readout_service']={
      'unitary':'Uread=(I-Qx) tensor I+Qx tensor X on admitted sector; identity elsewhere',
      'Hermitian_generator':'L=(pi/2) Qx tensor (I-X), extended by zero',
      'additional_interaction':'V(s)=exp(-is Hraw) L exp(+is Hraw), 0<=s<=1',
      'complete_evolution':'exp(-i Hraw) Uread; autonomous source runs on every branch',
      'duration_in_source_clock_units':1,'position_record_time':'start of arrival readout interval',
      'common_free_factor_commutes_through_subsequent_energy_instrument':True,
      'source_work_absolute_bound_on_admitted_sector':'15',
      'synthesis_and_controller_workspace_in_logical_pointer_count':False}
    return value
