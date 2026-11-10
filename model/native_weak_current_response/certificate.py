"""Fresh exact same-source dressed-current contraction and complete records."""
from pathlib import Path
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from . import bounds
from model.particle_common_contract.provenance import import_closure
PARENTS={
 'model/native_weak_current_interface/results.json':'db5ae5c7b9f43762a4cdb7128761d8f2014cf3c60d0c9f324569d27d1d7b3f66',
 'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json':'e135f7b15460427068dd18541b0f7990a3d259d298467604283841c68f75ab34',
 'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json':'7335d80394ab25ab949ca0e359e70a8eca55636256255b684c68086be2a8c7a6'}

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def validate(pins):
    for name,sha in pins.items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('repository-relative pin required')
        if digest(ROOT/name)!=sha:raise ValueError('frozen source changed: '+name)

def raw_certificate():
    pins={};saved={}
    for name,sha in PARENTS.items():
        if digest(ROOT/name)!=sha:raise ValueError('frozen prerequisite changed: '+name)
        value=json.loads((ROOT/name).read_text());validate(value['input_sha256']);saved[name]=value
        for path,expected in value['input_sha256'].items():
            if path in pins and pins[path]!=expected:raise ValueError('inconsistent source closure')
            pins[path]=expected
        pins[name]=sha
    from model.native_weak_current_interface.model import dressed_weak_density
    from model.coherent_defect_wave.chiral_joint.interacting_charged_response.certificate import certificate as charged
    from model.coherent_defect_wave.chiral_joint.interacting_ground_preparation.certificate import certificate as prep
    if dressed_weak_density()!=saved['model/native_weak_current_interface/results.json']['exact_Higgs_dressed_weak_current_interface']:
        raise ValueError('fresh original weak-current operator differs')
    for name,producer in ((list(PARENTS)[1],charged),(list(PARENTS)[2],prep)):
        if producer()!=saved[name]:raise ValueError('fresh full source/preparation differs: '+name)
    value=bounds.calculate();value['prerequisite_receipts']=PARENTS
    seeds=set(pins)|{'model/native_charge_resolution/bounds.py',
      'model/coherent_defect_wave/chiral_joint/atomic_momentum_response/model.py',
      'model/coherent_defect_wave/chiral_joint/atomic_momentum_response/PROOF.md'}
    for relative in import_closure(seeds):pins[relative]=digest(ROOT/relative)
    value['input_sha256']=dict(sorted(pins.items()))
    value['owned_source_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}
    value['fresh_replay']={'original_CAR_scalar_weak_operator_and_all_twelve_Gauss_commutators':True,
      'complete_local_source_and_all_magnetic_charged_target_images':True,
      'exact_full_source_charged_resolvent_and_error':True,
      'actual_physical_two_source_preparation_including_all_failures':True,
      'full_eigensolver_used':False,'large_full_source_pointer_matrix_executed':False}
    return value

RESEARCH_SHA='7a3b9756fc1a2e4db5d867492aef7a455676003ff4499293589714d4757f339d'

def compare_research(value):
    path=HERE/'research_reference.json'
    if digest(path)!=RESEARCH_SHA:raise ValueError('immutable weak-response research receipt changed')
    expected=json.loads(path.read_text());validate(expected['input_sha256'])
    ignored={'schema','input_sha256','owned_source_sha256'}
    science=lambda record:{k:v for k,v in record.items() if k not in ignored}
    if science(expected)!=science(value):raise ValueError('fresh weak-current science differs from research result')
    return {'reference_sha256':RESEARCH_SHA,'all_scientific_fields_equal':True,
      'fresh_complete_native_current_source_preparation_and_records_replayed':True,
      'ignored_only_provenance_fields':sorted(ignored)}

def certificate():
    value=raw_certificate();comparison=compare_research(value)
    value.pop('owned_source_sha256');value['schema']='oph.native_weak_current_response.v1'
    seeds={p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    seeds.add((HERE/'research_reference.json').relative_to(ROOT).as_posix())
    for relative in import_closure(seeds):value['input_sha256'][relative]=digest(ROOT/relative)
    value['input_sha256']=dict(sorted(value['input_sha256'].items()))
    value['independent_research_comparison']=comparison
    return value

def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');args=parser.parse_args()
    path=HERE/'results.json'
    if args.check:
        saved=json.loads(path.read_text());validate(saved['input_sha256'])
    value=certificate()
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    elif value!=saved:raise ValueError('fresh same-H current records differ')
    print('PASS: full native weak-current moments, inelastic finite records and ordered charged insertion')

if __name__=='__main__':main()
