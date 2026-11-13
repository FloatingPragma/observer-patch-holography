"""Bind actual joint native producers without changing physical premise statuses."""
from pathlib import Path
from hashlib import sha256
import argparse,importlib,json
from tools.joint_observer_results import NAMES,TITLES,project
ROOT=Path(__file__).resolve().parents[1]
def generate():
    data={};claims=[]
    for name,title in zip(NAMES,TITLES,strict=True):
        path=f'model/{name}/results.json';value=json.loads((ROOT/path).read_text())
        owner=importlib.import_module(f'model.{name}.__main__')
        if value.get('source_sha256')!=owner.pins()or not value.get('verdict','').startswith('PASS'):
            raise ValueError('joint source custody or acceptance changed: '+name)
        data[name]=value
        claims.append({'id':'joint_observer_'+name,'label':title,'status':'verified_in_stated_finite_domain',
          'receipt':path,'receipt_sha256':sha256((ROOT/path).read_bytes()).hexdigest(),'proof':f'model/{name}/PROOF.md',
          'physical_rest_mass_identification_added':False,'independent_Nature_match_added':False,'official_OpenWave_award_added':False})
    p=project(data)
    registers=('tracking/premise_register.json','claims/physical_identification_registry.json','claims/selection_ledger.json')
    return {'schema':'oph.joint_native_observer_claims.v1',
      'mission':'One declared executable construction reproduces all compatible physical observer records.',
      'architecture_uniqueness_required':False,'claims':claims,
      'canonical_source_sha256':{x:sha256((ROOT/x).read_bytes()).hexdigest()for x in registers},
      'common_program_sha256':p['common_program_sha256'],'shared_integer_extension':p['shared_integer_extension'],
      'proved_finite_interfaces':p['proved_finite_interfaces'],'domains':p['domains'],
      'physical_premises_closed':[],'physical_identification_ids_closed':[],
      'independent_Nature_matches_added':0,'new_independent_GeV_rest_masses':[],
      'complete_physical_mission_proved':False,'official_status':p['official_status'],
      'next_acceptance':'Identify physical particles and their rest spectrum on the joint charged field source, close the local clock and physical interaction normalization, and transfer the complete interacting observer program to its declared spatial refinement.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    text=json.dumps(generate(),sort_keys=True,indent=2,ensure_ascii=False)+'\n';out=ROOT/'claims/joint_observer.json'
    if a.write:out.write_text(text)
    else:assert out.read_text()==text
    print('PASS_JOINT_NATIVE_OBSERVER_CLAIMS_WITH_PHYSICAL_BOUNDARIES')
