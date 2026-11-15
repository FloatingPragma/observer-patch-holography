"""Bind ten scoped native results without promoting physical assumptions."""
from pathlib import Path
from hashlib import sha256
import argparse,importlib,json
from tools.native_execution_results import NAMES,TITLES,SHORTCUTS,project
ROOT=Path(__file__).resolve().parents[1]

def generate():
    data={};claims=[]
    for index,(name,title,shortcut)in enumerate(zip(NAMES,TITLES,SHORTCUTS,strict=True),1):
        prefix='model/openwave_'+name;receipt=ROOT/prefix/'results.json'
        value=json.loads(receipt.read_text());owner=importlib.import_module('model.openwave_'+name+'.__main__')
        assert value['source_sha256']==owner.pins() and value['verdict'].startswith('PASS'),name
        data[name]=value
        claims.append({'id':'native_execution_'+name,'task':index,'label':title,
            'status':'verified_in_stated_finite_domain','finite_theorem_shortcut':shortcut,
            'receipt':prefix+'/results.json','receipt_sha256':sha256(receipt.read_bytes()).hexdigest(),
            'proof':prefix+'/PROOF.md',
            'construction_relation':'literal repair process and explicit declared coherent-service boundary',
            'independent_Nature_match_added':False,'physical_particle_mass_identification_added':False,
            'official_criterion_awarded_by_this_work':False})
    project(data)
    canonical=('tracking/premise_register.json','claims/physical_identification_registry.json','claims/selection_ledger.json')
    return {'schema':'oph.native_execution_claims.v1',
        'mission':'One specified executable construction reproduces compatible observer records for all physical observables.',
        'architecture_uniqueness_required':False,'claims':claims,
        'canonical_source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest()for p in canonical},
        'physical_premises_closed':[],'independent_physical_masses_GeV_added':[],
        'official_status':data['native_platform_gate']['official_status'],
        'complete_physical_mission_proved':False,
        'next_acceptance':'Declare and execute the coherent memory/charge extension needed by the intended observer; prove its common-source operator and record transfer before claiming a physical field, mass or continuum.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    text=json.dumps(generate(),sort_keys=True,indent=2,ensure_ascii=False)+'\n';out=ROOT/'claims/native_execution.json'
    if a.write:out.write_text(text)
    else:assert out.read_text()==text
    print('PASS_TEN_NATIVE_EXECUTION_CLAIMS_WITH_SCOPED_PHYSICAL_BOUNDARY')
