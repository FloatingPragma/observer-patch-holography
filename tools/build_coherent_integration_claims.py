"""Bind the exact shared-hardware interfaces without changing physical statuses."""
from pathlib import Path
from hashlib import sha256
import argparse,importlib,json
from tools.coherent_integration_results import NAMES,TITLES,project
ROOT=Path(__file__).resolve().parents[1]
def generate():
 data={};claims=[]
 for name,title in zip(NAMES,TITLES,strict=True):
  path=f'model/openwave_{name}/results.json';value=json.loads((ROOT/path).read_text());owner=importlib.import_module(f'model.openwave_{name}.__main__')
  assert value['source_sha256']==owner.pins() and value['verdict'].startswith('PASS'),name
  data[name]=value;claims.append({'id':'coherent_integration_'+name,'label':title,'status':'verified_in_stated_finite_domain','receipt':path,'receipt_sha256':sha256((ROOT/path).read_bytes()).hexdigest(),'proof':f'model/openwave_{name}/PROOF.md','independent_Nature_match_added':False,'physical_rest_mass_identification_added':False,'official_OpenWave_award_added':False})
 p=project(data);paths=('tracking/premise_register.json','claims/physical_identification_registry.json','claims/selection_ledger.json')
 return {'schema':'oph.coherent_native_integration_claims.v1','mission':'One declared executable construction reproduces all compatible physical observer records.','architecture_uniqueness_required':False,'claims':claims,'canonical_source_sha256':{x:sha256((ROOT/x).read_bytes()).hexdigest()for x in paths},'exact_interfaces':p['exact_interfaces'],'shared_declared_hardware':p['shared_declared_hardware'],'named_programs':p['named_programs'],'uncompensated_clock_counterexample':p['uncompensated_clock_counterexample'],'joint_spatial_charge_clock_interacting_field_implemented':False,'scope_of_joint_flag':'These seven named component programs only; the separate joint observer register specifies the declared charged-field extension and complete tracked-clock readouts.','physical_premises_closed':[],'new_independent_GeV_rest_masses':[],'complete_physical_mission_proved':False,'official_status':{'validated':0,'untested':31},'next_acceptance':'Establish a shared stationary particle generator and physical identification on the joint charged-field construction; retain complete records, calibration ancestry and continuum obligations.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();text=json.dumps(generate(),sort_keys=True,indent=2,ensure_ascii=False)+'\n';out=ROOT/'claims/coherent_integration.json'
 if a.write:out.write_text(text)
 else:assert out.read_text()==text
 print('PASS_COHERENT_INTEGRATION_WITH_EXACT_INTERFACES_AND_PHYSICAL_BOUNDARIES')
