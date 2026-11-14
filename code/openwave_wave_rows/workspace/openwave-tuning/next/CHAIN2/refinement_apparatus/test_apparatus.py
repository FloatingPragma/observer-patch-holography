"""Executable development tests and formal resource/error evaluations."""
from fractions import Fraction as F
from pathlib import Path
import copy,hashlib,json,math
import numpy as np
import apparatus as b
HERE=Path(__file__).resolve().parent

def require_failure(call):
 try:call()
 except (ValueError,ArithmeticError):return True
 raise AssertionError('invalid input or incomplete law admitted')

def main():
 pins=b.check_pins();rows=[];natives=[]
 for size in (2,4):
  q,d=b.member(size);initial=np.zeros((1<<d,2),complex);initial[0,0]=1;initial[1<<(8+2*q),1]=1
  run=b.complete_run(b.preparation_template(size),d,initial);b.require_complete(run)
  target=b.expected_isometry(size);error=float(np.linalg.norm(run['success']-target,2))
  gram=float(np.linalg.norm(run['success'].conj().T@run['success']+sum(run['failure_grams'])-np.eye(2),2))
  assert error<1e-11 and gram<1e-11
  # Every two-column superposition and entangled reference inherits this norm.
  entangled=run['success']/math.sqrt(2)
  assert abs(np.linalg.norm(entangled)**2-1)<1e-11
  # Native source gate word on all clean data/law columns, no graph above4.
  program=b.a.reversible_program(size);memory=program['layout']['memory_bits'];law=program['layout']['law']
  columns=list(range(1<<d))+[x|(1<<law)for x in range(1<<d)]
  outputs=b.a.execute_bitplanes(program['gates'],columns,memory)
  expected=np.array([program['routing_permutation'][x&((1<<d)-1)]|(x&(1<<law))for x in columns])
  assert np.array_equal(outputs,expected)
  natives.append({'size':size,'clean_data_law_columns':len(columns),'reversible_gates':len(program['gates']),
                  'all_column_outputs_match':True,'memory_bits':memory})
  source=b.a.certificate(size,F(1,size*size))
  settings={}
  for setting in ('plus','minus','plus_i','minus_i'):
   prepared=b.prepare(size,setting);evolved=b.compressed_walk(prepared['success'],size,size)
   probability=np.abs(evolved[:,0])**2
   continuum=b.continuum_state(size,size,setting);comparison=np.abs(continuum)**2
   density=probability.reshape(size**3,256).sum(axis=1)
   mean=probability@b.current_table(size)
   initial_density=np.abs(prepared['success'][:,0]).reshape(size**3,256)**2
   laws={}
   for law in ('C','E'):
    terminal=b.complete_address_law(probability,law)
    total=float(sum(terminal['address_probabilities'])+sum(terminal['first_detector_export_failures']))
    assert abs(total-sum(probability))<1e-10
    laws[law]={'address_total':float(sum(terminal['address_probabilities'])),
               'first_detector_export_failures':terminal['first_detector_export_failures'],
               'fine_records_per_bit':[m for m,v in terminal['fine_record_factorization']],
               'density':terminal['address_probabilities'].reshape(size**3,256).sum(axis=1).tolist(),
               'current_mean':(terminal['address_probabilities']@b.current_table(size)).tolist(),
               'padded_output_probability':float(terminal['address_probabilities'].reshape(-1,16,16)[:,12:,:].sum()+terminal['address_probabilities'].reshape(-1,16,16)[:,:12,12:].sum())}
   change=float(np.sum(abs(density-initial_density.sum(axis=1)))/2)
   settings[setting]={'density':density.tolist(),'density_change_from_preparation_TV':change,'current_mean':mean.tolist(),
      'full_address_continuum_comparison_TV_diagnostic':float(np.sum(abs(probability-comparison))/2),
      'continuum_theorem_applies_to_this_small_member':False,'laws':laws,
      'preparation_export_failure_probability':float(sum(np.trace(g).real for g in prepared['failure_grams']))}
  # Phase convention is checked against the frozen144 Bloch source on the same small member.
  one=b.compressed_walk(target[:,[1]],size,1)[:,0]
  mode=b.continuum.bloch_walk((2*np.pi/size,0,0))@(np.ones(144)/12)
  expected_mode=np.zeros(1<<d,complex)
  for cell in range(size**3):
   x=cell//size**2
   for port in range(12):expected_mode[256*cell+16*port:256*cell+16*port+12]=np.exp(2j*np.pi*x/size)*mode[12*port:12*port+12]/math.sqrt(size**3)
  assert np.linalg.norm(one-expected_mode)<1e-11
  for left,right in (('plus','minus'),('plus_i','minus_i')):
   # Equal native selector phase controls implement input dephasing when the phase flag is forgotten.
   mixture=(np.array(settings[left]['density'])+np.array(settings[right]['density']))/2
   assert np.max(abs(mixture-1/size**3))<1e-11
  assert max(x['density_change_from_preparation_TV']for x in settings.values())>.01
  assert np.linalg.norm(np.array(settings['plus']['density'])-settings['minus']['density'])>.01
  if size==4:assert np.linalg.norm(np.array(settings['plus_i']['density'])-settings['minus_i']['density'])>.01
  rows.append({'size':size,'steps':size,'h_exact':str(F(1,size)),'time_exact':'1','K_bound':7,
               'hK_window_satisfied':False,'preparation_isometry_error_diagnostic':error,
               'complete_gram_error_diagnostic':gram,'source_certificate':source,'emitted_program_budget':b.development_budget(size,size),'settings':settings})
  print('development member',size,'passed',flush=True)
 # Complete M6 small-template tests with every export complement.
 templates={**b.a.primitive_templates(),'reflection':b.a.reflection_template(),'uniform':b.a.uniform_template(),'pointer':b.detector_template()}
 local=[]
 for name,word in templates.items():
  bits=max(t for k,ts,v in word for t in ts)+1
  run=b.complete_run(word,bits,np.eye(1<<bits));b.require_complete(run)
  gram=np.linalg.norm(run['success'].conj().T@run['success']+sum(run['failure_grams'])-np.eye(1<<bits),2)
  assert gram<1e-11
  actual=run['success'];expected=None
  if name=='flip':expected=np.array([[0,1],[1,0]])
  elif name=='hadamard':expected=np.array([[1,1],[1,-1]])/math.sqrt(2)
  elif name in ('cx','ccx'):
   perm=list(range(1<<bits));left,right=(1,3)if name=='cx'else(3,7);perm[left],perm[right]=perm[right],perm[left];expected=np.eye(1<<bits)[:,perm]
  elif name=='reflection':
   u=np.r_[np.ones(12)/math.sqrt(12),np.zeros(4)];expected=2*np.outer(u,u)-np.eye(16)
  elif name=='pointer':expected=b.pointer.target()
  if expected is not None:assert np.linalg.norm(actual-expected,2)<1e-11
  elif name=='uniform':assert np.linalg.norm(actual[:,0]-np.r_[np.ones(12)/math.sqrt(12),np.zeros(4)])<1e-11
  local.append({'template':name,'all_input_columns':1<<bits,'pulses':len(word),'complete_gram_residual':float(gram),
                'exact_operator_error_bound':str(b.a.template_bound(word,64)),
                'ideal_target_operator_residual':float(np.linalg.norm(actual-expected,2))if expected is not None else None})
 # Native pointer checks all matrix units under C and E.
 fam,run=b.pointer_blocks();pointer=[]
 for law in ('C','E'):
  for outcome in (0,1):
   instrument=sum(np.kron(k.conj(),k)for meta,ks in fam[law]if meta['bit']==outcome for k in ks)
   P=np.diag([int(outcome==0),int(outcome==1)])
   error=float(np.max(abs(instrument-np.kron(P,P))))
   assert error<1e-11;pointer.append({'law':law,'bit':outcome,'full_matrix_unit_error':error})
 # Clock phase generated from actual source pulse and read by actual pointer.
 clock=[]
 for size in (2,4):
  for quadrature in ('X','Y'):
   word=b.H(0)+[('Z',(0,),b.a.angle(0,0,-F(1,size*16)))]*size
   if quadrature=='Y':word+=[('Z',(0,),b.a.angle(F(1,4)))]
   word+=b.H(0);r=b.complete_run(word,1,np.array([1,0]));state=r['success']
   for law in ('C','E'):
    probs=np.zeros(2)
    for meta,ks in fam[law]:probs[meta['bit']]+=sum(np.linalg.norm(k@state)**2 for k in ks)
    contrast=float(probs[0]-probs[1]);target=math.cos(1/8)if quadrature=='X'else math.sin(1/8)
    assert abs(contrast-target)<1e-11
    clock.append({'size':size,'quadrature':quadrature,'law':law,'probabilities':probs.tolist(),'contrast':contrast,
                  'clock_phase_radians_exact':'1/8','physical_seconds':False})
 formal=[b.family_certificate(n,n)for n in (16,64,256,1024,4096)]
 assert all(F(r['delta_preparation_exact'])<=F(r['h_exact'])**2 and F(r['delta_W_exact'])<=F(r['h_exact'])**2 for r in formal)
 assert all(F(formal[i+1]['uncapped_TV_upper_exact'])<F(formal[i]['uncapped_TV_upper_exact'])for i in range(len(formal)-1))
 # Invalid members, preparation domains, missing failures and pulse mutations.
 controls=[]
 for fn in (lambda:b.member(3),lambda:b.family_certificate(4,4),lambda:b.admission(16,internal='port0'),
            lambda:b.admission(16,preparation_kind='installed_low_mode_vector'),lambda:b.admission(16,intermediate=True)):
  controls.append(require_failure(fn))
 bad=copy.deepcopy(run);bad['failure_complement_counts']=bad['failure_complement_counts'][:-1]
 controls.append(require_failure(lambda:b.require_complete(bad)))
 event=b.a.certified_pulse('Z',(0,),b.a.angle(F(1,4)),64);event['controls_exact'][0]='1'
 controls.append(require_failure(lambda:b.a.verify_pulse(event)))
  # External-reference selector is neither reset nor replaced by a prescribed preparation.
 external=next(b.native_program(2,1,external_selector=True));selector=8+2*b.member(2)[0]
 assert selector not in external['blank_wires']
 assert selector in next(b.native_program(2,1))['blank_wires']
 controls.append(True)
 # Native reset samples retain every fine letter and timeout.
 forced=b.execute_reset(b.circuit.decode(0),1,np.random.default_rng(1),True)
 assert not forced['success']
 trials=[b.sample_trial(2,2,'plus',law,19+i)for i,law in enumerate(('C','E'))]
 assert all(r['outcome']=='address'for r in trials)
 (HERE/'NATIVE_RECORDS.json').write_text(json.dumps(trials,indent=2)+'\n')
 # Exact compression identities fix signed current with the inherited P0.
 bloch=b.continuum.rational_bloch_data();P=np.asarray(bloch['P'],float)
 J=b.current_table(2)[:256].reshape(16,16,3)[:12,:12].reshape(144,3)
 signed=[float(np.linalg.norm(P@np.diag(J[:,axis])@P+np.asarray(bloch['H'][axis]),2))for axis in range(3)]
 assert max(signed)<1e-11
 assert b.check_pins()==pins
 result={'schema':'oph.chain2.refinement_apparatus.result.v1','model_id':b.MODEL,
  'checks':{'immutable_source_closure':True,'native_two_mode_arbitrary_reference_isometry':True,
   'full_M6_small_blocks_executed':True,'complete_native_C_E_pointer':True,'nonuniform_density_and_phase_controls':True,'dephased_selector_uniform_density_control':True,
   'source_all_clean_columns_sizes2_and4':True,'terminal_current_sign':True,'native_reset_success_and_timeout_executed':True,
   'missing_failure_and_domain_rejected':True,'general_family_error_tends_to_zero':True},
  'input_file_count':pins,'source_execution':natives,'development':rows,'local_blocks':local,'pointer':pointer,'clock':clock,
  'formal_bounds_only':formal,'current_compression_residuals':signed,'negative_controls':controls,
  'full_native_trial_record_file':'NATIVE_RECORDS.json','native_trials':len(trials),
  'forced_reset_timeout_attempts':sum(x.outcome.attempts for x in forced['blocks']),
  'execution_scope':'M6 preparation/clock/pointer and actual native reset trajectories; W certified block compression. The compressed-law samples have explicit native source TV enclosure. No larger graph execution.',
  'theorem':'fixed two-mode scalar preparation; all144 internal modes and padded addresses retained; terminal-only detector',
  'official_rows_changed':[]}
 from generate import build
 program=build();(HERE/'PROGRAM.json').write_text(json.dumps(program,separators=(',',':'))+'\n')
 result['program_sha256']=hashlib.sha256((HERE/'PROGRAM.json').read_bytes()).hexdigest()
 (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'checks':result['checks'],'formal_TV_bounds':[float(F(x['complete_record_TV_upper_exact']))for x in formal]}),flush=True)
if __name__=='__main__':main()
