"""Native two-Fourier-mode preparation and terminal C/E detector.

The declared pulse program uses exact rational controls. Development execution
uses complete M6 local blocks and a certified reversible-source compression.
The compression's full native error is explicit; it is not a pulse simulation.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib,json,math,sys
import numpy as np
HERE=Path(__file__).resolve().parent
INPUT=HERE/'inputs'
sys.path[:0]=[str(INPUT/'openwave-tuning/next/CHAIN2/local12_limit_audit'),str(INPUT/'openwave-tuning/common_dynamics/continuum'),str(INPUT/'openwave-tuning/common_dynamics/particles')]
import refinement_adapter as a
import continuum
import mass_coin
from model.native_interface import pointer
from model.native_clock import configuration
from model.observer import native
from model.observer.source_discrimination import circuit

MODEL='native-local12-two-fourier-mode-apparatus-v1'

def check_pins():
 pins=json.loads((HERE/'INPUTS.json').read_text())['files_sha256']
 for p,v in pins.items():
  if hashlib.sha256((INPUT/p).read_bytes()).hexdigest()!=v:raise ArithmeticError('frozen source bytes changed: '+p)
 return len(pins)

def member(size):
 if type(size)is not int or size<2 or size&(size-1):raise ValueError('dyadic cubic size N=2^q, q>=1 required')
 q=size.bit_length()-1
 return q,8+3*q

def remap(word,wires):return [(k,tuple(wires[j]for j in t),v)for k,t,v in word]
def H(wire):return remap(a.primitive_templates()['hadamard'],[wire])
def phase(wire,coefficient):
 # exp(i*pi*c*b) with its scalar phase retained.
 return [('I',(wire,),a.angle(F(coefficient,2))),('Z',(wire,),a.angle(-F(coefficient,2)))]
def controlled_phase(left,right,coefficient):
 c=F(coefficient,4)
 return remap(a.pauli('II',a.angle(c))+a.pauli('ZI',a.angle(-c))+a.pauli('IZ',a.angle(-c))+a.pauli('ZZ',a.angle(c)),[left,right])
def swap(left,right):return a.cnot(left,right)+a.cnot(right,left)+a.cnot(left,right)

def fourier(wires):
 """Chronological +DFT on little-endian register; no Fourier service."""
 word=[]
 for top in reversed(range(len(wires))):
  word+=H(wires[top])
  for low in reversed(range(top)):
   word+=controlled_phase(wires[low],wires[top],F(1,1<<(top-low)))
 for low in range(len(wires)//2):word+=swap(wires[low],wires[-1-low])
 return word

def preparation_template(size):
 """Isometry with selector x-bit0 supplied; other data bits start zero."""
 q,d=member(size);word=[]
 word+=remap(a.uniform_template(),[0,1,2,3])
 word+=remap(a.uniform_template(),[4,5,6,7])
 for bit in range(8,8+2*q):word+=H(bit)
 word+=fourier(list(range(8+2*q,d)))
 return word

def selector_template(size,setting):
 q,d=member(size);wire=8+2*q
 if setting not in ('plus','minus','plus_i','minus_i','zero','one'):raise ValueError('predeclared selector setting required')
 if setting=='zero':return []
 if setting=='one':return remap(a.primitive_templates()['flip'],[wire])
 return H(wire)+phase(wire,{'plus':F(0),'minus':F(1),'plus_i':F(1,2),'minus_i':F(-1,2)}[setting])

def detector_template():
 """Auxiliary low bits start code2, high bit is measured source bit."""
 return a.cnot(2,0)+a.cnot(0,1)

def emit(word,bits):return [a.certified_pulse(k,t,v,bits)for k,t,v in word]

@lru_cache(None)
def response(kind,angle,precision):
 event=a.certified_pulse(kind,tuple(range(2 if kind=='J'else 1)),angle,precision)
 a.verify_pulse(event)
 branches=a.execute_full_response(event);d=4 if kind=='J'else 2
 return (branches[0][:d],*(K[d:] for K in branches[1:]))

def apply_local(state,matrix,targets,bits):
 # Same target order as the proper-code loading: first target is high local bit.
 cols=state.shape[1]
 axes=[bits-1-j for j in targets]+[j for j in range(bits)if bits-1-j not in targets]+[bits]
 inverse=np.argsort(axes)
 shaped=state.reshape([2]*bits+[cols]).transpose(axes)
 out=matrix@shaped.reshape(1<<len(targets),-1)
 return out.reshape(shaped.shape).transpose(inverse).reshape(state.shape)

def apply_complement(state,matrix,targets,bits):
 cols=state.shape[1]
 axes=[bits-1-j for j in targets]+[j for j in range(bits)if bits-1-j not in targets]+[bits]
 arranged=state.reshape([2]*bits+[cols]).transpose(axes).reshape(1<<len(targets),-1)
 return (matrix@arranged).reshape(-1,cols)

def complete_run(word,bits,initial,precision=64):
 """Execute M6 success and every first-export complement, no renormalization.

 Each stored Gram row is computed from a literal prefix Kraus map. The full
 failure poststate is reconstructible by failure_branch with the same tape.
 """
 state=np.asarray(initial,complex).copy()
 if state.ndim==1:state=state[:,None]
 if state.shape[0]!=1<<bits:raise ValueError('complete padded register required')
 grams=[]
 for k,t,v in word:
  branches=response(k,tuple(v),precision)
  gram=np.zeros((state.shape[1],state.shape[1]),complex)
  for K in branches[1:]:
   failed=apply_complement(state,K,t,bits);gram+=failed.conj().T@failed
  grams.append(gram)
  state=apply_local(state,branches[0],t,bits)
 return {'success':state,'failure_grams':grams,'word':word,'initial':initial,'bits':bits,'precision':precision,'failure_complement_counts':[6-(4 if k=='J'else 2)for k,t,v in word]}

def failure_branch(run,index,complement):
 if not 0<=index<len(run['word']):raise ValueError('first failure index required')
 state=np.asarray(run['initial'],complex).copy()
 if state.ndim==1:state=state[:,None]
 for j,(k,t,v)in enumerate(run['word']):
  branches=response(k,tuple(v),run['precision'])
  if j==index:
   if not 1<=complement<len(branches):raise ValueError('export complement required')
   return apply_complement(state,branches[complement],t,run['bits'])
  state=apply_local(state,branches[0],t,run['bits'])

def prepare(size,setting='plus',precision=64):
 q,d=member(size);initial=np.zeros(1<<d,complex);initial[0]=1
 return complete_run(selector_template(size,setting)+preparation_template(size),d,initial,precision)

def expected_isometry(size):
 """Independent validation target; never a preparation command."""
 q,d=member(size);result=np.zeros((1<<d,2),complex)
 for x in range(size):
  for y in range(size):
   for z in range(size):
    c=(x*size+y)*size+z
    for port in range(12):
     for event in range(12):result[256*c+16*port+event]=np.array([1,np.exp(2j*np.pi*x/size)])/np.sqrt(144*size**3)
 return result

@lru_cache(None)
def pointer_blocks(precision=64):
 initial=np.eye(8,dtype=complex)[:,[2,6]]
 run=complete_run(detector_template(),3,initial,precision)
 families={}
 for law in ('C','E'):
  groups=defaultdict(lambda:defaultdict(lambda:np.zeros((2,2),complex)))
  metadata={}
  for code in range(4):
   source=circuit.decode(code);K=run['success'][[code,code+4]]
   for coin in (0,1):
    output,(_,current)=native.native(source,pointer.SEAM,coin)
    label=current if law=='C'else source[0]-source[1]
    key=(coin,label);bit=int(2*current+2*coin-1<0)
    groups[key][output]+=K/np.sqrt(2)
    metadata[key]={'coin':coin,'label':label,'current':current,'bit':bit}
  families[law]=[(metadata[k],list(g.values()))for k,g in sorted(groups.items())]
 return families,run

@lru_cache(None)
def detector_response(law='C',precision=64):
 """Complete sequential bit law incl individual first-export outcomes.

 For probabilities the pointer's effects are diagonal. This is checked from
 the actual M6 word; all off-diagonal residuals enter a numerical diagnostic.
 Fine records and failure poststates remain in pointer_blocks/failure_branch.
 """
 if law not in ('C','E'):raise ValueError('native C or E law required')
 families,run=pointer_blocks(precision);rows=[];offdiag=0.
 for meta,ks in families[law]:
  effect=sum(K.conj().T@K for K in ks)
  offdiag=max(offdiag,float(abs(effect[0,1])))
  rows.append(({'kind':'native',**meta},np.real(np.diag(effect))))
 for i,g in enumerate(run['failure_grams']):
  offdiag=max(offdiag,float(abs(g[0,1])))
  rows.append(({'kind':'export_failure','pulse':i},np.real(np.diag(g))))
 return rows,offdiag

def complete_address_law(probabilities,law='C',precision=64):
 p=np.asarray(probabilities,float);d=(len(p)-1).bit_length()
 if len(p)!=1<<d or min(p)<-1e-15:raise ValueError('complete nonnegative padded input law required')
 rows,offdiag=detector_response(law,precision)
 bit_response=np.zeros((2,2));failure=np.zeros(2)
 for meta,v in rows:
  if meta['kind']=='native':bit_response[meta['bit']]+=v
  else:failure+=v
 alive=p.copy();failures=[]
 for bit in range(d):
  indices=np.arange(len(p));one=(indices>>bit)&1
  failures.append(float(np.sum(alive*failure[one])))
  paired=alive.reshape(-1,2,1<<bit)
  alive=np.einsum('ab,nbc->nac',bit_response,paired).reshape(-1)
 return {'address_probabilities':alive,'first_detector_export_failures':failures,
         'offdiagonal_effect_diagnostic':offdiag,'offdiagonal_probability_TV_diagnostic':len(rows)*d*offdiag/2,
         'evaluation_scope':'diagonal-effect numerical compression; full coherent instrument is pointer_blocks, with certified detector error',
         'fine_record_factorization':rows}

@lru_cache(None)
def reset_base():
 receipt=json.loads((a.PROVIDER/'COMMON_WALK_RECEIPT.json').read_text())['preparation_reset']['certificate']
 ready=configuration.certificate(4,32,full=True,mode='reset',target=pointer.ready_state())
 return max(F(receipt['failure_upper']),F(ready['failure_upper'])),receipt,ready

def reset_plan(size,tolerance):
 q,d=member(size);memory=a.layout(size)['memory_bits']
 # All source memory (including selector for prescribed inputs), clock, and
 # one fresh detector auxiliary per data bit plus clock read.
 blank=memory+1;ready=d+1;count=blank+ready
 base,zero,orient=reset_base();r=1
 while count*base**r>F(tolerance):r+=1
 cap=max(zero['native_attempts_max'],orient['native_attempts_max'])
 return {'blank_carriers':blank,'ready_carriers':ready,'registers':count,
         'blocks_per_register':r,'uniform_q_exact':str(base),'timeout_bound_exact':str(count*base**r),
         'native_attempts_max':count*r*cap,'base_attempts_max':cap,
         'fine_records':'every attempt seam, coin, current and final configuration; first register timeout halts'}

def execute_reset(target,repeats,rng,force_failure=False):
 """Actual finite native controller; no desired-success output commands."""
 initial=native.states()[0];blocks=[]
 for block in range(repeats):
  def coins():
   while True:yield 0 if force_failure else int(rng.integers(2))
  execution=configuration.execute(initial,coins(),4,32,full=True,mode='reset',target=target)
  blocks.append(execution)
  if execution.outcome.success:return {'success':True,'blocks':blocks}
  initial=execution.output
 return {'success':False,'blocks':blocks}

def family_certificate(size,steps,tau0=F(1)):
 q,d=member(size);h=F(1,size);K=F(7)
 if h*K>F(1,2)or type(steps)is not int or steps<0 or F(tau0)<=0:raise ValueError('N>=16, nonnegative steps and positive common tau0 required')
 resource=a.resource_upper_bound(size)
 # A uniform per-gate envelope avoids constructing the large graph.
 gates=resource['reversible_gates_upper'];precision=max(64,(resource['native_pulses_upper']*size*size).bit_length()+40)
 prep=selector_template(size,'plus_i')+preparation_template(size)
 while True:
  primitive=max(a.template_bound(w,precision)for w in a.primitive_templates().values())
  source=a.template_bound(a.reflection_template(),precision)+gates*primitive
  prep_error=a.template_bound(prep,precision)
  pointer_error=a.template_bound(detector_template(),precision)
  clock_word=H(0)+[('Z',(0,),a.angle(0,0,-h/16))]*steps+[('Z',(0,),a.angle(F(1,4)))]+H(0)
  clock_error=a.template_bound(clock_word,precision)
  instrument=(d+1)*pointer_error+clock_error
  if max(source,prep_error,instrument)<=h*h:break
  precision*=2
 resets=reset_plan(size,h*h)
 continuum_error=F(5,2)*h*K+F(77,32)*steps*h*h*K*K
 total=continuum_error+prep_error+steps*source+instrument+F(resets['timeout_bound_exact'])
 return {'size':size,'h_exact':str(h),'tau0_exact':str(tau0),'steps':steps,'time_exact':str(steps*h*tau0),
  'bandwidth_bound_exact':'7','support':'{(0,0,0),(2*pi,0,0)} on unit torus','precision_bits':precision,
  'delta_preparation_exact':str(prep_error),'delta_W_exact':str(source),'n_delta_W_exact':str(steps*source),
  'clock_and_terminal_detector_error_exact':str(instrument),'continuum_error_exact':str(continuum_error),
  'complete_record_TV_upper_exact':str(min(F(1),total)),'uncapped_TV_upper_exact':str(total),
  'cell_integral_extra_exact':str(F(3,4)*h*K),'reset':resets,
  'native_pulses_upper':steps*resource['native_pulses_upper']+len(prep)+(d+1)*len(detector_template())+len(clock_word),
  'memory_bits':a.layout(size)['memory_bits']+1+2*(d+1),'development_executed':False,
  'arbitrary_reference':True,'physical_seconds':False}

def current_table(size):
 """Arrival-port displacement signs, fixed by source oriented links."""
 _,d=member(size);out=np.zeros((1<<d,3),int)
 for cell in range(size**3):
  for axis,links in enumerate(a.spatial.links()):
   for left,right in links:
    out[256*cell+16*left+11,axis]=-1
    out[256*cell+16*right+11,axis]=1
 return out

def compressed_walk(state,size,steps):
 """Certified native-template compression, not literal complete-W simulation."""
 _,d=member(size);value=np.asarray(state,complex).reshape(1<<d,-1).copy()
 program=a.reversible_program(size)
 uniform=np.r_[np.ones(12)/math.sqrt(12),np.zeros(4)]
 reflection=2*np.outer(uniform,uniform)-np.eye(16)
 for _ in range(steps):
  value=apply_local(value,reflection,(3,2,1,0),d)
  out=np.empty_like(value);out[program['routing_permutation']]=value;value=out
 return value

def continuum_state(size,steps,setting):
 """Comparison only, never fed into the native preparation or sampler."""
 q,d=member(size);h=1/size;e=np.ones(144)/12
 coef={'plus':(1,1),'minus':(1,-1),'plus_i':(1,1j),'minus_i':(1,-1j),'zero':(math.sqrt(2),0),'one':(0,math.sqrt(2))}[setting]
 Hx=np.array(continuum.rational_bloch_data()['H'][0],float)
 # Full144 exponent retains stationary modes. No acoustic projection.
 from scipy.linalg import expm
 mode=expm(steps*h*continuum.effective_generator((2*np.pi,0,0)))@e
 result=np.zeros(1<<d,complex)
 for cell in range(size**3):
  x=cell//size**2
  block=(coef[0]*e+coef[1]*np.exp(2j*np.pi*x/size)*mode)/math.sqrt(2*size**3)
  for port in range(12):result[256*cell+16*port:256*cell+16*port+12]=block[12*port:12*port+12]
 return result

def native_program(size,steps,setting='plus',quadrature='X',precision=64,external_selector=False):
 """Lazy chronological executable program, including all terminal operations."""
 q,d=member(size)
 if type(steps)is not int or steps<0 or quadrature not in ('X','Y'):raise ValueError('finite steps and X/Y clock setting required')
 memory=a.layout(size)['memory_bits'];clock=memory;h=F(1,size)
 plan=reset_plan(size,h*h)
 if type(external_selector)is not bool:raise ValueError('explicit selector mode required')
 selector=8+2*q
 blanks=[wire for wire in range(memory+1)if not(external_selector and wire==selector)]
 yield {'operation':'native_configuration_reset_batch','plan':plan,'blank_wires':blanks,'source_target_code':0,'pointer_target_code':2,'ready_auxiliary_pairs':[[clock+1+2*j,clock+2+2*j]for j in range(d+1)]}
 selector_word=[]if external_selector else selector_template(size,setting)
 for event in emit(selector_word+preparation_template(size)+H(clock),precision):yield {'operation':'M6',**event}
 for _ in range(steps):
  for event in a.pulse_stream(size,precision):yield {'operation':'M6',**event}
  yield {'operation':'M6',**a.clock_tick_pulse(size,h,precision,clock)}
 analysis=([]if quadrature=='X'else[('Z',(clock,),a.angle(F(1,4)))])+H(clock)
 for event in emit(analysis,precision):yield {'operation':'M6',**event}
 for j,bit in enumerate([*range(d),clock]):
  wires=[clock+1+2*j,clock+2+2*j,bit]
  for event in emit(remap(detector_template(),wires),precision):yield {'operation':'M6',**event}
  yield {'operation':'native_C_or_E','seam':[0,1],'auxiliary_wires':wires[:2],
         'record':['coin','native_label','current','bit'],'all_full_carrier_outputs':True}

def sample_trial(size=2,steps=2,setting='plus',law='C',seed=1,precision=64):
 """Sample generated M6 prep/pointer and actual native resets.

 W is evaluated by certified block compression. The receipt labels its TV
 enclosure; sampled records are not represented as exact full-W pulse draws.
 """
 from dataclasses import asdict
 if law not in ('C','E'):raise ValueError('C/E required')
 rng=np.random.default_rng(seed);q,d=member(size);plan=reset_plan(size,F(1,size*size));resets=[]
 for index in range(plan['registers']):
  target=circuit.decode(0 if index<plan['blank_carriers']else 2)
  result=execute_reset(target,plan['blocks_per_register'],rng)
  resets.append({'register':index,'target_code':0 if index<plan['blank_carriers']else 2,
                 'success':result['success'],'blocks':[asdict(x)for x in result['blocks']]})
  if not result['success']:return {'outcome':'reset_timeout','resets':resets,'first_failed_register':index}
 prepared=prepare(size,setting,precision)
 weights=[float(np.trace(g).real)for g in prepared['failure_grams']]
 success=float(np.vdot(prepared['success'],prepared['success']).real)
 pick=float(rng.random())*(success+sum(weights));acc=0.
 for index,p in enumerate(weights):
  acc+=p
  if pick<acc:return {'outcome':'preparation_export_failure','pulse':index,'resets':resets}
 state=compressed_walk(prepared['success'],size,steps)
 state/=np.linalg.norm(state) # conditional trajectory state, success probability retained above
 families,ptr=pointer_blocks(precision);records=[];address=0
 for bit in range(d):
  alternatives=[]
  for meta,ks in families[law]:
   for environment,K in enumerate(ks):
    out=apply_local(state,K,(bit,),d)
    alternatives.append((float(np.vdot(out,out).real),meta,environment,out))
  # All pointer first-export effects; halt preserves the prefix failure state.
  for index,g in enumerate(ptr['failure_grams']):
   transformed=apply_local(state,g,(bit,),d)
   probability=float(np.vdot(state,transformed).real)
   alternatives.append((probability,{'kind':'export_failure','pulse':index},None,None))
  total=sum(max(0,x[0])for x in alternatives);draw=float(rng.random())*total;running=0
  for weight,meta,environment,out in alternatives:
   running+=max(0,weight)
   if draw<=running:break
  if meta.get('kind')=='export_failure':return {'outcome':'detector_export_failure','read_bit':bit,'pulse':meta['pulse'],'resets':resets,'records':records}
  records.append({'read_bit':bit,**meta,'retained_environment_branch':environment})
  address|=meta['bit']<<bit;state=out/math.sqrt(weight)
 clock_state=complete_run(H(0)+[('Z',(0,),a.angle(0,0,-F(1,size*16)))]*steps+H(0),1,np.array([1,0]),precision)
 clock_probs=[]
 for meta,ks in families[law]:
  for K in ks:clock_probs.append((float(np.linalg.norm(K@clock_state['success'])**2),meta))
 clock_failure=float(sum(np.trace(g).real for g in clock_state['failure_grams']))
 for g in ptr['failure_grams']:clock_failure+=float((clock_state['success'].conj().T@g@clock_state['success'])[0,0].real)
 draw=float(rng.random())*(sum(p for p,m in clock_probs)+clock_failure);running=0.;clock_record=None
 for weight,meta in clock_probs:
  running+=weight
  if draw<=running:clock_record=meta;break
 return {'outcome':'address'if clock_record is not None else'clock_export_failure','address':address,
         'current':current_table(size)[address].tolist(),'clock_record':clock_record,'resets':resets,'records':records,
         'sampling':'actual reset and M6 apparatus; certified compressed W centre with separately stated native TV error'}

def admission(size,preparation_kind='native_two_mode_isometry',internal='uniform_e',intermediate=False):
 member(size)
 if preparation_kind!='native_two_mode_isometry'or internal!='uniform_e'or intermediate:
  raise ValueError('native two-mode scalar isometry and terminal-only observation required')
 return {'preparation':preparation_kind,'internal':internal,'P0':True,'bandwidth_K':7,'terminal_only':True}

def require_complete(run):
 expected=[6-(4 if k=='J'else 2)for k,t,v in run['word']]
 if run.get('failure_complement_counts')!=expected or len(run['failure_grams'])!=len(run['word']):
  raise ValueError('every first-export complement, including zero-weight branches, is required')
 return True

def development_budget(size,steps,precision=64):
 """Same emitted precision and full counted word as development programme."""
 q,d=member(size);h=F(1,size);counts=Counter(k for k,t in a.reversible_program(size)['gates'])
 source=a.template_bound(a.reflection_template(),precision)+sum(counts[k]*a.template_bound(w,precision)for k,w in a.primitive_templates().items()if k in counts)
 prep=a.template_bound(selector_template(size,'plus_i')+preparation_template(size),precision)
 clock=H(0)+[('Z',(0,),a.angle(0,0,-h/16))]*steps+[('Z',(0,),a.angle(F(1,4)))]+H(0)
 detector=(d+1)*a.template_bound(detector_template(),precision)+a.template_bound(clock,precision)
 resets=reset_plan(size,h*h);reset=F(resets['timeout_bound_exact'])
 return {'precision_bits':precision,'delta_preparation_exact':str(prep),'delta_W_exact':str(source),
         'n_delta_W_exact':str(steps*source),'clock_and_detector_exact':str(detector),
         'reset':resets,'native_vs_ideal_complete_law_TV_exact':str(min(F(1),prep+steps*source+detector+reset)),
         'source_compression_TV_exact':str(min(F(1),steps*source)),
         'scope':'full logical/record law certificate; binary64 development arithmetic residuals are diagnostic'}
