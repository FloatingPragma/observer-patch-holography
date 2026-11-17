#!/usr/bin/env python3
"""Frozen native three-tag source and persistent-payload observation test."""
from __future__ import annotations
import argparse, gzip, hashlib, importlib.util, itertools, json, math, random, sys, time
from collections import Counter
from dataclasses import asdict
from fractions import Fraction as Q
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'reverse-engineering-reality').is_dir())
RER=ROOT/'reverse-engineering-reality'
ROUTE=((0,1),(1,5),(2,5),(0,2));REVERSE=tuple(reversed(ROUTE))
INITIAL=(1,2,5);ORDERS=tuple(itertools.permutations(range(3)))
CAP=8;HORIZONS=tuple(range(7))
NATIVE_CALLS=0

def canonical(x):return json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load_module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def imports():
 global np,native,compiler,pointer,configuration,producer,old,work,internal,action
 import numpy as np
 sys.path.insert(0,str(RER))
 from model.observer import native
 from model.observer.source_scalar import compiler
 from model.native_interface import pointer
 from model.native_clock import configuration
 from model.common_source import native_repair_mass as producer
 from model.common_source import calibrated_matter as internal
 from model.common_source.chiral_candidate import action
 old=load_module('neutral_join_old',ROOT/'openwave/row_10_neutrino_oscillations/work/2026-11-14_codex_neutral_even_records/verify.py')
 wp=ROOT/'openwave/shared/rows08_10_reversible_work';sys.path.insert(0,str(wp))
 work=load_module('neutral_join_work',wp/'verify.py')

def pincheck():
 pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
 for name,h in pins.items():assert sha(ROOT/name)==h,('source changed',name)
 for line in (HERE/'PREREG.txt').read_text().splitlines():
  h,name=line.split('  ',1);assert sha(HERE/name)==h,('freeze changed',name)
 return sha(HERE/'PREREG.txt')

def endpoint(positions,seam,coin,mode='native'):
 global NATIVE_CALLS
 state=tuple(int(i in positions) for i in range(12))
 out,(_,j)=native.native(state,seam,coin);NATIVE_CALLS+=1
 if mode=='hold':out,j=state,0
 new=tuple(i for i,v in enumerate(out) if v)
 assert len(new)==3 and all(v in (0,1) for v in out)
 a,b=seam;assert all(out[i]==state[i]+j*(int(i==b)-int(i==a)) for i in range(12))
 return new,j

def transport_order(positions,new,order):
 oldset,newset=set(positions),set(new)
 if oldset==newset:return order
 origin,=oldset-newset;destination,=newset-oldset
 bypos=dict(zip(positions,order,strict=True));bypos[destination]=bypos.pop(origin)
 return tuple(bypos[p] for p in new)

def map_step(positions,new,mode='native'):
 return [i if mode=='erase' else ORDERS.index(transport_order(positions,new,o)) for i,o in enumerate(ORDERS)]

def compose(left,right):return [left[right[i]] for i in range(len(right))]
def matrix(mapping):
 out=np.zeros((8,8));out[6,6]=out[7,7]=1
 for i,j in enumerate(mapping):out[j,i]=1
 return out

def advance_path(moves,route=ROUTE):
 positions=INITIAL;mapping=list(range(6))
 for n in range(moves):
  seam=route[n%4];a,b=seam
  assert (a in positions)!=(b in positions)
  new=tuple(sorted((set(positions)-{a if a in positions else b})|{b if a in positions else a}))
  mapping=compose(map_step(positions,new),mapping);positions=new
 return positions,mapping

def chart():
 _,c=advance_path(4);ordered=[]
 for seed in (ORDERS.index((0,1,2)),ORDERS.index((0,2,1))):
  orbit=[seed,c[seed],c[c[seed]]];assert c[orbit[-1]]==seed
  ordered.extend(orbit)
 assert len(set(ordered))==6
 S=np.eye(8)[:,ordered+[6,7]]
 return ordered,S,c

def full_columns():
 positions=tuple(itertools.combinations(range(12),3));idx={p:i for i,p in enumerate(positions)}
 rows=[];checks=0
 for ei,seam in enumerate(work.SEAMS):
  branches={}
  for i,pos in enumerate(positions):
   for coin in (0,1):
    new,j=endpoint(pos,seam,coin)
    image=map_step(pos,new)
    for oi,oj in enumerate(image):
     key=coin,j;target=idx[new],oj
     seen=branches.setdefault(key,set());assert target not in seen;seen.add(target)
     rows.append([i,oi,ei,coin,j,idx[new],oj]);checks+=1
  assert sum(map(len,branches.values()))==2*220*6
 assert checks==79200
 return {'positions':[list(p) for p in positions],'orders':[list(o) for o in ORDERS],
         'seams':[list(e) for e in work.SEAMS],'rows':rows,'columns_checked':checks,
         'complete_raw_slots':[[c,j] for c in (0,1) for j in(-1,0,1)]}

def pointer_effects():
 local=pointer.local_matrix();effect=[]
 for bit,counts in enumerate(pointer.native_effects('C')):
  e=np.zeros((2,2),complex)
  for z,w in enumerate(counts):
   B=local[np.ix_([z,z+4],[pointer.READY_CODE,pointer.READY_CODE+4])]
   e+=w/2*(B.conj().T@B)
  effect.append(e)
 full=[]
 for label in range(8):
  e=np.array([[1]],complex)
  for bit in reversed(range(3)):e=np.kron(e,effect[(label>>bit)&1])
  full.append(e)
 assert np.linalg.norm(sum(full)-np.eye(8))<1e-9
 return np.array(full)

def frame_readout():
 raw,branches=old.source_records();d=producer.architecture();Y=np.array(d['Y'].evalf(30),float)
 vals,F=np.linalg.eigh(Y)
 # Fix numerical eigenvector signs without fitting any overlap.
 for j in range(3):
  pivot=int(np.argmax(abs(F[:,j])))
  if F[pivot,j]<0:F[:,j]*=-1
 frame=np.eye(8);frame[:3,:3]=frame[3:6,3:6]=F
 word=compiler.compile_real_orthogonal(frame)
 forward=compiler.execute(word,3);backward=compiler.execute(compiler.inverse(word),3)
 errors=[float(np.linalg.norm(forward-frame,2)),float(np.linalg.norm(backward-frame.T,2))]
 assert max(errors)<2e-10
 cert=pointer.certificate();effects=pointer_effects()
 basis_words={a:[compiler.pulse('X',math.pi/2,[bit]) for bit in range(3) if (a>>bit)&1] for a in range(8)}
 basis_full=[compiler.execute(basis_words[a],3) for a in range(8)]
 basis=np.column_stack([m[:,0] for m in basis_full])
 basis_errors=[]
 for a,m in enumerate(basis_full):
  target=np.zeros((8,8),complex)
  for j in range(8):target[j^a,j]=(1j)**a.bit_count()
  basis_errors.append(float(np.linalg.norm(m-target,2)))
 assert max(basis_errors)<1e-10
 prepared=forward@basis
 # Complete blank-prepare channel: measured bit followed by conditional native X.
 x=compiler.execute([compiler.pulse('X',math.pi/2,[0])],1)
 assert np.linalg.norm(abs(x)-np.array([[0,1],[1,0]]))<1e-10
 local=pointer.local_matrix();blank_error=0.
 for i in range(2):
  for j in range(2):
   out=np.zeros((2,2),complex);rho=np.zeros((2,2),complex);rho[i,j]=1
   for bit,counts in enumerate(pointer.native_effects('C')):
    for z,w in enumerate(counts):
     B=local[np.ix_([z,z+4],[pointer.READY_CODE,pointer.READY_CODE+4])]
     result=B@rho@B.conj().T
     if bit:result=x@result@x.conj().T
     out+=w/2*result
   target=np.array([[int(i==j),0],[0,0]],complex)
   blank_error=max(blank_error,float(np.max(abs(out-target))))
 assert blank_error<1e-9
 return {'raw_records':raw,'raw_branches':branches,'Y_exact':old.serial(d['Y']),
   'charged_eigenvalues':vals.tolist(),'word':word,'inverse_word':compiler.inverse(word),
   'word_errors':errors,'word_pulses':len(word),'pointer_certificate':{k:str(v) if isinstance(v,Q) else v for k,v in cert.items()},
   'frame':frame,'forward':forward,'backward':backward,'prepared':prepared,'basis_words':basis_words,'basis_matrices_first_columns':basis,'basis_full_matrices':basis_full,'basis_word_errors':basis_errors,'blank_superoperator_error':blank_error,'effects':effects,'pointer_matrix':pointer.local_matrix()}

def reset_audit():
 target=pointer.ready_state();bound=configuration.budget(full=True,mode='reset',target=target)
 probs=configuration.success_probabilities(full=True,mode='reset',target=target)
 assert 1-min(probs)<=bound['failure_upper']
 records=[];checks=0
 for i,initial in enumerate(native.states()):
  rng=random.Random(int.from_bytes(hashlib.sha256(f'neutral-three-tag/reset/{i}'.encode()).digest(),'big'))
  run=configuration.execute(initial,(rng.randrange(2) for _ in range(bound['native_attempts_max'])),full=True,mode='reset',target=target,law='C')
  st=list(initial)
  for seam,coin,j in run.records:
   a,b=seam;d=st[a]-st[b];expected=(d-(2*coin-1)*(abs(d)%2))//2
   assert j==expected;st[a]-=j;st[b]+=j;checks+=1
  assert tuple(st)==run.output
  records.append({'initial':list(initial),'execution':asdict(run)})
 return {'records':records,'literal_attempts':checks,'exact_success_probabilities':[str(x) for x in probs],
   'timeout_bound':str(bound['failure_upper']),'timeout_max_exact':str(1-min(probs)),
   'attempt_cap':bound['native_attempts_max'],'successful':sum(r['execution']['outcome']['success'] for r in records),
   'six_bit_blank_plus_terminal_error_bound':str(6*(bound['failure_upper']+pointer.certificate()['half_diamond_upper']))}

def flight(e,sign,coin,k,mode='native',forced=False):
 records=[]
 class OwnedWorld(work.World):
  def __init__(self,L,log):super().__init__(L,log);self.owner=['source',0]
  def swap(self,i,lane,j,other):
   super().swap(i,lane,j,other)
   a,b=['slot',i,lane],['slot',j,other]
   if self.owner==a:self.owner=b
   elif self.owner==b:self.owner=a
   self.log.append(['owner',self.owner[:]])
  def detect(self,i,enabled=True):
   kind=super().detect(i,enabled)
   if kind=='absorb' and self.owner[0]=='slot' and self.owner[1]==i:
    self.owner=['detector',i];self.log.append(['owner',self.owner[:]])
   return kind
 L=2*k if k else 2;w=OwnedWorld(L,records);a,b=work.SEAMS[e]
 if mode!='no_contact':w.x[0][a]=sign;w.x[0][b]=-sign
 w.reference();emitted=w.charge_event(0,e,coin)=='emit'
 if emitted:w.owner=['slot',0,0];records.append(['owner',w.owner[:]])
 positions=INITIAL;mapping=list(range(6));internal=[];timeout=None;steps=0
 seed_mode='native' if mode=='erase' else mode
 seed=f'neutral-three-tag/{e}/{sign}/{coin}/{k}/{seed_mode}/{int(forced)}'
 rng=random.Random(int.from_bytes(hashlib.sha256(seed.encode()).digest(),'big'))
 route=REVERSE if mode=='reverse' else ROUTE
 for cycle in range(k):
  if emitted and timeout is None:
   for stage,seam in enumerate(route):
    moved=False
    for attempt in range(CAP):
     c=rng.randrange(2)
     if forced and cycle==stage==0:
      c=1 if seam[0] in positions else 0
     new,j=endpoint(positions,seam,c,'hold' if mode=='hold' else 'native')
     delta=map_step(positions,new,'erase' if mode=='erase' else 'native')
     mapping=compose(delta,mapping)
     internal.append([cycle,stage,attempt,list(positions),list(seam),c,j,list(new),w.owner[:]])
     positions=new
     if j:steps+=1;moved=True;break
    if not moved:timeout=[cycle,stage,CAP];break
  w.cycle(cycle)
 if k==0:w.detect(0)
 if emitted:assert w.owner==['detector',k if k else 0]
 else:assert w.owner==['source',0]
 assert not any(w.errors.values())
 if emitted and timeout is None:
  expected=advance_path(4*k,route)[1] if mode!='erase' else list(range(6))
  assert positions==INITIAL and mapping==expected
 return {'case':[e,sign,coin,k,mode],'seed':seed,'emitted':emitted,'owner':w.owner,
   'timeout':timeout,'moves':steps,'position':list(positions),'tag_map':mapping,
   'internal':internal,'internal_fine_tape_probability':str(Q(1,2**len(internal))),'work_records':records,'work_native_calls':w.native_calls,
   'record_counts':dict(w.counts),'source_budget':w.energy(),'source_charge':w.charge()}

def laws(fr):
 ordered,S,c=chart();U=S.T@matrix(c)@S
 F=fr['frame'];forward=fr['prepared'];back=fr['backward'];effects=fr['effects']
 s=Q(255,256);rows=[];maxnorm=0.;maxerror=0.;maxappearance=0.;minreturn=1.
 ideal_dephase=lambda rho:sum(np.outer(F[:,a],F[:,a])*float(np.vdot(F[:,a],rho@F[:,a]).real) for a in range(8))
 for k in HORIZONS:
  branches=[(m,(1-s)*s**m,'timeout') for m in range(4*k)]+[(4*k,s**(4*k),'success')]
  assert sum(p for _,p,_ in branches)==1
  for a in range(6):
   outputs=[]
   for m,weight,flag in branches:
    P=S.T@matrix(advance_path(m)[1])@S
    psi=back@P@forward[:,a];ideal=F.T@P@F[:,a]
    probs=[float(np.vdot(psi,E@psi).real) for E in effects]
    maxnorm=max(maxnorm,abs(sum(probs)-1));maxerror=max(maxerror,max(abs(np.array(probs)-abs(ideal)**2)))
    if flag=='success' and k==1:maxappearance=max(maxappearance,max(p for b,p in enumerate(probs) if b!=a))
    if flag=='success' and k==3:minreturn=min(minreturn,probs[a])
    outputs.append({'moves':m,'weight':str(weight),'flag':flag,'probabilities':probs})
   rows.append({'horizon':k,'input':a,'branches':outputs})
 assert maxnorm<1e-9 and maxerror<1e-9 and maxappearance>1e-4 and minreturn>1-1e-9
 committed=[]
 for a in range(6):
  rho=np.outer(F[:,a],F[:,a])
  for k in range(1,7):
   rho=ideal_dephase(U@rho@U.T)
   probs=[float(np.vdot(F[:,b],rho@F[:,b]).real) for b in range(8)]
   committed.append({'input':a,'horizon':k,'probabilities':probs})
 loss=max(1-r['probabilities'][r['input']] for r in committed if r['horizon']==3)
 assert loss>1e-5
 planted=[]
 for a in range(6):
  psi=back@forward[:,a];p=[float(np.vdot(psi,E@psi).real) for E in effects]
  assert p[a]>1-1e-9;planted.append(p)
 assert any(p[0]<.5 for p in planted) and any(sum(p[:3])<.5 for p in planted)
 return {'chart':ordered,'cycle_matrix':U,'rows':rows,'committed':committed,'planted':planted,
  'max_normalization_defect':maxnorm,'max_ideal_difference':maxerror,'max_appearance':maxappearance,
  'min_coherent_return_at3':minreturn,'max_committed_return_loss_at3':loss,
  'success_weights':{str(k):str(s**(4*k)) for k in HORIZONS},
  'principal_phases_pi_units':['0','2/3','-2/3'],
  'squared_principal_phases_pi2_units':['0','4/9','4/9'],
  'two_distinct_nonzero_mass_squared_gaps_under_principal_identification':False}

def current_join(fr):
 import sympy as sp
 ell=internal._LH_tensor()*sp.Matrix([0,1]);nu=next(i for i in range(15) if ell[i]!=0)
 D=internal._oriented_components('charged_lepton')[1];charged,singlet=next(iter(D.todok()))
 J=np.array(internal.weak_raising_current().evalf(30),complex)
 _,S,c=chart();C=(S.T@matrix(c)@S)[:3,:3]
 omega=np.exp(2j*np.pi/3);Z=np.array([[omega**(-i*j)/math.sqrt(3) for j in range(3)] for i in range(3)])
 assert np.max(abs(C@Z-Z@np.diag([1,omega,omega**2])))<1e-12
 F=fr['frame'][:3,:3];W=np.eye(45,dtype=complex)
 W[3*charged:3*charged+3,3*charged:3*charged+3]=F
 W[3*nu:3*nu+3,3*nu:3*nu+3]=Z
 transformed=W.conj().T@J@W;sign=complex(internal.weak_raising_current(families=False)[nu,charged])
 B=transformed[3*nu:3*nu+3,3*charged:3*charged+3]
 error=float(np.max(abs(B-sign*Z.conj().T@F)));assert error<1e-12
 mixing=F.T@Z;mod=abs(mixing)**2
 G=np.eye(3)[[2,0,1],:];covariance_error=float(np.max(abs((G@F).T@(G@Z)-mixing)))
 assert covariance_error<1e-12
 return {'neutral_index':nu,'charged_index':charged,'current_block_error':error,
  'common_frame_change_error':covariance_error,'full_current':J,'full_transformed_current':transformed,'full_frame':W,
  'cycle_frame':Z,'mixing':mixing,'mixing_moduli':mod,
  'diagnostic_sin2':[float(mod[0,1]/(1-mod[0,2])),float(mod[0,2]),float(mod[1,2]/(1-mod[0,2]))],
  'charged_Y_to_tag_chart_is_declared_apparatus':True,'physical_weak_interaction_measured':False}

def jsonable(x):
 if isinstance(x,np.ndarray):return {'real':x.real.tolist(),'imag':x.imag.tolist()} if np.iscomplexobj(x) else x.tolist()
 if isinstance(x,dict):return {k:jsonable(v) for k,v in x.items()}
 if isinstance(x,(tuple,list)):return [jsonable(v) for v in x]
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,Q):return str(x)
 return x

def run():
 start=time.perf_counter();freeze=pincheck();imports();records={}
 records['columns']=full_columns();fr=frame_readout();records['frame']=jsonable(fr)
 records['reset']=reset_audit();records['laws']=jsonable(laws(fr));records['current']=jsonable(current_join(fr))
 cases=[]
 for e,sign,c,k in itertools.product(range(30),(-1,1),(0,1),HORIZONS):cases.append(flight(e,sign,c,k))
 for e,sign,c,mode in itertools.product(range(30),(-1,1),(0,1),('hold','erase','reverse','no_contact')):cases.append(flight(e,sign,c,1,mode))
 forced=flight(0,1,0,1,forced=True);assert forced['timeout']==[0,0,CAP]
 cases.append(forced)
 ordered,S,_=chart();effects=fr['effects']
 for item in cases:
  item['pointer_probabilities_conditioned_ready_auxiliaries']=[]
  if item['emitted']:
   P=S.T@matrix(item['tag_map'])@S
   for a in range(6):
    psi=fr['backward']@P@fr['prepared'][:,a]
    probs=[float(np.vdot(psi,E@psi).real) for E in effects]
    assert abs(sum(probs)-1)<1e-9
    item['pointer_probabilities_conditioned_ready_auxiliaries'].append(probs)
 native_cases={(r['case'][0],r['case'][1],r['case'][2],r['case'][3]):r for r in cases[:840]}
 erased_difference=0.;reverse_difference=0.
 for item in cases[840:1320]:
  e,sign,c,k,mode=item['case'];base=native_cases[e,sign,c,k]
  if mode=='erase':
   assert item['internal']==base['internal']
   if item['timeout'] is None and base['timeout'] is None:
    erased_difference=max(erased_difference,float(np.max(abs(np.array(item['pointer_probabilities_conditioned_ready_auxiliaries'])-base['pointer_probabilities_conditioned_ready_auxiliaries']))))
  if mode=='reverse' and item['timeout'] is None and base['timeout'] is None:
   reverse_difference=max(reverse_difference,float(np.max(abs(np.array(item['pointer_probabilities_conditioned_ready_auxiliaries'])-base['pointer_probabilities_conditioned_ready_auxiliaries']))))
 assert erased_difference>1e-4 and reverse_difference>1e-6
 records['changed_law_reader_differences']={'erase':erased_difference,'reverse':reverse_difference}
 records['flights']=cases
 reset_p=sum(map(Q,records['reset']['exact_success_probabilities']))/78
 assert 0<reset_p<=1
 timeout=Q(records['reset']['timeout_bound']);pointer_delta=pointer.certificate()['half_diamond_upper']
 ef=Q(2,10**10);ex=Q(1,10**10);df=ef*(2+ef)/2;dx=ex*(2+ex)/2
 gain=(1+ef)**4*(1+ex)**12*(1+Q(1,10**11))**12
 complete_bound=gain*(6*(timeout+pointer_delta)+6*dx+2*df)
 assert complete_bound<Q(1,10**6)
 records['complete_apparatus']={'auxiliary_initial_law':'fresh independent uniform full78 at each reset','exact_reset_success':str(reset_p),'all_six_reset_success':str(reset_p**6),'first_reset_failure_masses':[str(reset_p**m*(1-reset_p)) for m in range(6)],'failure_records':'preparation reset stage0..2; after successful preparation, terminal read reset stage0..2 with preceding native bits and remaining quantum state retained by factorized instrument','conditional_pointer_laws_are_not_unconditional':True,'frame_word_operator_envelope':str(ef),'basis_and_conditional_X_operator_envelope':str(ex),'frame_word_halfdiamond_term':str(df),'X_word_halfdiamond_term':str(dx),'complete_conditional_operator_envelope':str(complete_bound),'conservative_complete_bound':'1/1000000','numeric_word_envelopes_are_not_authenticated_hardware_error_bounds':True}
 assert sum(reset_p**m*(1-reset_p) for m in range(6))+reset_p**6==1
 summary={'schema':'oph.neutral_three_tag_flight.v1','prereg_sha256':freeze,
  'native_endpoint_calls':NATIVE_CALLS,'charged_frame_literal_calls':144,'full_tag_columns':79200,'primary_flights':840,'control_flights':480,'forced_timeout_flights':1,
  'sampled_primary_timeouts':sum(x['timeout'] is not None for x in cases[:840]),
  'total_internal_attempt_records':sum(len(x['internal']) for x in cases),
  'rwork_native_calls':sum(x['work_native_calls'] for x in cases),
  'reset_native_attempts':records['reset']['literal_attempts'],
  'source_matrix_errors':fr['word_errors'],
  'max_appearance':records['laws']['max_appearance'],
  'min_coherent_return_at3':records['laws']['min_coherent_return_at3'],
  'committed_return_loss_at3':records['laws']['max_committed_return_loss_at3'],
  'weak_current_frame_error':records['current']['current_block_error'],'erase_pointer_difference':erased_difference,'reverse_pointer_difference':reverse_difference,'complete_apparatus_bound':'1/1000000',
  'all_gates_pass':True,'full_PMNS_criterion_satisfied':False,'row_complete':False,
  'runtime_seconds':time.perf_counter()-start}
 with gzip.GzipFile(filename=str(HERE/'RECORDS.json.gz'),mode='wb',mtime=0) as f:f.write(canonical(records).encode())
 summary['records_sha256']=sha(HERE/'RECORDS.json.gz')
 (HERE/'RECEIPT.json').write_text(canonical(summary));print(canonical(summary))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
 if not args.write:p.error('--write required for new frozen evidence')
 if (HERE/'RECEIPT.json').exists():p.error('sealed/previous result exists; preserve it before a new attempt')
 run()
