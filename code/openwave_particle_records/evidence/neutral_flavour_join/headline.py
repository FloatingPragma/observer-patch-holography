#!/usr/bin/env python3
"""Independent retained-record replay; no source experiment is regenerated."""
from __future__ import annotations
import gzip, hashlib, itertools, json, math, time
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'reverse-engineering-reality').is_dir())
ORDERS=tuple(itertools.permutations(range(3)))
ROUTE=((0,1),(1,5),(2,5),(0,2));INITIAL=(1,2,5)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def arr(x):return np.array(x['real'])+1j*np.array(x['imag']) if isinstance(x,dict) else np.array(x)
def move(pos,seam,c,hold=False):
 x=[int(i in pos) for i in range(12)];a,b=seam;d=x[a]-x[b]
 j=0 if hold else (d-(2*c-1)*(abs(d)%2))//2;x[a]-=j;x[b]+=j
 return tuple(i for i,v in enumerate(x) if v),j

def tag_map(pos,new,erase=False):
 if erase or pos==new:return list(range(6))
 origin,=set(pos)-set(new);destination,=set(new)-set(pos);result=[]
 for order in ORDERS:
  tags=dict(zip(pos,order,strict=True));tags[destination]=tags.pop(origin)
  result.append(ORDERS.index(tuple(tags[p] for p in new)))
 return result

def compose(a,b):return [a[b[i]] for i in range(6)]
def pmatrix(p):
 M=np.zeros((8,8));M[6,6]=M[7,7]=1
 for a,b in enumerate(p):M[b,a]=1
 return M

def path(m):
 pos=INITIAL;p=list(range(6))
 for k in range(m):
  seam=ROUTE[k%4];a,b=seam;c=0 if a in pos else 1
  new,j=move(pos,seam,c);assert j;p=compose(tag_map(pos,new),p);pos=new
 return pos,p

def work_replay(records,seams):
 x=slots=None;owner=['source',0];expected_owner=None
 for r in records:
  kind=r[0]
  if kind=='reference':x=[v[:] for v in r[1]];slots=[[None if z is None else z[:] for z in row] for row in r[2]];continue
  if kind=='owner':
   assert r[1]==expected_owner,(r,expected_owner);owner=r[1][:];expected_owner=None;continue
  i=r[1]
  if kind in ('native','emit','capacity_hold'):
   _,i,e,c,j=r;a,b=seams[e];d=x[i][a]-x[i][b]
   if kind=='capacity_hold':assert abs(d)==2 and any(z is not None for z in slots[i]) and j==0
   else:
    assert j==(d-(2*c-1)*(abs(d)%2))//2;x[i][a]-=j;x[i][b]+=j
    if kind=='emit':
     assert abs(d)==2 and slots[i]==[None,None];slots[i]=[[e,c,j],[e,c,j]]
     if owner==['source',0] and i==0:expected_owner=['slot',0,0]
  elif kind=='swap':
   _,i,lane,j,other=r;slots[i][lane],slots[j][other]=slots[j][other],slots[i][lane]
   a,b=['slot',i,lane],['slot',j,other]
   expected_owner=b if owner==a else a if owner==b else owner[:]
  elif kind=='absorb':
   _,i,e,c,j=r;a,b=seams[e]
   assert slots[i]==[[e,c,-j],[e,c,-j]] and x[i][a]==x[i][b]==0
   x[i][a]-=j;x[i][b]+=j;slots[i]=[None,None]
   if owner[0]=='slot' and owner[1]==i:expected_owner=['detector',i]
  elif kind.startswith('detector_'):pass
  else:raise AssertionError(kind)
 assert expected_owner is None
 return owner,sum(sum(v*v for v in row) for row in x)+sum(z is not None for row in slots for z in row),sum(map(sum,x))

def main():
 started=time.perf_counter();receipt=json.loads((HERE/'RECEIPT.json').read_text())
 assert sha(HERE/'RECORDS.json.gz')==receipt['records_sha256']
 assert sha(HERE/'PREREG.txt')==receipt['prereg_sha256']
 for line in (HERE/'PREREG.txt').read_text().splitlines():
  h,p=line.split('  ',1);assert sha(HERE/p)==h
 for p,h in json.loads((HERE/'SOURCE_PINS.json').read_text()).items():assert sha(ROOT/p)==h,p
 data=json.loads(gzip.decompress((HERE/'RECORDS.json.gz').read_bytes()))
 col=data['columns'];pos=[tuple(p) for p in col['positions']];seams=[tuple(e) for e in col['seams']]
 seen={}
 for i,o,e,c,j,ni,no in col['rows']:
  new,current=move(pos[i],seams[e],c);assert new==pos[ni] and current==j
  assert tag_map(pos[i],new)[o]==no
  key=e,c,j;target=ni,no;bucket=seen.setdefault(key,set());assert target not in bucket;bucket.add(target)
 assert len(col['rows'])==79200
 frame=data['frame'];F=arr(frame['frame']);A=arr(frame['prepared']);B=arr(frame['backward']);effects=arr(frame['effects'])
 assert np.max(abs(F.T@F-np.eye(8)))<1e-12
 assert np.linalg.norm(sum(effects)-np.eye(8))<1e-9
 assert np.linalg.norm(arr(frame['forward'])-F,2)<2e-10 and np.linalg.norm(B-F.T,2)<2e-10
 # Independently rebuild the bit effects from the retained complete source matrix.
 local=arr(frame['pointer_matrix']);counts=frame['pointer_certificate']['exact_instrument']['native_effect_numerators_over_two']['C']
 bit_effects=[]
 for row in counts:
  effect=np.zeros((2,2),complex)
  for z,w in enumerate(row):
   block=local[np.ix_([z,z+4],[2,6])];effect+=w/2*block.conj().T@block
  bit_effects.append(effect)
 for label,effect in enumerate(effects):
  rebuilt=np.array([[1]],complex)
  for bit in reversed(range(3)):rebuilt=np.kron(rebuilt,bit_effects[(label>>bit)&1])
  assert np.max(abs(rebuilt-effect))<1e-14
 chart=data['laws']['chart'];S=np.eye(8)[:,chart+[6,7]]
 U=S.T@pmatrix(path(4)[1])@S;assert np.array_equal(U,arr(data['laws']['cycle_matrix']))
 maxappearance=0.;minreturn=1.;laws=data['laws']['rows']
 for row in laws:
  a=row['input'];k=row['horizon'];assert sum(Q(b['weight']) for b in row['branches'])==1
  for branch in row['branches']:
   m=branch['moves'];P=S.T@pmatrix(path(m)[1])@S;psi=B@P@A[:,a]
   probs=[float(np.vdot(psi,E@psi).real) for E in effects]
   assert np.max(abs(np.array(probs)-branch['probabilities']))<1e-13
   assert abs(sum(probs)-1)<1e-9
   ideal=F.T@P@F[:,a];assert np.max(abs(np.array(probs)-abs(ideal)**2))<1e-9
   s=Q(255,256);expected=s**(4*k) if branch['flag']=='success' else (1-s)*s**m
   assert Q(branch['weight'])==expected
   if branch['flag']=='success' and k==1:maxappearance=max(maxappearance,max(p for b,p in enumerate(probs) if b!=a))
   if branch['flag']=='success' and k==3:minreturn=min(minreturn,probs[a])
 assert maxappearance>1e-4 and minreturn>1-1e-9
 committed={}
 for a in range(6):
  rho=np.outer(F[:,a],F[:,a])
  for k in range(1,7):
   rho=U@rho@U.T;rho=sum(np.outer(F[:,b],F[:,b])*float(np.vdot(F[:,b],rho@F[:,b]).real) for b in range(8))
   committed[a,k]=[float(np.vdot(F[:,b],rho@F[:,b]).real) for b in range(8)]
 for row in data['laws']['committed']:assert np.max(abs(np.array(committed[row['input'],row['horizon']])-row['probabilities']))<1e-13
 loss=max(1-committed[a,3][a] for a in range(6));assert loss>1e-5
 internal_count=0;work_count=0
 for flight in data['flights']:
  e,sign,c,k,mode=flight['case'];posi=INITIAL;mapping=list(range(6));route=tuple(reversed(ROUTE)) if mode=='reverse' else ROUTE
  for cycle,stage,attempt,before,seam,coin,j,after,owner in flight['internal']:
   assert tuple(before)==posi and tuple(seam)==route[stage]
   new,current=move(posi,tuple(seam),coin,mode=='hold');assert j==current and tuple(after)==new
   assert owner[0]=='slot'
   mapping=compose(tag_map(posi,new,mode=='erase'),mapping);posi=new;internal_count+=1
  assert list(posi)==flight['position'] and mapping==flight['tag_map']
  assert Q(flight['internal_fine_tape_probability'])==Q(1,2**len(flight['internal']))
  if flight['emitted']:
   P=S.T@pmatrix(mapping)@S
   for a,retained in enumerate(flight['pointer_probabilities_conditioned_ready_auxiliaries']):
    psi=B@P@A[:,a];probs=[float(np.vdot(psi,E@psi).real) for E in effects]
    assert np.max(abs(np.array(probs)-retained))<1e-13
  else:assert flight['pointer_probabilities_conditioned_ready_auxiliaries']==[]
  owner,budget,charge=work_replay(flight['work_records'],seams)
  assert owner==flight['owner'] and budget==flight['source_budget'] and charge==0
  if mode=='no_contact':assert owner==['source',0]
  else:assert owner==['detector',k if k else 0]
  work_count+=len(flight['work_records'])
 reset_checks=0
 for item in data['reset']['records']:
  x=item['initial'][:]
  for seam,coin,j in item['execution']['records']:
   a,b=seam;d=x[a]-x[b];assert j==(d-(2*coin-1)*(abs(d)%2))//2;x[a]-=j;x[b]+=j;reset_checks+=1
  assert x==item['execution']['output']
 ap=data['complete_apparatus'];p=Q(ap['exact_reset_success'])
 assert sum(map(Q,ap['first_reset_failure_masses']))+Q(ap['all_six_reset_success'])==1
 assert Q(ap['complete_conditional_operator_envelope'])<Q(1,10**6)
 primary={tuple(r['case'][:4]):r for r in data['flights'][:840]};differences={'erase':0.,'reverse':0.}
 for row in data['flights'][840:1320]:
  e,sign,c,k,mode=row['case'];base=primary[e,sign,c,k]
  if mode=='erase':assert row['internal']==base['internal']
  if mode in differences and row['timeout'] is None and base['timeout'] is None:
   differences[mode]=max(differences[mode],float(np.max(abs(np.array(row['pointer_probabilities_conditioned_ready_auxiliaries'])-base['pointer_probabilities_conditioned_ready_auxiliaries']))))
 assert differences['erase']>1e-4 and differences['reverse']>1e-6
 cur=data['current'];W=arr(cur['full_frame']);J=arr(cur['full_current']);trans=W.conj().T@J@W
 assert np.max(abs(trans-arr(cur['full_transformed_current'])))<1e-13
 assert np.max(abs(W.conj().T@W-np.eye(45)))<1e-12
 assert data['laws']['squared_principal_phases_pi2_units']==['0','4/9','4/9']
 report={'status':'PASS_FINITE_NATIVE_THREE_TAG_JOIN','tag_columns':79200,'flights':len(data['flights']),
   'internal_attempts_replayed':internal_count,'work_records_replayed':work_count,'reset_attempts_replayed':reset_checks,
   'max_appearance':maxappearance,'minimum_return_at3':minreturn,'committed_return_loss':loss,
   'full_PMNS_criterion_satisfied':False,'seconds':time.perf_counter()-started}
 print(json.dumps(report,sort_keys=True,indent=2))
if __name__=='__main__':main()
