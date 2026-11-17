"""Independent density/crossblock, selector and finite-error verification."""
from pathlib import Path
from fractions import Fraction as F
import json,math,cmath
import numpy as np
B=Path(__file__).resolve().parents[2]/'work/2026-11-15_codex_subdivided_rotation'
EDGES=((0,1),(1,2),(0,2));I=np.eye(3);counts={'error_flags':0,'error_crossblocks':0,'coherence_stages':0,'selector_words':0}

def fine(edges,blank=False):
 out=[]
 for a,b in edges:
  for c in (0,1):
   block=np.zeros((3,3,3))
   for p in range(3):
    q=[int(x==p) for x in range(3)];d=q[a]-q[b];j=0 if blank else (d-(2*c-1)*(abs(d)%2))//2
    q[a]-=j;q[b]+=j;block[j+1,q.index(1),p]=1
   out.extend(block)
 return out

def eraser(maps):
 n=len(maps);count=n//3;prefix=np.zeros((3,3));ops=[]
 for k,a in enumerate(maps):
  if k:ops.append((prefix-k*a)/math.sqrt(count*k*(k+1)))
  prefix+=a
 return [prefix/math.sqrt(count*n)]+ops

def phase(edges,delta):
 a,b=eraser(fine(edges)),eraser(fine(edges,True))
 return [(y+s*cmath.exp(1j*delta)*x)/2 for x,y in zip(a,b) for s in (1,-1)]

def routes(slot,m):
 v,k=divmod(slot,m)
 return [e for e in EDGES if v in e]*(m-k)+[e for e in EDGES if (v+1)%3 in e]*k

def block(a,b):
 out=np.zeros((6,6),complex);out[:3,:3]=a;out[3:,3:]=b;return out

error=json.loads((B/'ERROR_RECEIPT.json').read_text());results=[]
for case in error['cases']:
 m,t,delta=case['M'],case['turns'],case['delta_radians'];prep=phase((EDGES[0],),delta)
 sigma=np.zeros((3,3),complex);initial=np.ones(3)/math.sqrt(3)
 for index,map0 in enumerate(fine(EDGES)):
  v=map0@initial/math.sqrt(6);ps=[float(np.vdot(op@v,op@v).real) for op in prep]
  stored=case['records'][index]
  assert stored['kind']=='preparation' and stored['warmup_flag']==index
  assert max(abs(np.array(ps)-stored['joint_flag_probabilities']))<1e-14
  selected=prep[1]@v;sigma+=np.outer(selected,selected.conj());counts['error_flags']+=len(ps)
 reach=float(np.trace(sigma).real);rho=np.kron(np.ones((2,2))/2,sigma/reach)
 targets=[eraser(fine(routes(slot,m))) for slot in range(3*m)];filters=phase(EDGES,delta)
 for record in case['records'][18:]:
  assert abs(record['reach']/reach-1)<3e-10
  if record['kind']=='marker':
   assert abs(record['timeout']/reach-1/256)<1e-13;reach*=255/256;continue
  raw=np.asarray(record['input_amplitudes']);z=raw[...,0]+1j*raw[...,1]
  ensemble=np.concatenate((z[0].T,z[1].T),axis=0)/math.sqrt(record['reach'])
  assert np.max(abs(ensemble@ensemble.conj().T-rho))<5e-10;counts['error_crossblocks']+=1
  if record['kind']=='target':
   step=record['step'];k=step%m;rot=(2*m+step)%(3*m);ref=(2*m+min(k,m-k))%(3*m)
   ops=[block(a,b) for a,b in zip(targets[rot],targets[ref])]
  else:ops=[block(a,a) for a in filters]
  outputs=[op@rho@op.conj().T for op in ops];ps=np.array([np.trace(x).real for x in outputs])
  assert max(abs(ps-np.array(record['flag_probabilities'])/record['reach']))<2e-10
  assert abs(sum(ps)-1)<1e-12;counts['error_flags']+=len(ps)
  accepted=record['accept_flag'];rho=outputs[accepted]/ps[accepted];reach*=ps[accepted]
 visibility=2*np.trace(rho[:3,3:]).real;f=np.array([1.,-1.,0.])/math.sqrt(2)
 fidelity=np.vdot(f,rho[:3,:3]@f).real/np.trace(rho[:3,:3]).real
 assert abs(visibility-case['visibility'])<1e-10 and abs(fidelity-case['rotated_contrast_fidelity'])<1e-10
 assert abs(reach/case['success']-1)<5e-10
 results.append({'M':m,'turns':t,'delta':delta,'visibility':visibility,'fidelity':fidelity})

coherent=json.loads((B/'COHERENCE_RECEIPT.json').read_text());cresults=[]
P=np.array([[F(int(i==j))-F(1,3) for j in range(3)] for i in range(3)],dtype=object)
for case in coherent['cases']:
 m,t=case['M'],case['turns'];f=[1,-1,0];prep=F(5,144)
 left=right=cross=np.array([[prep*F(f[i]*f[j],4) for j in range(3)] for i in range(3)],dtype=object)
 for record in case['records']:
  reach=sum(left[i,i]+right[i,i] for i in range(3));assert reach==F(record['reach'])
  if record['kind']=='marker':
   assert F(record['timeout'])==reach/256;left,right,cross=[x*F(255,256) for x in(left,right,cross)];continue
  assert [[str(x) for x in z.flat] for z in(left,right,cross)]==record['input_joint_density_blocks']
  if record['kind']=='committed_target':
   a=[np.array(x,int) for x in fine(routes(record['rotating_slot'],m))]
   b=[np.array(x,int) for x in fine(routes(record['reference_slot'],m))]
   left=sum((x@left@x.T for x in a),np.zeros((3,3),dtype=object))/F(48*m*m)
   right=sum((x@right@x.T for x in b),np.zeros((3,3),dtype=object))/F(48*m*m)
   cross=sum((x@cross@y.T for x,y in zip(a,b)),np.zeros((3,3),dtype=object))/F(48*m*m)
  else:left,right,cross=[P@x@P/F(48) for x in(left,right,cross)]
  after=sum(left[i,i]+right[i,i] for i in range(3));assert reach-after==F(record['first_rejection']);counts['coherence_stages']+=1
 total=sum(left[i,i]+right[i,i] for i in range(3));vis=2*sum(cross[i,i] for i in range(3))/total
 assert total==F(case['success']) and vis==F(case['exact_visibility'])
 cresults.append({'M':m,'turns':t,'visibility':float(vis)})

selectors=json.loads((B/'SELECTOR_RECEIPT.json').read_text())
assert len(selectors['selector_tables'])==72
for tab in selectors['selector_tables']:
 m,slot=tab['M'],tab['slot'];v,k=divmod(slot,m)
 assert len(tab['all_words'])==2*m
 for number,word in enumerate(tab['all_words']):
  bits=[(number>>j)&1 for j in reversed(range(m.bit_length()))]
  assert word['native_coin_word']==bits and word['vacuum_current_word']==[0]*len(bits)
  ticket,bit=divmod(number,2);star=v if ticket<m-k else (v+1)%3
  edge=[e for e in EDGES if star in e][bit]
  assert tuple(word['decoded_address'])==edge and word['star']==star and word['ticket']==ticket and word['edge_bit']==bit
  assert F(word['probability'])==F(1,2*m);counts['selector_words']+=1
confirmation=json.loads((B/'CONFIRM_RECEIPT.json').read_text())
for c,cost in zip(confirmation['cases'],selectors['endpoint_costs']):
 m,law=c['M'],c['law'];expected=F(2)
 for record in c['records']:
  per=(F(255,128) if law=='native' else F(1)) if record['kind']=='marker' else (m.bit_length()+1 if record['kind']=='target' else 1)
  expected+=F(record['reach'])*per
 assert expected==F(cost['expected_native_endpoint_opportunities_per_attempt'])
print(json.dumps({'verdict':'PASS','same_group_only':True,'counts':counts,'finite_error_cases':results,'committed_current_control':cresults,'three_way_uniform_selector':'Explicit supplied coherent apparatus, not counted as a finite exact fair-bit sampler.'},indent=2))
