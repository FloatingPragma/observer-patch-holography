"""Independent finite-state/crossblock review; no producer imports or experiments."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,math,sys,time
import numpy as np
B=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[2]/'work/2026-11-16_codex_common_refinement_family'
k=int(sys.argv[2]) if len(sys.argv)>2 else 2
start=time.monotonic();v=json.loads(gzip.decompress((B/f'k{k}/EVIDENCE.json.gz').read_bytes()));r=v['base_execution'];Q=256*2**k
checks=0;counts={'source_catalogs':0,'quantum_stages':0,'reference_blocks':0,'ramsey_records':0,'classical_events':0,'marker_histories':0}
def need(ok,label):
 global checks
 checks+=1
 if not ok:raise AssertionError(label)
def array(x):return np.array(x['real'])+1j*np.array(x['imag']) if isinstance(x,dict) else np.array(x)
def near(a,b,label,tol=3e-10):need(np.max(np.abs(np.asarray(a)-np.asarray(b)))<tol,label)
def event(q,e,c,law='native'):
 need(c in (0,1),'valid coin');q=list(q);a,b=e;d=q[a]-q[b]
 j=(d-(2*c-1)*(abs(d)%2))//2 if law=='native' else d if law=='swap' else c*d
 q[a]-=j;q[b]+=j;return q,j
TRI=((0,1),(1,2),(0,2));ENERGY=((0,1),(0,1),(0,2),(0,2),(1,2))
def source(edges,law='native',blank=False):
 raw=[]
 for e in edges:
  for c in (0,1):
   block=np.zeros((3,3,3))
   for p in range(3):
    out,j=event([int(i==p) for i in range(3)],e,c,law) if not blank else ([int(i==p) for i in range(3)],0)
    block[j+1,out.index(1),p]=1
   raw.extend(block/math.sqrt(2*len(edges)))
 total=len(raw);out=[sum(raw)/math.sqrt(total)];partial=raw[0].copy()
 for k in range(1,total):out.append((partial-k*raw[k])/math.sqrt(k*(k+1)));partial+=raw[k]
 return np.array(out)
def contrast(edges,law='native'):
 return np.array([(b+s*a)/2 for a,b in zip(source(edges,law),source(edges,law,True)) for s in (1,-1)])
def route(slot,m):
 p,k=divmod(slot,m);return [e for e in TRI if p in e]*(m-k)+[e for e in TRI if (p+1)%3 in e]*k
R=np.eye(33)
for p in range(3):
 for t,e in enumerate(ENERGY):
  for c in (0,1):
   out,j=event([int(i==p) for i in range(3)],e,c);dest=11*out.index(1)+2*t+(c^bool(j))
   R[:,11*p+2*t+c]=0;R[dest,11*p+2*t+c]=1
eps=2.**-k;a=np.r_[np.full(10,eps/math.sqrt(10)),math.sqrt(1-eps**2)];A=np.kron(np.eye(3),a[:,None])
P=A.T@R@A;W=R@(2*A@A.T-np.eye(33));near(P,np.eye(3)-eps**2*(np.eye(3)-np.array([[.6,.2,.2],[.2,.7,.1],[.2,.1,.7]])),'native transfer')
for key,x in [('A',A),('R',R),('W',W)]:near(x,array(r['clock'][key]),'clock '+key)
waves=[];live=A.astype(complex)
for _ in range(Q):waves.append(live);live=W@live
qpe=np.fft.fft(np.array(waves),axis=0)/Q
qpe_gram=np.einsum('kxa,kxb->kab',qpe.conj(),qpe)
positive_map=np.einsum('kxa,kyb->xayb',qpe[1:Q//2],qpe[1:Q//2].conj(),optimize=True)
u=a.copy();u[0]-=1;H=np.eye(11)-2*np.outer(u,u)/(u@u)
decoder=np.array([np.kron(np.eye(3),H[k:k+1]) for k in range(11)])
ops={}
for key,encoded in r['catalog'].items():
 if key=='qpe':value=qpe
 elif key=='decode':value=decoder
 elif key=='cooling':value=source(ENERGY)
 elif key=='cooling_contrast':value=contrast(ENERGY)
 elif key.endswith(':plane'):value=contrast(TRI,key.split(':')[0])
 else:
  m,law,_,slot=key.split(':');value=source(route(int(slot),int(m)),law)
 near(value,array(encoded),'source catalog '+key);near(np.einsum('kba,kbc->ac',value.conj(),value),np.eye(3 if key!='decode' else 33),'complete '+key)
 ops[key]=value;counts['source_catalogs']+=1
def channel(rho,kraus,accept):
 if kraus is qpe:
  need(list(accept)==list(range(1,Q//2)),'fixed complete positive counter selection')
  selected=np.einsum('xayb,ab->xy',positive_map,rho);mass=np.trace(selected).real
  return selected/mass,mass,np.einsum('kab,ba->k',qpe_gram,rho).real
 out=np.einsum('kxa,ab,kyb->kxy',kraus,rho,kraus.conj());prob=np.trace(out,axis1=1,axis2=2).real
 selected=out[accept].sum(axis=0);mass=np.trace(selected).real;return selected/mass,mass,prob
rho=np.diag([0.,1.,0.]);prefix=0.
for prep_index,row in enumerate(r['preparation_records']):
 near(rho,array(row['input_density']),'actual preparation input');near(prefix,row['log10_reach_before'],'preparation reach')
 accepted=row['accepted'];accepted=accepted if isinstance(accepted,list) else [accepted]
 rho,p,probs=channel(rho,ops[row['catalog']],accepted);near(probs,row['probabilities'],'complete preparation probabilities');prefix+=math.log10(p)
 if prep_index==8:near(rho,array(r['cooled_position_density']),'cooled final')
 if row['catalog']=='qpe':near(rho,array(r['positive_counter_density']),'finite counter state');positive=rho.copy()
need(len(r['preparation_records'])==11,'eleven complete preparation stages')
near(rho,array(r['decoded_position_density']),'decoded actual mixed input');decoded=rho.copy()
need(np.linalg.matrix_rank(decoded,tol=1e-9)>1,'input not replaced by ideal ray')
f=np.array([0.,1.,-1.])/math.sqrt(2);g=np.array([-2.,1.,1.])/math.sqrt(6);basis=np.stack([f,g],axis=1)
lambdas=np.diag(basis.T@P@basis);J=(R@A@basis-A@basis@np.diag(lambdas))@np.diag(1/np.sqrt(1-lambdas**2))
plus=(A@basis-1j*J)/math.sqrt(2);minus=(A@basis+1j*J)/math.sqrt(2)
for chart in (plus,minus):
 near(chart.conj().T@chart,np.eye(2),'signed chart isometry');near(decoder[0]@chart,basis/math.sqrt(2),'both frequency signs survive decode')
if k==2:near(float(np.trace(decoder[0]@positive@decoder[0].conj().T).real),.5002919412895217,'independent prior finite decode')
# Endpoint audit validates the actual classical source and the decoder output;
# the producer's separate audit reconstructs the current-only inverse itself.
for prep in r['charge_records']:
 q=prep['initial_audit_only'].copy()
 for section in ('body','probe','route'):
  for aa,bb,c,j in prep[section]:
   q,jj=event(q,(aa,bb),c,prep['law'] if section=='body' else 'native');need(jj==j,'classical current');counts['classical_events']+=1
  if section=='probe':need(q[:24]==prep['reading']['after_probe'],'actual selected sector')
 if prep['outcome']=='target_basis_unit':need(q[prep['decoded_target']]==1,'actual routed unit')
work_edges=((0,1),(0,2),(1,2),(3,4),(3,5),(4,5),(0,3),(1,4),(2,5))
work_edges=[oriented for e in work_edges for oriented in (e,e[::-1])]
need(len(r['transport']['word'])==21,'actual 21-event transfer word')
for p0 in range(6):
 p,c=p0,0
 for step,primitive_id in enumerate(r['transport']['word']):
  out,j=event([int(x==p) for x in range(6)],work_edges[primitive_id],c);y=out.index(1);cc=c^bool(j)
  need([primitive_id,p,c,j,y,cc]==r['transport']['traces'][p0][step],'native transport record');p,c=y,cc
 need([p,c]==[(p0+3)%6,0],'whole-state swap and restored work coin')
T=np.kron(np.array([[0,1],[1,0]]),np.eye(33));near(T@np.kron(np.eye(2),W),np.kron(np.eye(2),W)@T,'full clock transport intertwiner')
perm=np.eye(3)[[1,2,0]];spin_input=perm@decoded@perm.T;results=[]
for case,corrected in zip(r['cases'],v['cases']):
 m,t,d,law=(case[x] for x in ('M','turns','direction','law'));rho=np.kron(np.ones((2,2))/2,spin_input);logp=0.;positions={'rot':2*m,'ref':2*m}
 for row in case['records']:
  near(logp,row['log10_reach_before'],'rotation reach')
  if row['kind']=='marker':
   success=0.
   for history in row['histories']:
    need(len(history['coins'])==len(history['currents']) and 1<=len(history['coins'])<=8,'marker stopped tape')
    need(history['probability']==2**-len(history['coins']),'marker branch mass')
    need(all(j==0 for j in history['currents'][:-1]),'marker stops at first move')
    marker=[int(x==positions[row['arm']]) for x in range(3*m)]
    for coin,current in zip(history['coins'],history['currents']):
     marker,jj=event(marker,row['edge'],coin);need(jj==current,'common marker current')
    need(marker.index(1)==history['decoded'],'common current-decoded marker')
    if history['outcome']=='hop':need(history['currents'][-1]!=0,'actual marker move');success+=history['probability']
    else:need(len(history['coins'])==8 and history['currents'][-1]==0,'marker timeout')
    counts['marker_histories']+=1
   step=row['step'];offset=step if row['arm']=='rot' else min(step%m,m-step%m)
   positions[row['arm']]=(2*m+d*offset)%(3*m)
   need(all(h['decoded']==positions[row['arm']] for h in row['histories'] if h['outcome']=='hop'),'actual common marker schedule')
   logp+=math.log10(success);continue
  near(rho,array(row['input_density']),'actual coherent six-dimensional input')
  if row['kind']=='target':left=ops[f"{m}:{law}:target:{row['rot']}"];right=ops[f"{m}:{law}:target:{row['ref']}"]
  else:left=right=ops[law+':plane']
  joint=np.zeros((len(left),6,6),complex);joint[:,:3,:3]=left;joint[:,3:,3:]=right
  rho,p,prob=channel(rho,joint,[row['accepted']]);near(prob,row['probabilities'],'complete common-flag probabilities');logp+=math.log10(p);counts['quantum_stages']+=1;counts['reference_blocks']+=2
 near(rho,array(case['spin_final_density']),'same finite spin output')
 change=np.kron(np.eye(2),perm.T);energy=change@rho@change.T;near(energy,array(case['reencode_input']),'physical energy/spin permutation')
 blocks=energy.reshape(2,3,2,3);out=np.zeros((2,33,2,33),complex)
 for i in range(2):
  for j in range(2):out[i,:,j,:]=np.einsum('xayb,ab->xy',positive_map,blocks[i,:,j,:]);counts['reference_blocks']+=1
 final=out.reshape(66,66);p=np.trace(final).real;final/=p;logp+=math.log10(p)
 near(final,array(case['final_density']),'full reencoded outer crossblocks')
 marginal=blocks[0,:,0,:]+blocks[1,:,1,:];_,_,probs=channel(marginal,qpe,list(range(1,Q//2)))
 near(probs,case['all_reencode_bin_probabilities'],'all reencoding bins')
 bb=final.reshape(2,33,2,33);path=[np.trace((bb[0,:,0,:]+bb[1,:,1,:]+s*(bb[0,:,1,:]+bb[1,:,0,:]))/2).real for s in (1,-1)]
 near(path,case['path_probabilities'],'both path outcomes');near(2*path[0]-1,case['visibility'],'signed interference')
 selected=0 if t==2 or law=='swap' else 1;sign=1 if selected==0 else -1
 arc=(bb[0,:,0,:]+bb[1,:,1,:]+sign*(bb[0,:,1,:]+bb[1,:,0,:]))/(2*path[selected]);near(arc,array(case['clock_density']),'actual selected final arc state')
 power=np.eye(33);unit=np.linalg.matrix_power(W,2**k);moments=[]
 for record in case['ramsey_records']:
  z=np.trace(arc@power);moments.append(z);near([(1+z.real)/2,(1+z.imag)/2],[record['x_plus'],record['y_plus']],'raw Ramsey record law');power=unit@power;counts['ramsey_records']+=1
 z=np.array(moments);real=float(z[1].real);eta_bound=.002
 interval=[math.acos(max(-1,min(1,(real+eta_bound)/(1-eta_bound)))),math.acos(max(-1,min(1,(real-eta_bound)/(1-eta_bound))))]
 expected=math.acos(float(f@P@f))/eps
 near(interval,[case['stable_reader']['lower'],case['stable_reader']['upper']],'stable source-blind interval')
 need(interval[0]<=expected<=interval[1] and interval[1]-interval[0]<.01,'fixed native family interval gate')
 need(max(abs(z[n].real-math.cos(n*expected)) for n in range(2,17))<=.004+1e-10,'heldout real moments')
 fidelity=float(np.vdot(plus[:,0],arc@plus[:,0]).real);near(fidelity,case['positive_fidelity'],'finite positive mode fidelity')
 eta=float(1-fidelity-np.vdot(minus[:,0],arc@minus[:,0]).real)
 before=np.zeros((66,66),complex);before[:33,:33]=arc;after=T@before@T.T;near(after[33:,33:],arc,'whole actual returned state transported')
 near(logp,case['log10_conditional_success'],'all quantum reach');near(logp+prefix,case['log10_total_with_preparation'],'conditional prefix cost')
 factor=F(1)
 for n,count in corrected['selector_counts'].items():
  n=int(n);b=(n-1).bit_length();factor*=(1-F(2**b-n,2**b)**8)**count
 need(factor==F(corrected['selector_success_factor']),'exact all-selector survival');near(logp+prefix+math.log10(float(factor)),corrected['log10_success_conditional_on_classical_prefix'],'v2 complete conditional cost')
 pair=v['actual_returned_pair_records'][len(results)];pairz=z[1]*z[1]
 near([pairz.real,pairz.imag],pair['pair_ramsey'],'actual product Ramsey record')
 pair_eta=1-fidelity**2-float(1-fidelity-eta)**2
 near(pair_eta,pair['audited_pair_off_band_weight'],'complete pair off-band, including cross-sign zero modes')
 need(-1e-10<=eta<=.002+1e-10 and -1e-10<=pair_eta<=.01+1e-10,'frozen single/pair certificates')
 pair_interval=[math.acos(max(-1,min(1,(pairz.real+.01)/.99))),math.acos(max(-1,min(1,(pairz.real-.01)/.99)))]
 near(pair_interval,[pair['stable_pair_reader']['lower'],pair['stable_pair_reader']['upper']],'bounded pair reader')
 need(pair_interval[0]<=2*expected<=pair_interval[1] and pair_interval[1]-pair_interval[0]<.03,'pair interval gate')
 results.append({'M':m,'turns':t,'direction':d,'law':law,'visibility':float(2*path[0]-1),'positive_fidelity':fidelity,'opposite_sign_weight':float(1-fidelity-eta),'off_band_weight':eta,'pair_off_band_weight':pair_eta,'stable_interval':interval,'pair_interval':pair_interval,'source_frequency_comparison':expected,'producer_four_pole_diagnostic':case['reader']})
need(len(results)==16,'full case grid')
print(json.dumps({'k':k,'verdict':'PASS','checks':checks,'counts':counts,'cases':results,'input_evidence_sha256':hashlib.sha256((B/f'k{k}/EVIDENCE.json.gz').read_bytes()).hexdigest(),'same_group_only':True,'producer_imports':False,'seconds':time.monotonic()-start,'finite_shot_precision_claim':False},indent=2))
