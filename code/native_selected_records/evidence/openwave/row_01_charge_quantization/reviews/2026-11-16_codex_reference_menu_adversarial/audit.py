"""Independent finite-state/crossblock review; no producer imports or experiments."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,math,sys,time
import numpy as np
B=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[2]/'work/2026-11-16_codex_common_reference_menu'
k=int(sys.argv[2]) if len(sys.argv)>2 else 2
reference=sys.argv[3] if len(sys.argv)>3 else 'stationary'
evidence=B/f'k{k}_{reference}/EVIDENCE.json.gz'
start=time.monotonic();v=json.loads(gzip.decompress(evidence.read_bytes()));r=v['base_execution'];Q=256*2**k
correction=json.loads((B/'EXACT_PATH_CORRECTION.json').read_text())
amended=next(m for m in correction['members'] if m['k']==k and m['reference']==reference)
checks=0;counts={'source_catalogs':0,'quantum_stages':0,'reference_blocks':0,'ramsey_records':0,'classical_events':0,'marker_histories':0,'stationary_holds':0}
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
  if row['kind']=='marker_hold':
   need(reference=='stationary' and row['arm']=='ref' and row['position']==2*m and not row['native_event'] and row['probability']==1.,'literal reference hold without endpoint event')
   need(row['physical_position']==r['graphs'][0 if m==4 else 1]['rings'][1][2*m],'physical stationary marker')
   counts['stationary_holds']+=1
   continue
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
   need(reference!='stationary' or row['arm']=='rot','stationary reference cannot move')
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
 direct=float(2*np.trace(bb[0,:,1,:]).real/np.trace(final).real)
 normcross=float(np.trace(bb[0,:,1,:]).real/math.sqrt(np.trace(bb[0,:,0,:]).real*np.trace(bb[1,:,1,:]).real))
 near(direct,case['visibility_direct_crossblock'],'direct crossblock contrast',1e-12);near(normcross,case['normalized_crossblock'],'normalized crossblock',1e-12)
 if law=='native':need((direct<(-.9 if reference=='matched' else -1e-7) if t==1 else direct>(.9 if reference=='matched' else 1e-12)),'native signed fringe gate')
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
 factor=(1-F(3,8)**8)**11*(1-F(1,4)**8)**(3*m*t)
 need(factor==F(corrected['selector_success_factor']),'exact selector survival')
 analyzer=logp+prefix+math.log10(float(factor))
 near(analyzer,corrected['log10_analyzer_probability'],'analyzer reach from every accepted source event')
 # Reconstruct the exact competing-law path probability from literal transpositions.
 physical=list(path);exact_law=None
 if law=='swap':
  amplitudes=[F(1),F(1)]
  for step in range(1,3*m*t+1):
   slots=[(2*m+d*step)%(3*m),2*m if reference=='stationary' else (2*m+d*min(step%m,m-step%m))%(3*m)]
   for arm,slot in enumerate(slots):
    edges=route(slot,m);S=[[F(0) for _ in range(3)] for _ in range(3)]
    for a0,b0 in edges:
     for j in range(3):S[b0 if j==a0 else a0 if j==b0 else j][j]+=F(1,len(edges))
    E=[[F(int(i==j))-F(1,3) for j in range(3)] for i in range(3)]
    def mul(x,y):return [[sum(x[i][h]*y[h][j] for h in range(3)) for j in range(3)] for i in range(3)]
    Z=mul(E,mul(S,E));Z2=mul(Z,Z);lam=sum(Z2[i][i] for i in range(3))/2
    need(lam>0 and Z2==[[lam*x for x in row] for row in E],'source-derived positive optical square')
    amplitudes[arm]*=lam*lam
  aa,bb=amplitudes;den=2*(aa*aa+bb*bb);exact_law=[(aa+bb)**2/den,(aa-bb)**2/den]
  need(sum(exact_law)==1 and min(exact_law)>=0,'exact physical nonnegative complete swap law')
  if reference=='matched':need(exact_law==[F(1),F(0)],'exact dark-port law')
  physical=[float(x) for x in exact_law]
  scalar=np.einsum('xayb,ab->xy',positive_map,decoded);scalar/=np.trace(scalar).real
  near(arc,scalar,'actual finite swap output equals common CP map, not ideal ray',1e-10)
 need(min(physical)>=0 and abs(sum(physical)-1)<1e-12,'complete physical binary law')
 near(path,physical,'raw centers close to physical law',1e-12)
 repaired=amended['rows'][len(results)]
 near(physical,repaired['physical_path_law'],'amendment physical path outcomes',1e-12)
 need(repaired['selected_path']==selected,'selected physical path flag')
 if exact_law is not None:need([str(x) for x in exact_law]==repaired['exact_control_law']['path_probabilities'],'exact rational probabilities')
 clock_log=analyzer+math.log10(physical[selected]);pair_log=2*clock_log
 near([analyzer,clock_log,pair_log],[repaired['log10_probability_reaching_analyzer'],repaired['log10_probability_selected_clock'],repaired['log10_probability_independent_identical_pair']],'complete single/pair cost',1e-9)
 pair_eta=1-fidelity**2-float(1-fidelity-eta)**2
 near(pair_eta,corrected['pair_off_band'],'complete pair off-band including cross-sign zero modes')
 need(-1e-10<=eta<=.002+1e-10 and -1e-10<=pair_eta<=.01+1e-10,'frozen certificates on actual returned state')
 pair_moments={n:complex(np.trace(arc@np.linalg.matrix_power(unit,n)))**2 for n in (1,8)}
 for n,zp in pair_moments.items():
  pr=corrected['pair_records'][str(n)]
  near([(1+zp.real)/2,(1-zp.real)/2,(1+zp.imag)/2,(1-zp.imag)/2],[pr[a] for a in ('x_plus','x_minus','y_plus','y_minus')],'complete pair Ramsey law')
 def cosine_band(real):
  lo=max(-1,min(1,(real-4e-5-.01)/.99));hi=max(-1,min(1,(real+4e-5+.01)/.99))
  return math.acos(hi),math.acos(lo)
 lo,hi=cosine_band(pair_moments[1].real);a8,b8=cosine_band(pair_moments[8].real)
 components=[]
 # Each monotone half-wave on [0,pi] is enumerated, without predicted frequency.
 for h in range(8):
  a0,b0=(h*math.pi+a8,h*math.pi+b8) if h%2==0 else ((h+1)*math.pi-b8,(h+1)*math.pi-a8)
  aa,bb=max(lo,a0/8),min(hi,b0/8)
  if aa<=bb:components.append([aa,bb])
 near(components,corrected['pair_reader']['frequency_components'],'all compatible stride branches')
 need(len(components)==1 and components[0][0]<=2*expected<=components[0][1],'unique bounded pair comparison')
 need(corrected['R_events_per_first_binary_setting']==2*2**k and corrected['R_events_per_stride8_binary_setting']==16*2**k,'actual pair event cost per setting')
 scale=511*math.sqrt(5)/2;energy_bounds=[[scale*x for x in c] for c in components]
 near(energy_bounds,corrected['pair_energy_interval_keV'],'one limiting calibration')
 if k==6:need(energy_bounds[0][1]<1024.4508487023415,'held-out exclusion after reader')
 results.append({'M':m,'turns':t,'direction':d,'law':law,'visibility':float(2*path[0]-1),'direct_visibility':direct,'normalized_crossblock':normcross,'physical_path_law':physical,'positive_fidelity':fidelity,'off_band_weight':eta,'pair_off_band_weight':pair_eta,'pair_energy_interval':energy_bounds,'log10_analyzer':analyzer,'log10_selected_clock':clock_log,'log10_pair':pair_log})
need(len(results)==16,'full case grid')
print(json.dumps({'k':k,'verdict':'PASS','checks':checks,'counts':counts,'cases':results,'reference':reference,'input_evidence_sha256':hashlib.sha256(evidence.read_bytes()).hexdigest(),'same_group_only':True,'producer_imports':False,'seconds':time.monotonic()-start,'finite_shot_precision_claim':False},indent=2))
