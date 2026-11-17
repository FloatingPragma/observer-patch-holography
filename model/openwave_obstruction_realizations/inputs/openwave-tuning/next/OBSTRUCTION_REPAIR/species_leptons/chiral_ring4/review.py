"""Independent source reconstruction and full-state residual review."""
from pathlib import Path
from itertools import combinations
import json,math
import numpy as np
from scipy import sparse
ROOT=Path(__file__).parent
M=16;L=4
parity=np.array([(-1)**w.bit_count()for w in range(1<<M)])
X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1,-1]).astype(complex)

def one(m):
 A=np.zeros((16,16),complex);C=np.zeros_like(A)
 P=(np.eye(2)+X)/2;Q=(np.eye(2)-X)/2
 for x in range(4):
  a=4*x
  A[a:a+2,a:a+2]+=8*Z;A[a+2:a+4,a+2:a+4]+=8*Z
  T=-4*(Z+1j*Y)
  A[a:a+2,a+2:a+4]+=T;A[a+2:a+4,a:a+2]+=T.conj().T
  A[a:a+2,a+2:a+4]+=m*P@Z@Q
  A[a+2:a+4,a:a+2]+=m*Q@Z@P
  for wall in range(2):
   a=4*x+2*wall;b=4*((x+1)%4)+2*wall
   T=-4*(Z+1j*X)
   if x<3:A[a:a+2,b:b+2]+=T;A[b:b+2,a:a+2]+=T.conj().T
   else:C[a:a+2,b:b+2]+=T
 return A,C

def index(words):
 idx=np.full(1<<M,-1,dtype=int);idx[words]=np.arange(len(words));return idx

def lift(H,words):
 idx=index(words);rr=[];cc=[];vv=[]
 for a,b in zip(*np.nonzero(H)):
  a,b=int(a),int(b)
  sel=np.flatnonzero((words>>b)&1)
  removed=words[sel]^(1<<b)
  valid=((removed>>a)&1)==0;sel=sel[valid];removed=removed[valid]
  sg=parity[words[sel]&((1<<b)-1)]*parity[removed&((1<<a)-1)]
  rr.extend(idx[removed|(1<<a)]);cc.extend(sel);vv.extend(sg*H[a,b])
 return sparse.csr_matrix((vv,(rr,cc)),shape=(len(words),len(words)))

def fields(words,flux):
 N=int(words[0]).bit_count();q=np.array([[int(w&((1<<(4*(x+1)))-1)).bit_count()+2*(3-x)-N for x in range(4)]for w in words])
 e=np.sum((flux[:,None,None]+q[None,:,:])**2,axis=2)/4
 return q,e

def action(A,C,e,V):
 f,n=e.shape;k=V.shape[1];v=V.reshape(f,n,k)
 out=np.empty_like(v)
 for j in range(f):
  out[j]=A@v[j]+e[j,:,None]*v[j]
  if j:out[j]+=C@v[j-1]
  if j+1<f:out[j]+=C.conj().T@v[j+1]
 return out.reshape(f*n,k)

def create(low,high,coef):
 idx=index(high);out=sparse.lil_matrix((len(high),len(low)),dtype=complex)
 for a,c in enumerate(coef):
  if c==0:continue
  cols=np.flatnonzero(((low>>a)&1)==0)
  out[idx[low[cols]|(1<<a)],cols]=c*parity[low[cols]&((1<<a)-1)]
 return out.tocsr()

rows=json.loads((ROOT/'RESULT.json').read_text());pot=json.loads((ROOT/'POTENTIAL_BOUNDS.json').read_text())
assert len(rows)==4
reviews=[]
for row in rows:
 name=row['case'];m=row['mass'];K=row['K']
 p=np.load(ROOT/(name+'_PREPARATION.npz'));r=np.load(ROOT/(name+'_RESOLVENTS.npz'))
 words=p['words'];flux=p['flux'];psi=p['amplitudes'];cw=r['words'];sol=r['solutions'];rhs=r['rhs']
 assert len(words)==math.comb(16,8) and len(cw)==math.comb(16,7)
 assert len(flux)==2*K+1 and len(psi)==len(words)*len(flux)
 assert all(int(w).bit_count()==8 for w in words)
 assert np.max(abs(p['probabilities'].ravel()-abs(psi)**2))<1e-14
 A,C=one(m);AN=lift(A,words);CN=lift(C,words);AC=lift(A,cw);CC=lift(C,cw)
 q,e=fields(words,flux);qc,ec=fields(cw,flux)
 assert np.array_equal(q,p['Gauss_prefix_charge']) and np.array_equal(qc,r['Gauss_prefix_charge'])
 HP=action(AN,CN,e,psi[:,None]).ravel()
 energy=float(np.vdot(psi,HP).real/np.vdot(psi,psi).real)
 assert abs(energy-row['variational_energy'])<1e-10
 rr=action(AC,CC,ec,sol)-row['variational_energy']*sol-rhs
 rr=rr.reshape(len(flux),len(cw),4)
 outside=np.vstack([CC.conj().T@sol.reshape(len(flux),len(cw),4)[0],CC@sol.reshape(len(flux),len(cw),4)[-1]])
 for a,test in enumerate(row['complete_resolvent_checks']):
  total=math.sqrt(float(np.linalg.norm(rr[:,:,a])**2+np.linalg.norm(outside[:,a])**2))
  assert total<=test['full_infinite_residual_upper']
  assert test['full_infinite_solution_error_upper']>=total/row['charged_infinite_resolvent_gap_lower']
 # Reconstruct both currents. The particle-hole sign uses a closed
 # parity formula rather than a sequence of annihilations.
 high=np.array([sum(1<<a for a in v)for v in combinations(range(16),9)],dtype=int)
 highidx=index(high);phase=[];target=[]
 for word in cw:
  w=int(word);n_pairs=sum(((w>>(2*j))&3)==3 for j in range(8))
  exponent=sum((a^1)for a in range(16)if(w>>a)&1)-n_pairs
  phase.append((-1)**exponent)
  comp=((1<<16)-1)^w
  target.append(highidx[((comp&0x5555)<<1)|((comp&0xAAAA)>>1)])
 P=sparse.csr_matrix((phase,(target,np.arange(len(cw)))),shape=(len(high),len(cw)))
 E=np.zeros((16,2),complex);E[:2]=(np.eye(2)+X)/2;E[2:4]=(np.eye(2)-X)/2
 computed=[]
 for a in range(2):
  v=(create(words,high,E[:,a])@psi.reshape(len(flux),len(words)).T).T
  computed.append((P.conj().T@v[::-1].T).T.ravel())
 for a in range(2):
  computed.append((create(cw,words,E[:,a]).conj().T@psi.reshape(len(flux),len(words)).T).T.ravel())
 current_error=float(np.max(abs(np.column_stack(computed)-rhs)))
 assert current_error<1e-12
 CAR=rhs[:,:2].conj().T@rhs[:,:2]+(rhs[:,2:].conj().T@rhs[:,2:]).T
 assert np.max(abs(CAR-np.eye(2)))<1e-10
 assert len(row['every_sector'])==17
 for sector in row['every_sector']:
  N=sector['N'];b=pot['all_number_sectors'][N]['best_Barta_infinite_rotor_lower']-(0 if N in[0,16]else 4*abs(m))
  assert abs(b-sector['infinite_rotor_lower'])<1e-12
  if N!=8:assert b>row['variational_energy_upper']
 assert not row['exact_neutral_vacuum_state_certified']
 assert not row['true_vacuum_propagator_certified']
 reviews.append({'case':name,'source_module_imported':False,'independent_full_operator_action_pass':True,
  'independent_particle_hole_current_error':current_error,'all_infinite_residual_bounds_pass':True,
  'all_sectors_compared':17,'full_preparation_records':len(psi),'full_charged_solution_coefficients':int(sol.size)})

# Independent literal full one-body checks of every momentum and band.
A,C=one(0);banderr=0
for j in range(64):
 theta=2*np.pi*j/64;literal=np.linalg.eigvalsh(A+C*np.exp(1j*theta)+C.conj().T*np.exp(-1j*theta))
 theory=[]
 for k in ((theta+2*np.pi*np.arange(4))/4):
  W=1-np.cos(k);ep=math.sqrt(2*W+.5+.5*math.sqrt(1+4*W*W));em=math.sqrt(max(0,(3*W*W+2*W)/(ep*ep)))
  theory.extend([-8*ep,-8*em,8*em,8*ep])
 banderr=max(banderr,float(max(abs(literal-np.sort(theory)))))
assert banderr<1e-11
receipt={'cases':reviews,'literal_full_band_formula_error':banderr,'momentum_samples':64,
 'all16CAR_modes_retained':True,'sparse_neutral_order_not_promoted':True,
 'warm_seed_failed_second_Ritz_control_retained':(ROOT/'warm_seed_control/SCOPE.json').exists()}
(ROOT/'REVIEW.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
