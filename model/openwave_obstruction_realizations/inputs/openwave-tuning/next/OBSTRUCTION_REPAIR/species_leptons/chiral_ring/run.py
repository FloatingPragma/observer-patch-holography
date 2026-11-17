"""All ring sectors, complete small spectra and full-space large resolvents."""
from fractions import Fraction as F
import json,math,time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh,spsolve
import source as s
ROOT=Path(__file__).parent
def save(n,x):(ROOT/n).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def norm(x):return float(np.linalg.norm(x))
def low(a,d,c):return (a+d-math.sqrt((a-d)**2+4*c*c))/2

def complete(H):
 D,Q=eigh(H.toarray(),driver='evd');n=len(D);u=np.finfo(float).eps/2
 R=H@Q-Q*D;gram=Q.conj().T@Q-np.eye(n)
 gamma=16*n*u/(1-16*n*u)
 eta=norm(gram)+gamma*norm(Q)**2+1e-12
 residual=norm(R)+gamma*(norm(H.data)*norm(Q)+norm(Q*D))+1e-10
 e=residual/math.sqrt(1-eta)
 return D,Q,{'complete':True,'finite_eigenvalue_error_bound':e,'orthogonality_bound':eta,'residual_bound':residual}

def eigen(H,K):
 if K<=8:return complete(H)
 n=len(H.toarray()) if H.shape[0]<64 else H.shape[0]
 if n<64:return complete(H)
 # Fixed deterministic seed vector, independent of physical targets.
 v=np.cos(np.arange(n)+1)+np.sin((np.arange(n)+1)/3)
 D,Q=eigsh(H,k=min(24,n-2),which='SA',tol=1e-11,maxiter=20000,v0=v)
 order=np.argsort(D);D=D[order];Q=Q[:,order]
 return D,Q,{'complete':False,'eigenpairs':len(D),'ground_residual':norm(H@Q[:,0]-D[0]*Q[:,0]),
  'finite_ground_order_certified':False}

intervals=json.loads((ROOT/'A0_INTERVALS.json').read_text())['intervals']
def tail(N,K,m,E0,error):
 a=sum(float(F(x[0]))for x in intervals[:N])-2*abs(m)-1e-9
 c=s.coupling_norm_bound(N)+1e-9
 d=s.electric_outside_lower(N,K)+a-2*c
 return {'A_lower':a,'coupling_upper':c,'electric_outside_lower':s.electric_outside_lower(N,K),
  'discarded_lower':d,'infinite_ground_lower':low(E0-error,d,c),
  'finite_ground_upper':E0+error,'input_finite_ground_certified':bool(error>0)}

results=[];spectra=[];records=[]
for label,m in [('zero',0.),('register_electron',6.0300795404423236e-5)]:
 for K in [4,8,16]:
  started=time.time();sectors=[];cached={}
  for N in range(5):
   H,A,C,E,q=s.block(N,K,m);D,Q,cert=eigen(H,K);cached[N]=(H,D,Q,cert)
   row={'N':N,'dimension':H.shape[0],'ground':float(D[0]),'finite_check':cert}
   row['tail']=tail(N,K,m,float(D[0]),cert.get('finite_eigenvalue_error_bound',0))
   sectors.append(row)
  for N in range(5,9):
   row=dict(sectors[8-N]);row['N']=N;row['particle_hole_of']=8-N
   sectors.append(row)
  assert min(sectors,key=lambda x:x['ground'])['N']==4
  H,D,Q,cert=cached[4];E0=float(D[0]);state=Q[:,0];n=math.comb(8,4)
  wf=state.reshape(2*K+1,n);cos=float(np.vdot(wf[1:],wf[:-1]).real)
  flux=np.arange(-K,K+1);words=s.loop.words(8,4)
  _,_,_,E,q=s.block(4,K,m)
  electric=((flux[:,None]+q[None,:])**2+flux[:,None]**2)/4
  electric_mean=float(np.sum(abs(wf)**2*electric))
  density=np.array([[np.vdot(wf[j],wf[k])for k in range(2*K+1)]for j in range(2*K+1)])
  interval=2*np.pi/3;diff=flux[None,:]-flux[:,None]
  kernel=np.full(diff.shape,interval/np.pi,dtype=float)
  mask=diff!=0;kernel[mask]=np.sin(interval*diff[mask])/(np.pi*diff[mask])
  localization_probability=float(np.sum(density*kernel).real)
  Hm,Dm,Qm,certm=cached[3];Hp,*_=s.block(5,K,m)
  PH=s.particle_hole(3,K)
  Rp=Hp-PH@Hm@PH.conj().T
  assert max(abs(Rp.data),default=0)<1e-12
  Dp=Dm;Qp=PH@Qm
  vp=np.column_stack([sparse.kron(sparse.eye(2*K+1),s.creation(4,E[:,a]))@state for a in range(4)])
  vm=np.column_stack([sparse.kron(sparse.eye(2*K+1),s.creation(3,E[:,a]).conj().T)@state for a in range(4)])
  Ip=sparse.eye(Hp.shape[0],format='csr');Im=sparse.eye(Hm.shape[0],format='csr')
  solP=spsolve(Hp-E0*Ip,vp);solM=spsolve(Hm-E0*Im,vm)
  G=-vp.conj().T@solP+(vm.conj().T@solM).T
  MG=-np.linalg.inv(G)
  coeff={name:np.array([[np.trace(P@MG[2*x:2*x+2,2*y:2*y+2])/2 for y in range(2)]for x in range(2)])for name,P in [('I',np.eye(2)),('X',s.loop.X),('Y',s.loop.Y),('Z',s.loop.Z)]}
  CAR=vp.conj().T@vp+(vm.conj().T@vm).T
  if K<=8:
   ap=Qp.conj().T@vp;am=Qm.conj().T@vm
   rp=np.einsum('ia,ib->iab',ap.conj(),ap);rm=np.einsum('ib,ia->iab',am.conj(),am)
   assert np.max(abs(sum(rp)+sum(rm)-np.eye(4)))<1e-10
   def pack(vals,res):return [{'gap':float(e-E0),'residue_real':r.real.tolist(),'residue_imag':r.imag.tolist()}for e,r in zip(vals,res)]
   spectra.append({'mass':label,'K':K,'particle':pack(Dp,rp),'hole':pack(Dm,rm),'complete':True})
  else:
   spectra.append({'mass':label,'K':K,'retained_low_gaps':(Dm-E0).tolist(),'complete_spectrum':False,
    'full_space_resolvent_used':True,'particle_dimension':Hp.shape[0],'hole_dimension':Hm.shape[0]})
  lows={r['N']:r['tail']['infinite_ground_lower']for r in sectors};highs={r['N']:r['tail']['finite_ground_upper']for r in sectors}
  certified=K<=8 and all(lows[N]>highs[4]for N in range(9)if N!=4)
  item={'mass_label':label,'mass':m,'K':K,'all_sectors':sectors,'vacuum_number':4,
   'infinite_sector_certified':certified,'vacuum_energy':E0,'holonomy_cos_mean':cos,
   'electric_energy_mean':electric_mean,'probability_low_W_at_most_half':localization_probability,
   'edge_flux_probability':float(norm(wf[0])**2+norm(wf[-1])**2),
   'charged_reference_gap':float(Dm[0]-E0),'coarse_gap_interval':[lows[3]-highs[4],highs[3]-lows[4]],
   'gap_interval_certified':K<=8,'CAR_sum_error':float(np.max(abs(CAR-np.eye(4)))),
   'particle_hole_intertwiner_error':float(max(abs(Rp.data),default=0)),
   'full_particle_resolvent_residual':norm((Hp-E0*Ip)@solP-vp),
   'full_hole_resolvent_residual':norm((Hm-E0*Im)@solM-vm),
   'minus_G_inverse_coefficients':{a:{'real':v.real.tolist(),'imag':v.imag.tolist()}for a,v in coeff.items()},
   'Green_real':G.real.tolist(),'Green_imag':G.imag.tolist(),
   'physical_lepton_pole':False,'elapsed_seconds':time.time()-started}
  results.append(item)
  records.append({'mass':label,'K':K,'probabilities':[{'flux':int(flux[j//n]),'occupation':int(words[j%n]),'probability':float(p)}for j,p in enumerate(abs(state)**2)]})
  save('RESULT.json',results);save('BOUNDARY_SPECTRA.json',spectra);save('VACUUM_RECORDS.json',records)
  print(json.dumps({'mass':label,'K':K,'E0':E0,'gap':item['charged_reference_gap'],'cos':cos,
   'P_low_W_le_half':localization_probability,'beta_matrix':coeff['Z'].real.tolist(),
   'CAR_error':item['CAR_sum_error'],'sector_certified':certified,'seconds':item['elapsed_seconds']}),flush=True)
