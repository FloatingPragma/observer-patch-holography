"""Full-ball versus exact G1 reduction, before attempting larger balls."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
from pathlib import Path
import sys,json,itertools,math
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from audit import ALPHA,BETA,dagger
from group import full_small_bases

L=20;R=1.5;W=4.;radius2=9
grid=np.stack(np.meshgrid(*([np.arange(L)]*3),indexing='ij'),-1).reshape(-1,3)
k=2*np.pi*(grid+.5)/L;t=np.sin(k/2)**2
d=np.sin(k)*(1+2*t/3);lap=np.sum(4*t*(1+t/3),axis=1)
H=sum(d[:,j,None,None]*ALPHA[j] for j in range(3))+R**3*lap[:,None,None]**2*BETA
e,u=np.linalg.eigh(H);P=(u*(e<0)[:,None,:])@dagger(u)
ker=np.fft.ifftn(P.reshape(L,L,L,4,4),axes=(0,1,2))
cells=np.array([x for x in itertools.product(range(-3,4),repeat=3) if sum(z*z for z in x)<=radius2])
N=len(cells);delta=cells[:,None,:]-cells[None,:,:];ix=delta%L
C=(ker[ix[:,:,0],ix[:,:,1],ix[:,:,2]]*np.exp(1j*np.pi*delta.sum(-1)/L)[:,:,None,None]).transpose(0,2,1,3).reshape(4*N,4*N)
kw=(k+np.pi)%(2*np.pi)-np.pi;f=np.exp(-W*W*np.sum(kw*kw,axis=1)/4)
modes=[]
for proj in (np.eye(4)-P,P):
    raw=f[:,None,None]*proj;gram=np.einsum('kas,kat->st',raw.conj(),raw)
    ge,gu=np.linalg.eigh(gram);normalizer=(gu/np.sqrt(ge))@dagger(gu)
    mom=raw@normalizer
    position=np.fft.ifftn(mom.reshape(L,L,L,4,4),axes=(0,1,2))*math.sqrt(L**3)
    pos=cells%L
    X=position[pos[:,0],pos[:,1],pos[:,2]]*np.exp(1j*np.pi*cells.sum(1)/L)[:,None,None]
    modes.append(X.reshape(4*N,4))
Xp,Xh=modes;change=Xp@dagger(Xp)-Xh@dagger(Xh)
B1,B2,orbits=full_small_bases(cells)
PG=B1@dagger(B1)+B2@dagger(B2)
C1=dagger(B1)@C@B1;C2=dagger(B2)@C@B2
G1=dagger(B1)@change@B1;G2=dagger(B2)@change@B2
def ent(M):
    ev=np.clip(np.linalg.eigvalsh(M),1e-13,1-1e-13)
    return float(-np.sum(ev*np.log(ev)+(1-ev)*np.log(1-ev)))
def first(C,change):
    ev,u=np.linalg.eigh(C);ev=np.clip(ev,1e-13,1-1e-13)
    return float(np.sum(np.log((1-ev)/ev)*np.diag(dagger(u)@change@u).real))
full_delta=ent(C+change)-ent(C);reduced_delta=2*(ent(C1+G1)-ent(C1))
out=dict(status='SMALL_BALL_CUBIC_REDUCTION_CONTROL',L=L,R=R,width=W,cells=N,
         full_modes=4*N,G1_multiplicity=B1.shape[1],G1_irrep_dimension=2,
         all_G1_copies_retained=True,parity_not_discarded=True,
         projector_idempotence_error=float(np.max(abs(PG@PG-PG))),
         orthonormal_error=float(np.max(abs(dagger(B1)@B1-np.eye(B1.shape[1])))),
         cross_component_error=float(np.max(abs(dagger(B1)@B2))),
         covariance_commutator=float(np.max(abs(C@PG-PG@C))),
         equal_irrep_blocks_error=float(np.max(abs(C1-C2))),
         cross_irrep_component_covariance=float(np.max(abs(dagger(B1)@C@B2))),
         matter_support_outside_G1=float(np.max(abs(change-PG@change@PG))),
         particle_support_outside_G1=float(np.max(abs(Xp-PG@Xp))),
         hole_support_outside_G1=float(np.max(abs(Xh-PG@Xh))),
         equal_matter_blocks_error=float(np.max(abs(G1-G2))),
         full_entropy_change=full_delta,reduced_entropy_change=reduced_delta,
         entropy_change_error=abs(full_delta-reduced_delta),
         full_first_derivative=first(C,change),reduced_first_derivative=2*first(C1,G1),
         geometry_response_reduced_to_G1=False)
out['pass']=bool(max(out[x] for x in ['projector_idempotence_error','orthonormal_error','cross_component_error','covariance_commutator','equal_irrep_blocks_error','cross_irrep_component_covariance','matter_support_outside_G1','particle_support_outside_G1','hole_support_outside_G1','equal_matter_blocks_error'])<1e-10 and out['entropy_change_error']<1e-8 and abs(out['full_first_derivative']-out['reduced_first_derivative'])<1e-8)
(HERE/'SMALL_CONTROL.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2),flush=True)
if not out['pass']:raise SystemExit(1)
