"""Independent full finite-carrier projector/charge/energy control."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
from pathlib import Path
import sys,json,itertools
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from audit import ALPHA,BETA,dagger
L=4;N=L**3;R=1.5;W=2.
grid=np.array(list(itertools.product(range(L),repeat=3)))
k=2*np.pi*(grid+.5)/L;t=np.sin(k/2)**2
d=np.sin(k)*(1+2*t/3);lap=np.sum(4*t*(1+t/3),axis=1)
Hk=sum(d[:,j,None,None]*ALPHA[j] for j in range(3))+R**3*lap[:,None,None]**2*BETA
E,U=np.linalg.eigh(Hk);Pm=(U*(E<0)[:,None,:])@dagger(U);Pp=np.eye(4)-Pm
kw=(k+np.pi)%(2*np.pi)-np.pi;f=np.exp(-W**2*np.sum(kw*kw,axis=1)/4)
F=np.exp(1j*grid@k.T)/np.sqrt(N)
H=np.einsum('xk,kab,yk->xayb',F,Hk,F.conj()).reshape(4*N,4*N)
P=np.einsum('xk,kab,yk->xayb',F,Pm,F.conj()).reshape(4*N,4*N)
all_modes=[];energies=[]
for proj,sign in [(Pp,1),(Pm,-1)]:
    raw=f[:,None,None]*proj
    gram=np.einsum('kas,kat->st',raw.conj(),raw)
    e,u=np.linalg.eigh(gram);normalizer=(u/np.sqrt(e))@dagger(u)
    modes=raw@normalizer
    X=np.einsum('xk,kas->xas',F,modes).reshape(4*N,4)
    all_modes.append(X);energies.append(sign*np.trace(dagger(X)@H@X).real)
Xp,Xh=all_modes;C=P-Xh@dagger(Xh)+Xp@dagger(Xp)
ev=np.linalg.eigvalsh(C)
out=dict(L=L,source_R=R,source_difference_order=4,full_modes=4*N,
         particle_gram_error=float(np.max(abs(dagger(Xp)@Xp-np.eye(4)))),
         hole_gram_error=float(np.max(abs(dagger(Xh)@Xh-np.eye(4)))),
         cross_gram_error=float(np.max(abs(dagger(Xp)@Xh))),
         particle_in_empty_error=float(np.max(abs(P@Xp))),
         hole_in_filled_error=float(np.max(abs(P@Xh-Xh))),
         projector_error=float(np.max(abs(C@C-C))),
         hermiticity_error=float(np.max(abs(C-dagger(C)))),
         trace_charge_change=float(np.trace(C-P).real),vacuum_rank=int(round(np.trace(P).real)),
         excited_rank=int(np.sum(ev>.5)),min_eigenvalue=float(ev.min()),max_eigenvalue=float(ev.max()),
         excitation_energy_by_modes=float(sum(energies)),
         excitation_energy_by_complete_covariance=float(np.trace(H@(C-P)).real))
out['pass']=bool(max(out[x] for x in ['particle_gram_error','hole_gram_error','cross_gram_error','particle_in_empty_error','hole_in_filled_error','projector_error','hermiticity_error'])<1e-10 and abs(out['trace_charge_change'])<1e-10 and out['excited_rank']==out['vacuum_rank'] and abs(out['excitation_energy_by_modes']-out['excitation_energy_by_complete_covariance'])<1e-8)
(HERE/'COMPLETE_CONTROL.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2),flush=True)
if not out['pass']:raise SystemExit(1)
