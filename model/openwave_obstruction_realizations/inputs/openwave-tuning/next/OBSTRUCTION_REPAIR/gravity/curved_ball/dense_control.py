"""Full finite-matrix check of the curved-ball entropy derivative."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
import numpy as np
import itertools,json
from pathlib import Path
from run import run,ALPHA,BETA,HERE

L=4;N=L**3;R=1;q=2*np.pi/L
grid=np.array(list(itertools.product(range(L),repeat=3)))
k=2*np.pi*(grid+.5)/L
F=np.exp(1j*grid@k.T)/np.sqrt(N)
lap=2*np.sum(1-np.cos(k),axis=1)
kin=sum(np.sin(k[:,j])[:,None,None]*ALPHA[j] for j in range(3))
K=np.einsum('xk,kab,yk->xayb',F,kin,F.conj()).reshape(4*N,4*N)
Lap=(F*lap)@F.conj().T
cells=np.array([c for c in itertools.product(range(-1,2),repeat=3) if sum(x*x for x in c)<=1])
idc=(cells%L)@np.array([L*L,L,1])
sel=(4*idc[:,None]+np.arange(4)).ravel()
means=np.mean(np.cos(q*cells),axis=0)
sigma=np.sum(means[None,:]-np.cos(q*grid),axis=1)


def hamiltonian(eps):
    v=np.exp(-eps*sigma)
    vv=np.repeat(v,4)
    B=v[:,None]*Lap*v[None,:]
    return .5*(vv[:,None]*K+K*vv[None,:])+np.kron(B@B,BETA)


def entropy(eps):
    H=hamiltonian(eps)
    E,U=np.linalg.eigh(H)
    P=(U*(E<0))@U.conj().T
    ev=np.clip(np.linalg.eigvalsh(P[np.ix_(sel,sel)]),1e-13,1-1e-13)
    return float(-np.sum(ev*np.log(ev)+(1-ev)*np.log(1-ev)))


S=np.diag(sigma)
dH=-.5*(np.repeat(sigma,4)[:,None]*K+K*np.repeat(sigma,4)[None,:])-np.kron(Lap@Lap@S+2*Lap@S@Lap+S@Lap@Lap,BETA)
eps=1e-5
tangent_error=float(np.max(np.abs((hamiltonian(eps)-hamiltonian(-eps))/(2*eps)-dH)))
bloch=run(4,1)
slopes=[dict(epsilon=e,derivative=(entropy(e)-entropy(-e))/(2*e)) for e in (1e-4,1e-5)]
out=dict(hamiltonian_tangent_absolute_error=tangent_error,bloch_derivative=bloch['curved_entropy_response'],
         full_matrix_entropy_derivatives=slopes,absolute_difference=abs(slopes[-1]['derivative']-bloch['curved_entropy_response']))
out['pass']=bool(tangent_error<1e-4 and out['absolute_difference']<1e-6)
(HERE/'DENSE_CONTROL.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
if not out['pass']:raise SystemExit(1)
