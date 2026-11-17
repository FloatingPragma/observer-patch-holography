"""Actual finite source entropy checks for both independent source vertices.

This is a small Galerkin component check, not a continuum/locality claim.
All retained source modes are diagonalized; omitted Galerkin complement
is explicitly filled by the local beta reference in the observer control.
"""
import itertools,json
import numpy as np
from scipy.linalg import expm
from kernels import s,cell,vertices

L=3;h=.5;R=1.;twist=.5
labels=np.array(list(itertools.product(range(-1,2),repeat=3)))
p=2*np.pi*(labels+twist)/(L*h);N=len(p);dimension=4*N
D=np.zeros((dimension,dimension),complex)
BETA=np.kron(np.eye(N),s.BETA)
Sigma=np.zeros_like(D)
for i in range(N):
    D[4*i:4*i+4,4*i:4*i+4]=np.einsum('j,jab->ab',p[i],s.ALPHA)
    for j in range(N):
        diff=labels[i]-labels[j]
        if abs(diff[0])==1 and np.all(diff[1:]==0):Sigma[4*i:4*i+4,4*j:4*j+4]=s.I4/2

def flat(D):
    de,U=np.linalg.eigh(D)
    D4=(U*de**4)@U.conj().T
    f=(U/(np.sqrt(1+R**6*de**6)))@U.conj().T
    H=(D+R**3*BETA@D4)@f
    return (H+H.conj().T)/2

def occupied(H):
    e,U=np.linalg.eigh(H)
    return (U*(e<0))@U.conj().T

H0=flat(D);P0=occupied(H0)
cells=np.array([c for c in itertools.product(range(-1,2),repeat=3) if sum(x*x for x in c)<=1])
w,_=cell.window(p,h)
overlap=w[:,None]*np.exp(-1j*p@(cells*h).T)/np.sqrt(N)
F=np.einsum('pa,ij->piaj',overlap,s.I4).reshape(dimension,4*len(cells))
PBfull=(np.eye(dimension)-BETA)/2
PBregional=np.kron(np.eye(len(cells)),cell.PB)
def regional(P):return PBregional+F.conj().T@(P-PBfull)@F
C0=regional(P0);ce,U=np.linalg.eigh(C0);ce=np.clip(ce,1e-13,1-1e-13)
K=(U*np.log((1-ce)/ce))@U.conj().T
eps=1e-5;rows=[]
for kind in ['coframe','lapse']:
    states=[]
    for t in [eps,-eps]:
        if kind=='coframe':
            v=expm(-t*Sigma/2);H=flat(v@D@v)
        else:
            lapse=expm(t*Sigma);H=(lapse@H0+H0@lapse)/2
        states.append(occupied(H))
    finiteP=(states[0]-states[1])/(2*eps)
    analytic=np.zeros_like(P0)
    for i in range(N):
        for j in range(N):
            factor=Sigma[4*i,4*j].real
            if factor==0:continue
            dm,dn=vertices(p[j],p[i]-p[j],R,0.)
            analytic[4*i:4*i+4,4*j:4*j+4]=factor*(dm if kind=='coframe' else dn)
    finiteS=(s.entropy(regional(states[0]))-s.entropy(regional(states[1])))/(2*eps)
    first=float(np.trace(K@F.conj().T@analytic@F).real)
    rows.append(dict(kind=kind,projector_derivative_error=float(np.linalg.norm(finiteP-analytic)),finite_entropy_derivative=finiteS,analytic_entropy_derivative=first,entropy_derivative_error=abs(first-finiteS)))
out=dict(status='COMPLETE_FINITE_GALERKIN_SOURCE_VERTEX_CONTROL',source_modes=dimension,observer_modes=4*len(cells),epsilon=eps,omitted_complement='Declared local beta reference, not a claim that the true continuum complement vanishes.',max_cell_gram_eigenvalue=float(np.linalg.eigvalsh(F.conj().T@F)[-1]),rows=rows)
with open('VERTEX_DENSE_CONTROL.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2),flush=True)
assert max(r['projector_derivative_error'] for r in rows)<1e-6
assert max(r['entropy_derivative_error'] for r in rows)<1e-6
