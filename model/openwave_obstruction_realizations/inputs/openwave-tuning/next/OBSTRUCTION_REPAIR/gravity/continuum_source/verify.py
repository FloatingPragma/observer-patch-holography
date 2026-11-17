"""Finite Clifford/entropy witnesses for the analytic continuum proof."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
from pathlib import Path
import json,math
import numpy as np
HERE=Path(__file__).resolve().parent
rng=np.random.default_rng(42819);n=8;R=.7
q,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
w,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
lam=np.geomspace(.4,8,n);A=(q*lam)@w.conj().T
D=np.block([[np.zeros((n,n)),A],[A.conj().T,np.zeros((n,n))]])
beta=np.diag([1.]*n+[-1.]*n);D2=D@D;D4=D2@D2
H=D+R**3*beta@D4
E,U=np.linalg.eigh(H);P=(U*(E<0))@U.conj().T;Pb=(np.eye(2*n)-beta)/2
Delta=P-Pb;dv,du=np.linalg.eigh(D2);tt=np.sqrt(dv);xx=R**3*tt**3
stable_squared=1/(2*np.sqrt(1+xx*xx)*(np.sqrt(1+xx*xx)+xx))
predicted=(du*stable_squared)@du.conj().T
selected=np.array([0,1,3,4,6,n,n+1,n+3,n+4,n+6]);other=np.array([x for x in range(2*n) if x not in selected])
C=P[np.ix_(selected,selected)];T=P[np.ix_(selected,other)]
ce=np.clip(np.linalg.eigvalsh(C),1e-15,1-1e-15)
S=float(-np.sum(ce*np.log(ce)+(1-ce)*np.log(1-ce)))
p=.75;Cp=2**p*(1+1/(math.e*(1-p)))
upper=float(Cp*np.sum(np.linalg.svd(T,compute_uv=False)**(2*p)))
samples=np.r_[np.geomspace(1e-14,.5,10000),1-np.geomspace(1e-14,.5,10000)]
hs=-samples*np.log(samples)-(1-samples)*np.log(1-samples)
scalar_ratio=float(np.max(hs/(Cp*(samples*(1-samples))**p)))
tables=[]
for cutoff in (4,8,16,32):
    grid=np.stack(np.meshgrid(*([np.arange(-cutoff,cutoff)+.5]*3),indexing='ij'),-1)
    t=2*np.pi*np.linalg.norm(grid,axis=-1);x=R**3*t**3
    sv=1/np.sqrt(2*np.sqrt(1+x*x)*(np.sqrt(1+x*x)+x))
    tables.append(dict(momentum_cube_cutoff=cutoff,absolute_D_gap=float(t.min()),
                       S1_partial_sum=float(4*np.sum(sv)),S1p5_partial_sum=float(4*np.sum(sv**1.5))))
out=dict(status='FINITE_ALGEBRA_CONTROL_NOT_A_SUBSTITUTE_FOR_PROOF',
         anticommutation_error=float(np.max(abs(beta@D+D@beta))),
         square_identity_relative_error=float(np.max(abs(H@H-D2-R**6*D4@D4))/np.max(abs(H@H))),
         exact_projector_modulus_error=float(np.max(abs(Delta@Delta-predicted))),
         compressed_projector_identity_error=float(np.max(abs(C-C@C-T@T.conj().T))),
         regional_entropy=S,proved_entropy_upper_bound=upper,
         maximum_sampled_scalar_bound_ratio=scalar_ratio,torus_schatten_diagnostic=tables,
         functional_calculus_and_Weyl_proof_in='PROOF.md',
         area_asymptotic_proved=False,shape_differentiability_proved=False,stress_renormalized=False)
out['pass']=bool(max(out['anticommutation_error'],out['square_identity_relative_error'],out['exact_projector_modulus_error'],out['compressed_projector_identity_error'])<1e-10 and S<=upper and scalar_ratio<=1)
(HERE/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2),flush=True)
if not out['pass']:raise SystemExit(1)
