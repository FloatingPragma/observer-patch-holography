"""Independent complete finite-matrix check; no regional Gaussian-mixture shortcut."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
from pathlib import Path
import sys,itertools,json,argparse
import numpy as np
from scipy.linalg import expm
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from audit import ALPHA,BETA,dagger

parser=argparse.ArgumentParser();parser.add_argument('--order',type=int,choices=[2,4],default=2)
args=parser.parse_args()
L=5;N=L**3;R=1.
grid=np.array(list(itertools.product(range(L),repeat=3)))
coords=grid-L//2
k=2*np.pi*(grid+.5)/L
F=np.exp(1j*coords@k.T)/np.sqrt(N)
if args.order==4:
    t=np.sin(k/2)**2
    la=np.sum(4*t*(1+t/3),axis=1)
    d=np.sin(k)*(1+2*t/3)
else:
    la=2*np.sum(1-np.cos(k),axis=1);d=np.sin(k)
Ak=sum(d[:,j,None,None]*ALPHA[j] for j in range(3))
D=np.einsum('xk,kab,yk->xayb',F,Ak,F.conj()).reshape(4*N,4*N)
Lap=np.kron((F*la)@F.conj().T,np.eye(4))
W=Lap-D@D
Beta=np.kron(np.eye(N),BETA)
sigma=np.cos(2*np.pi*coords[:,0]/L)+.4*np.cos(2*np.pi*coords[:,1]/L)
sig=np.repeat(sigma,4)
selected=np.sum(coords**2,axis=1)<=1
sel=np.flatnonzero(np.repeat(selected,4))
mean=float(np.mean(sigma[selected]))
# xi=-mean*x on the ball plus its stencil neighbourhood. Periodic outside
# continuation is retained explicitly; no assertion of exact continuum locality.
A=np.zeros((N,N),complex)
for j in range(3):
    deriv=(F*(1j*np.sin(k[:,j])))@F.conj().T
    xi=-mean*coords[:,j]
    A+=.5*(xi[:,None]*deriv+deriv*xi[None,:])
A=np.kron(A,np.eye(4))

def hamiltonian(eps):
    v=np.exp(-eps*sig)
    de=.5*(v[:,None]*D+D*v[None,:])
    qe=de@de+v[:,None]*W*v[None,:]
    return de+R**3*Beta@qe@qe

def projector(H):
    e,u=np.linalg.eigh(H)
    return (u*(e<0))@dagger(u)

def entropy(C):
    e=np.clip(np.linalg.eigvalsh(C[np.ix_(sel,sel)]),1e-13,1-1e-13)
    return float(-np.sum(e*np.log(e)+(1-e)*np.log(1-e)))

H0=hamiltonian(0)
P=projector(H0)
ev,vec=np.linalg.eigh(P[np.ix_(sel,sel)])
ec=np.clip(ev,1e-13,1-1e-13)
K=(vec*np.log((1-ec)/ec))@dagger(vec)
dD=-.5*(sig[:,None]*D+D*sig[None,:])
dQ=D@dD+dD@D-sig[:,None]*W-W*sig[None,:]
dH=dD+R**3*Beta@(Lap@dQ+dQ@Lap)
e,u=np.linalg.eigh(H0);oc=(e<0).astype(float)
df=oc[:,None]-oc[None,:];de=e[:,None]-e[None,:]
ratio=np.zeros_like(df);mask=df!=0;ratio[mask]=df[mask]/de[mask]
dP=u@(ratio*(dagger(u)@dH@u))@dagger(u)
dCshape=A@P-P@A
metric=float(np.trace(K@dP[np.ix_(sel,sel)]).real)
shape=float(np.trace(K@dCshape[np.ix_(sel,sel)]).real)
slopes=[]
for eps in (1e-4,5e-5):
    plus=expm(eps*A);minus=expm(-eps*A)
    Pp=projector(hamiltonian(eps));Pm=projector(hamiltonian(-eps))
    slopes.append(dict(epsilon=eps,metric=(entropy(Pp)-entropy(Pm))/(2*eps),
                 shape=(entropy(plus@P@dagger(plus))-entropy(minus@P@dagger(minus)))/(2*eps),
                 joint=(entropy(plus@Pp@dagger(plus))-entropy(minus@Pm@dagger(minus)))/(2*eps)))
epsilon=1e-5
hamerr=float(np.max(abs((hamiltonian(epsilon)-hamiltonian(-epsilon))/(2*epsilon)-dH)))
# Exact coordinate-change covariance checks sign and the CAR-algebra convention.
U=expm(.01*A)
Pcoord=projector(dagger(U)@H0@U)
gauge_error=float(np.max(abs(U@Pcoord@dagger(U)-P)))
out=dict(L=L,source_difference_order=args.order,full_single_particle_modes=4*N,regional_modes=len(sel),
         min_eigenvalue_W=float(np.linalg.eigvalsh(W).min()),
         flat_Q_minus_Lap=float(np.max(abs(D@D+W-Lap))),
         generator_antihermitian_error=float(np.max(abs(A+dagger(A)))),
         generator_is_second_quantizable_number_preserving=True,
         operator_tangent_error=hamerr,metric_entropy_tangent=metric,
         geometric_algebra_entropy_tangent=shape,joint_entropy_tangent=metric+shape,
         finite_entropy_differences=slopes,coordinate_covariance_error=gauge_error,
         physical_diffeomorphism_equivalence_proved=False)
out['pass']=bool(hamerr<1e-5 and gauge_error<1e-11 and
                 max(abs(slopes[-1]['metric']-metric),abs(slopes[-1]['shape']-shape),
                     abs(slopes[-1]['joint']-metric-shape))<1e-5)
name='DENSE_CONTROL.json' if args.order==2 else 'DENSE_CONTROL_ORDER4.json'
(HERE/name).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2),flush=True)
if not out['pass']:raise SystemExit(1)
