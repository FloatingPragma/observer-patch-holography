"""Independent real-space dC construction and modular-floor sensitivity."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
from pathlib import Path
import sys,json,itertools
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from audit import ALPHA,BETA,dagger
HERE=Path(__file__).resolve().parent
L=24;radius=4;R=1
n=np.arange(L)
grid=np.stack(np.meshgrid(n,n,n,indexing='ij'),-1).reshape(-1,3)
k=2*np.pi*(grid+.5)/L
kin=sum(np.sin(k[:,j])[:,None,None]*ALPHA[j] for j in range(3))
lap=2*np.sum(1-np.cos(k),axis=1)
H=kin+lap[:,None,None]**2*BETA
E,U=np.linalg.eigh(H);occ=E<0
P=(U*occ[:,None,:])@dagger(U)
cells=np.array([c for c in itertools.product(range(-radius,radius+1),repeat=3) if sum(x*x for x in c)<=radius**2])
N=len(cells);delta=cells[:,None,:]-cells[None,:,:];ix=delta%L
phase=np.exp(1j*np.pi*np.sum(delta,axis=-1)/L)


def regional_from_fourier(Pq,q):
    kernel=np.fft.ifftn(Pq.reshape(L,L,L,4,4),axes=(0,1,2))
    fac=phase*np.exp(1j*cells@q)[:,None]
    return (kernel[ix[:,:,0],ix[:,:,1],ix[:,:,2]]*fac[:,:,None,None]).transpose(0,2,1,3).reshape(4*N,4*N)


def dc(shift):
    idx=(grid+np.array(shift))%L;ids=idx[:,0]*L*L+idx[:,1]*L+idx[:,2]
    V=-.5*(kin+kin[ids])-(lap+lap[ids])[:,None,None]**2*BETA
    df=occ[ids].astype(float)[:,:,None]-occ.astype(float)[:,None,:]
    de=E[ids][:,:,None]-E[:,None,:]
    F=np.zeros_like(de);m=df!=0;F[m]=df[m]/de[m]
    dP=U[ids]@(F*(dagger(U[ids])@V@U))@dagger(U)
    Cq=regional_from_fourier(dP,2*np.pi*np.array(shift)/L)
    return .5*(Cq+Cq.conj().T)


C=regional_from_fourier(P,np.zeros(3));ev,vec=np.linalg.eigh(C)
means=np.mean(np.cos(2*np.pi*cells/L),axis=0)
dC=means.sum()*dc((0,0,0))-sum(dc(tuple(int(j==a) for j in range(3))) for a in range(3))
diagonal=np.real(np.diag(vec.conj().T@dC@vec))
near=np.flatnonzero(np.all((grid==0)|(grid==L-1),axis=1))
weights=np.zeros(4*N)
for idx in near:
    waves=np.einsum('x,ab->xab',np.exp(1j*cells@k[idx]),U[idx]).reshape(4*N,4)
    projected=vec.conj().T@waves
    weights+=(abs(projected[:,2])**2-abs(projected[:,1])**2)/(E[idx,2]-E[idx,1])/len(near)
records=[]
for floor in (1e-11,1e-12,1e-13,1e-14,1e-15):
    ec=np.clip(ev,floor,1-floor);logit=np.log((1-ec)/ec)
    records.append(dict(floor=floor,entropy_response=float(logit@diagonal),bulk_coefficient=float(logit@weights)))
ref=json.loads((HERE/'RESULT_L24_r4_R1.json').read_text())
out=dict(records=records,minimum_restricted_eigenvalue=float(ev.min()),maximum_restricted_eigenvalue=float(ev.max()),
         independent_response_difference=abs(records[2]['entropy_response']-ref['curved_entropy_response']),
         independent_bulk_difference=abs(records[2]['bulk_coefficient']-ref['same_ball_bulk_coefficient']))
out['checks']={
    'independent_realspace_response_agrees':out['independent_response_difference']<1e-7,
    'independent_bulk_agrees':out['independent_bulk_difference']<1e-7,
    'response_floor_spread_below_point_one_percent':(max(x['entropy_response'] for x in records)-min(x['entropy_response'] for x in records))/abs(records[2]['entropy_response'])<.001,
    'bulk_floor_spread_below_point_one_percent':(max(x['bulk_coefficient'] for x in records)-min(x['bulk_coefficient'] for x in records))/abs(records[2]['bulk_coefficient'])<.001,
}
out['pass']=all(out['checks'].values())
(HERE/'FLOOR_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
if not out['pass']:raise SystemExit(1)
