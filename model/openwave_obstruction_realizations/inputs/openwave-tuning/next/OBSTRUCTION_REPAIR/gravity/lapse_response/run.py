"""Actual lapse entropy derivative on the same improved quartic source."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
from pathlib import Path
import sys,json,itertools,math,argparse,time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from audit import ALPHA,BETA,dagger

def source(L,R):
    grid=np.stack(np.meshgrid(*([np.arange(L)]*3),indexing='ij'),-1).reshape(-1,3)
    k=2*np.pi*(grid+.5)/L;t=np.sin(k/2)**2
    d=np.sin(k)*(1+2*t/3);lap=np.sum(4*t*(1+t/3),axis=1)
    H=sum(d[:,j,None,None]*ALPHA[j] for j in range(3))+R**3*lap[:,None,None]**2*BETA
    E=np.sqrt(np.sum(d*d,axis=1)+R**6*lap**4)
    P=(np.eye(4)-H/E[:,None,None])/2
    kernel=np.fft.ifftn(P.reshape(L,L,L,4,4),axes=(0,1,2))
    return dict(L=L,R=R,grid=grid,k=k,E=E,P=P,kernel=kernel,H=H)

def ball(src,r2):
    L=src['L'];r=int(math.sqrt(r2));P=src['P'];E=src['E'];grid=src['grid']
    cells=np.array([x for x in itertools.product(range(-r,r+1),repeat=3) if sum(z*z for z in x)<=r2])
    N=len(cells);delta=cells[:,None,:]-cells[None,:,:];ix=delta%L
    C=(src['kernel'][ix[:,:,0],ix[:,:,1],ix[:,:,2]]*np.exp(1j*np.pi*delta.sum(-1)/L)[:,:,None,None]).transpose(0,2,1,3).reshape(4*N,4*N)
    ev,u=np.linalg.eigh(C);ec=np.clip(ev,1e-13,1-1e-13)
    K=((u*np.log((1-ec)/ec))@dagger(u)).reshape(N,4,N,4);del C,u
    base=K.transpose(0,2,1,3);di=ix.reshape(-1,3);records=[]
    for shift in [(0,0,0),(1,0,0),(2,0,0)]:
        q=2*np.pi*np.array(shift)/L;ixq=(grid+np.array(shift))%L;ids=ixq[:,0]*L*L+ixq[:,1]*L+ixq[:,2]
        ph=np.exp(-1j*np.pi*delta.sum(-1)/L)*np.exp(1j*cells@q)[None,:]
        acc=np.zeros((L,L,L,4,4),complex)
        np.add.at(acc,(di[:,0],di[:,1],di[:,2]),(base*ph[:,:,None,None]).reshape(-1,4,4))
        kap=np.fft.fftn(acc,axes=(0,1,2)).reshape(L**3,4,4);del acc
        fac=(E[ids]-E)/(2*(E[ids]+E))
        dP=fac[:,None,None]*(P[ids]-P)
        response=np.einsum('kab,kba->',kap,dP)/L**3
        overlap=np.einsum('kab,kba->k',np.eye(4)-P[ids],P).real
        susceptibility=float(np.sum((E[ids]-E)**2/(2*(E[ids]+E))*overlap)/L**3)
        rec=dict(shift=list(shift),q2=float(q@q),lapse_entropy_kernel=[float(response.real),float(response.imag)],
                 energy_susceptibility=susceptibility,
                 chi_over_q4=None if not q@q else susceptibility/(q@q)**2)
        records.append(rec);print('kernel',r2,rec,flush=True)
    return dict(radius_squared=r2,cells=N,modes=records)

def vertex_control(src):
    L=src['L'];ids0=np.arange(min(23,L**3));qidx=(src['grid'][ids0]+[1,0,0])%L
    ids1=qidx[:,0]*L*L+qidx[:,1]*L+qidx[:,2]
    H0=src['H'][ids0];H1=src['H'][ids1]
    E0,U0=np.linalg.eigh(H0);E1,U1=np.linalg.eigh(H1)
    f0=(E0<0).astype(float);f1=(E1<0).astype(float)
    df=f1[:,:,None]-f0[:,None,:];de=E1[:,:,None]-E0[:,None,:]
    div=np.zeros_like(de);mask=df!=0;div[mask]=df[mask]/de[mask]
    general=U1@(div*(dagger(U1)@((H0+H1)/2)@U0))@dagger(U0)
    fac=(src['E'][ids1]-src['E'][ids0])/(2*(src['E'][ids1]+src['E'][ids0]))
    exact=fac[:,None,None]*(src['P'][ids1]-src['P'][ids0])
    return float(np.max(abs(general-exact)))

def dense_control(src,analytic):
    L=src['L'];N=L**3;coords=src['grid']-L//2
    F=np.exp(1j*coords@src['k'].T)/math.sqrt(N)
    H=np.einsum('xk,kab,yk->xayb',F,src['H'],F.conj()).reshape(4*N,4*N)
    lapse=np.repeat(np.cos(2*np.pi*coords[:,0]/L),4)
    selected=np.flatnonzero(np.repeat(np.sum(coords*coords,axis=1)<=1,4))
    def entropy(eps):
        n=np.exp(eps*lapse);He=.5*(n[:,None]*H+H*n[None,:])
        e,u=np.linalg.eigh(He);P=(u*(e<0))@dagger(u)
        c=np.clip(np.linalg.eigvalsh(P[np.ix_(selected,selected)]),1e-13,1-1e-13)
        return float(-np.sum(c*np.log(c)+(1-c)*np.log(1-c)))
    slopes=[dict(epsilon=e,actual_entropy_derivative=(entropy(e)-entropy(-e))/(2*e)) for e in (1e-4,5e-5)]
    expected=analytic['modes'][1]['lapse_entropy_kernel'][0]
    return dict(full_modes=4*N,analytic_derivative=expected,finite_entropy_differences=slopes,
                last_error=abs(slopes[-1]['actual_entropy_derivative']-expected))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--control',action='store_true');args=ap.parse_args()
    L=5 if args.control else 60;R=1.5;radii=[1] if args.control else [16,20.25]
    started=time.process_time();src=source(L,R)
    vertex_error=vertex_control(src)
    if not args.control:del src['H']
    records=[ball(src,r) for r in radii]
    out=dict(status='ACTUAL_SOURCE_LAPSE_ENTROPY_COMPONENT',L=L,R=R,source_difference_order=4,
             lapse_coupling='H_N={N,H}/2; N=exp(n)',exact_vertex_error=vertex_error,balls=records,
             selected_lapse=False,source_clock_identified_with_electron=False,global_stationarity=False,
             cpu_seconds=time.process_time()-started)
    if args.control:out['dense_control']=dense_control(src,records[0])
    name='CONTROL.json' if args.control else 'RESULT.json'
    (HERE/name).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
    if vertex_error>1e-10 or (args.control and out['dense_control']['last_error']>1e-6):raise SystemExit(1)

if __name__=='__main__':main()
