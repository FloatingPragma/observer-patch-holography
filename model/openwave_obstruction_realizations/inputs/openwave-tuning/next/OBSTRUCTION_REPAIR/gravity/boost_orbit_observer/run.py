import argparse,importlib.util,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,svdvals

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('flattening_source',HERE.parent/'spectral_flattening/source.py')
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
parser=argparse.ArgumentParser()
parser.add_argument('--count',type=int,default=256)
parser.add_argument('--box',type=float,default=64.)
parser.add_argument('--output',default='RESULT.json')
args=parser.parse_args()
start=time.process_time()
N=args.count;box=args.box;h=box/N
x=(np.arange(N)-N/2)*h
p=2*np.pi*(np.arange(N)-N/2+.5)/box
F=np.exp(1j*x[:,None]*p[None,:])/np.sqrt(N)
Xp=F.conj().T@(x[:,None]*F)
Kraw=.5*(p[:,None]+p[None,:])*Xp
same_sign=np.sign(p[:,None])*np.sign(p[None,:])>0
K=Kraw*same_sign
K=(K+K.conj().T)/2
e,Q=eigh(K,check_finite=False)
q=np.zeros((N,3));q[:,0]=p
V=s.rotation(q)
Vd=V.conj().transpose(0,2,1)
_,P,E=s.symbol(q)
Pa=(s.I4+s.ALPHA[0])/2
Pb=s.I4-Pa
cache={}

def matrix(z):
    if z not in cache:
        cache[z]=(Q*np.exp(-1j*z*e)[None,:])@Q.conj().T
    return cache[z]

def spin(A,f):
    return np.einsum('pij,pjm->pim',A,f,optimize=True)

def boost(f,z,flat=True):
    g=spin(Vd,f) if flat else f
    gp=np.einsum('ij,pjm->pim',Pa,g,optimize=True)
    gm=np.einsum('ij,pjm->pim',Pb,g,optimize=True)
    shape=g.shape
    g=(matrix(z)@gp.reshape(N,-1)+matrix(-z)@gm.reshape(N,-1)).reshape(shape)
    return spin(V,g) if flat else g

def to_position(f):
    return (F@f.reshape(N,-1)).reshape(f.shape)

def embed(indices):
    return (F.conj().T[:,indices,None,None]*s.I4[None,None,:,:]).transpose(0,2,1,3).reshape(N,4,4*len(indices))

def entropy_of_frame(f):
    C=f.reshape(4*N,-1).conj().T@spin(P,f).reshape(4*N,-1)
    return s.entropy(C)

random=np.random.default_rng(1109)
test=random.normal(size=(N,4,3))+1j*random.normal(size=(N,4,3))
test/=np.linalg.norm(test,axis=(0,1),keepdims=True)
a=boost(boost(test,.2),-.07)
b=boost(test,.13)
unitarity=float(np.linalg.norm(np.einsum('pim,pin->mn',test.conj(),test)-np.einsum('pim,pin->mn',b.conj(),b)))
vacuum_comm=float(np.linalg.norm(boost(spin(P,test),.3)-spin(P,boost(test,.3))))
group_error=float(np.linalg.norm(a-b))

wave=[]
ae,avec=np.linalg.eigh(s.ALPHA[0])
for z in [.2,.5]:
    for p0,width in [(1.,.3),(3.,.5)]:
        for col in [0,3]:
            alpha=float(ae[col]);v=avec[:,col]
            f=np.exp(-.5*((p-p0)/width)**2)[:,None]*v[None,:]
            fnorm=np.linalg.norm(f)
            initial=spin(V,(f/fnorm)[:,:,None])
            actual=boost(initial,z)
            expected=spin(V,(np.exp(alpha*z/2)*np.exp(-.5*((np.exp(alpha*z)*p-p0)/width)**2)[:,None]*v[None,:]/fnorm)[:,:,None])
            wave.append(dict(rapidity=z,p0=p0,width=width,alpha=alpha,analytic_wavepacket_error=float(np.linalg.norm(actual-expected)),analytic_norm=float(np.linalg.norm(expected))))

rows=[]
for ell in [2.,4.,8.]:
    region=np.flatnonzero(abs(x)<=ell)
    f=embed(region)
    S0=entropy_of_frame(f)
    for z in [.3,.6]:
        transformed=boost(f,z)
        transformed_dirac=boost(f,z,flat=False)
        S1=entropy_of_frame(transformed)
        gram=transformed.reshape(N*4,-1).conj().T@transformed.reshape(N*4,-1)
        overlaps=[]
        for margin in sorted(set([1.,2.,4.,ell])):
            d=ell+np.exp(abs(z))*ell+margin
            if d+np.exp(abs(z))*ell>=box/2:
                overlaps.append(dict(margin=margin,centre_distance=d,excluded='Geometric interval reaches the finite box boundary; no locality value claimed.'))
                continue
            phase=np.exp(-1j*p*d)[:,None,None]
            a=to_position(phase*transformed)[region].reshape(len(region)*4,-1)
            b=to_position(phase*transformed_dirac)[region].reshape(a.shape)
            overlaps.append(dict(margin=margin,centre_distance=d,flat_overlap_operator_norm=float(svdvals(a,check_finite=False)[0]),dirac_overlap_operator_norm=float(svdvals(b,check_finite=False)[0]),flat_overlap_mean_square=float(np.linalg.norm(a)**2/a.shape[1])))
        row=dict(ell=ell,rapidity=z,speed=float(np.tanh(z)),regional_modes=f.shape[-1],entropy_unboosted=S0,entropy_boosted=S1,entropy_difference=S1-S0,frame_gram_error=float(np.linalg.norm(gram-np.eye(gram.shape[0]))),overlaps=overlaps)
        rows.append(row)
        print(json.dumps(dict(ell=ell,rapidity=z,entropy_difference=S1-S0,cpu_seconds=time.process_time()-start)),flush=True)

out=dict(status='COMPLETE_FINITE_BOOST_ORBIT_APPARATUS',configuration=vars(args),spacing=h,retained_modes=4*N,finite_generator_closure='Pinch {x,p}/2 into positive/negative momentum blocks; no mode removed. Complete self-adjoint finite boost, with continuum geometric accuracy assessed separately.',raw_to_pinched_relative_frobenius=float(np.linalg.norm(K-Kraw)/np.linalg.norm(Kraw)),unitarity_error=unitarity,composition_error=group_error,vacuum_commutator_error=vacuum_comm,wavepacket_controls=wave,regions=rows,cpu_seconds=time.process_time()-start,full_3d_geometric_locality=False,stationarity_or_physical_G=False)
with open(args.output,'w') as f:json.dump(out,f,indent=2)
assert unitarity<1e-10 and group_error<1e-10 and vacuum_comm<1e-10
