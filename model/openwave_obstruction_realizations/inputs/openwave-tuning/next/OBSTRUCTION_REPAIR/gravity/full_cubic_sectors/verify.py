import importlib.util,itertools,json,time
from pathlib import Path
import numpy as np
from group import full_small_bases

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('continuum_flat_source',HERE.parent/'spectral_flattening/source.py')
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
start=time.process_time();L=21;h=.5;R=1.;mass=1e-4
indices=np.stack(np.meshgrid(*([np.arange(L)]*3),indexing='ij'),-1)
axis=2*np.pi*np.fft.fftfreq(L,d=h)
p=np.stack(np.meshgrid(*([axis]*3),indexing='ij'),-1)
H,P,E=s.symbol(p,R,mass)
cells=np.array([c for c in itertools.product(range(-3,4),repeat=3) if sum(x*x for x in c)<=9])
N=len(cells);delta=cells[:,None,:]-cells[None,:,:];ix=delta%L

def regional(momentum):
    kernel=np.fft.ifftn(momentum,axes=(0,1,2))
    return kernel[ix[:,:,0],ix[:,:,1],ix[:,:,2]].transpose(0,2,1,3).reshape(4*N,4*N)

C=regional(P)
epsilon=.002
_,Pchanged,_=s.symbol(np.exp(-epsilon)*p,R,mass)
Cchanged=regional(Pchanged)
# A genuine cubic lapse variation in the actual flattened clock.
deltaC=np.zeros_like(C)
for direction in range(3):
    for sign in [-1,1]:
        shift=np.zeros(3,int);shift[direction]=sign
        iq=(indices+shift)%L
        Pq=P[iq[...,0],iq[...,1],iq[...,2]]
        Eq=E[iq[...,0],iq[...,1],iq[...,2]]
        derivative=(Eq-E)[...,None,None]*(Pq-P)/(2*(Eq+E)[...,None,None])
        kernel=np.fft.ifftn(derivative,axes=(0,1,2))
        phase=np.exp(2j*np.pi*(cells@shift)/L)
        block=kernel[ix[:,:,0],ix[:,:,1],ix[:,:,2]]*phase[:,None,None,None]
        deltaC+=block.transpose(0,2,1,3).reshape(4*N,4*N)/6

bases,orbits=full_small_bases(cells)
assembled=np.hstack([b for group in bases.values() for b in group])
gram=assembled.conj().T@assembled

def first(C,dC):
    e,u=np.linalg.eigh(C);e=np.clip(e,1e-13,1-1e-13)
    return float(np.sum(np.log((1-e)/e)*np.diag(u.conj().T@dC@u).real))

rows=[];entropy=0.;entropy_changed=0.;lapse_response=0.;equal_error=0.;off_error=0.
for name,components in bases.items():
    B=components[0];block=B.conj().T@C@B;changed=B.conj().T@Cchanged@B;dblock=B.conj().T@deltaC@B
    d=len(components)
    for i,Bi in enumerate(components):
        equal_error=max(equal_error,float(np.max(abs(Bi.conj().T@C@Bi-block))),float(np.max(abs(Bi.conj().T@deltaC@Bi-dblock))))
        for j,Bj in enumerate(components):
            if i!=j:off_error=max(off_error,float(np.max(abs(Bi.conj().T@C@Bj))))
    entropy+=d*s.entropy(block);entropy_changed+=d*s.entropy(changed);lapse_response+=d*first(block,dblock)
    rows.append(dict(irrep=name,dimension=d,multiplicity=B.shape[1],entropy=d*s.entropy(block),uniform_metric_entropy_change=d*(s.entropy(changed)-s.entropy(block)),lapse_entropy_derivative=d*first(block,dblock)))

out=dict(status='ALL_CUBIC_SPINOR_SECTORS_CONTROL',L=L,h=h,R=R,mass=mass,mass_scope='Small positive gap removes the isolated finite torus zero-mode ambiguity; no species mass fitted.',cells=N,full_modes=4*N,represented_modes=sum(r['dimension']*r['multiplicity'] for r in rows),sectors=rows,complete_basis_error=float(np.linalg.norm(assembled@assembled.conj().T-np.eye(4*N))),gram_error=float(np.linalg.norm(gram-np.eye(4*N))),equal_component_error=equal_error,off_component_error=off_error,max_orbit_complete_error=max(o['complete_error'] for o in orbits),lapse_hermiticity=float(np.linalg.norm(deltaC-deltaC.conj().T)),full_entropy=s.entropy(C),reduced_entropy=entropy,entropy_error=abs(s.entropy(C)-entropy),full_metric_entropy_change=s.entropy(Cchanged)-s.entropy(C),reduced_metric_entropy_change=entropy_changed-entropy,metric_entropy_change_error=abs((s.entropy(Cchanged)-s.entropy(C))-(entropy_changed-entropy)),full_lapse_derivative=first(C,deltaC),reduced_lapse_derivative=lapse_response,lapse_derivative_error=abs(first(C,deltaC)-lapse_response),cpu_seconds=time.process_time()-start)
out['pass']=out['represented_modes']==4*N and max(out[k] for k in ['complete_basis_error','gram_error','equal_component_error','off_component_error','lapse_hermiticity','entropy_error','metric_entropy_change_error','lapse_derivative_error'])<1e-8
with open('CONTROL.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
assert out['pass']
