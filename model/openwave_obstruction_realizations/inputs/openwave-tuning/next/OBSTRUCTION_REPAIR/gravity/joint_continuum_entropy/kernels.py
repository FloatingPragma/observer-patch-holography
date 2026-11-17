import importlib.util,itertools
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
def sibling(name,path):
    spec=importlib.util.spec_from_file_location(name,HERE.parent/path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

cell=sibling('continuous_cell_source','continuum_cell_regions/source.py')
s=cell.s

def vertices(p,q,R,mass):
    b=p+q
    A=np.einsum('...j,jab->...ab',b,s.ALPHA)
    B=np.einsum('...j,jab->...ab',p,s.ALPHA)
    r2=np.sum(p*p,axis=-1);rq2=np.sum(b*b,axis=-1)
    M=mass+R**3*r2*r2;Mq=mass+R**3*rq2*rq2
    E=np.sqrt(r2+M*M);Eq=np.sqrt(rq2+Mq*Mq)
    T=(B+M[...,None,None]*s.BETA)/E[...,None,None]
    Tq=(A+Mq[...,None,None]*s.BETA)/Eq[...,None,None]
    deltaD=-(A+B)/2
    deltaQ=-.5*((r2+rq2)[...,None,None]*s.I4+2*A@B)
    V=deltaD+R**3*(r2+rq2)[...,None,None]*(s.BETA@deltaQ)
    metric=-.5*(V-Tq@V@T)/(Eq+E)[...,None,None]
    Ef=np.sqrt(r2+mass*mass);Eqf=np.sqrt(rq2+mass*mass)
    lapse=-(Eqf-Ef)[...,None,None]*(Tq-T)/(4*(Eqf+Ef)[...,None,None])
    return metric,lapse

def response_kernels(L,h,R,mass,aliases,twist,shift):
    axis=2*np.pi*(np.fft.fftfreq(L)*L+twist)/(L*h)
    p0=np.stack(np.meshgrid(*([axis]*3),indexing='ij'),-1)
    q=2*np.pi*np.asarray(shift)/(L*h)
    metric=np.zeros((L,L,L,4,4),complex);lapse=np.zeros_like(metric)
    for alias in itertools.product(range(-aliases,aliases+1),repeat=3):
        p=p0+2*np.pi*np.array(alias)/h
        w,_=cell.window(p,h);wq,_=cell.window(p+q,h)
        dm,dn=vertices(p,q,R,mass)
        metric+=(w*wq)[...,None,None]*dm
        lapse+=(w*wq)[...,None,None]*dn
    return np.fft.ifftn(metric,axes=(0,1,2)),np.fft.ifftn(lapse,axes=(0,1,2))

def regional_response(kernel,source,cells,shift):
    d=cells[:,None,:]-cells[None,:,:];ix=d%source['L']
    phase=np.exp(2j*np.pi*source['twist']*d.sum(-1)/source['L'])
    phase=phase*np.exp(2j*np.pi*(cells@np.asarray(shift))/source['L'])[:,None]
    blocks=kernel[ix[...,0],ix[...,1],ix[...,2]]*phase[...,None,None]
    full=blocks.transpose(0,2,1,3).reshape(4*len(cells),4*len(cells))
    return (full+full.conj().T)/2
