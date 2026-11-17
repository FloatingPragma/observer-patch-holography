import importlib.util,itertools
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('flat_source_for_cells',HERE.parent/'spectral_flattening/source.py')
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
GAMMA=np.array([*s.ALPHA,s.BETA]);PB=(s.I4-s.BETA)/2

def window(p,h):
    z=p*h/2
    factors=np.sinc(z/np.pi)
    w=np.prod(factors,axis=-1)
    radial_derivative=np.zeros_like(w)
    for j in range(3):
        others=[i for i in range(3) if i!=j]
        radial_derivative+=2*factors[...,j]*(np.cos(z[...,j])-factors[...,j])*np.prod(factors[...,others]**2,axis=-1)
    return w,radial_derivative

def coefficients(p,R,mass):
    r2=np.sum(p*p,axis=-1);M=mass+R**3*r2*r2
    E=np.sqrt(r2+M*M)
    coeff=np.concatenate([-.5*p/E[...,None],(.5*r2/(E*(E+M)))[...,None]],axis=-1)
    deltaM=-4*R**3*r2*r2
    deltaE=(-r2+M*deltaM)/E
    metric=np.concatenate([.5*p*(1/E+deltaE/E**2)[...,None],(-.5*(deltaM/E-M*deltaE/E**2))[...,None]],axis=-1)
    return coeff,metric

def prepare(L,h,R,mass,aliases,twist=0.):
    axis=2*np.pi*(np.fft.fftfreq(L)*L+twist)/(L*h)
    p0=np.stack(np.meshgrid(*([axis]*3),indexing='ij'),-1)
    cov=np.zeros((L,L,L,4),complex);metric=np.zeros_like(cov);shape=np.zeros_like(cov)
    gradients=np.zeros((3,L,L,L,4),complex)
    for alias in itertools.product(range(-aliases,aliases+1),repeat=3):
        p=p0+2*np.pi*np.array(alias)/h
        w,dw2=window(p,h);c,dc=coefficients(p,R,mass)
        cov+=w[...,None]**2*c
        metric+=w[...,None]**2*dc
        shape+=(-3*w*w-dw2)[...,None]*c
        for j in range(3):gradients[j]+=1j*p[...,j,None]*w[...,None]**2*c
    fft=lambda x:np.fft.ifftn(x,axes=(0,1,2))
    return dict(L=L,h=h,R=R,mass=mass,aliases=aliases,twist=twist,cov=fft(cov),metric=fft(metric),shape=fft(shape),gradients=np.array([fft(x) for x in gradients]))

def region_matrix(source,cells,key='cov'):
    d=cells[:,None,:]-cells[None,:,:];ix=d%source['L']
    coeff=source[key][ix[...,0],ix[...,1],ix[...,2]]
    if key=='shape':
        for j in range(3):coeff=coeff-source['h']*d[...,j,None]*source['gradients'][j][ix[...,0],ix[...,1],ix[...,2]]
    coeff=coeff*np.exp(2j*np.pi*source.get('twist',0.)*d.sum(-1)/source['L'])[...,None]
    blocks=np.einsum('...a,aij->...ij',coeff,GAMMA)
    if key=='cov':blocks+=np.eye(len(cells))[:,:,None,None]*PB
    return blocks.transpose(0,2,1,3).reshape(4*len(cells),4*len(cells))

def literal_scaled_region(L,h,R,mass,aliases,cells,scale,twist=0.):
    """Actual normalized indicator cells with width/centres multiplied by scale.

    Constant beta covariance stays PB*identity exactly. The remaining
    continuum-torus Fourier sum is evaluated without using a Ward subtraction.
    """
    axis=2*np.pi*(np.fft.fftfreq(L)*L+twist)/(L*h)
    p0=np.stack(np.meshgrid(*([axis]*3),indexing='ij'),-1).reshape(-1,3)
    unique,inverse=np.unique((cells[:,None,:]-cells[None,:,:]).reshape(-1,3),axis=0,return_inverse=True)
    coeff=np.zeros((len(unique),4),complex)
    # Chunk separations so no large all-momentum/all-separation array is held.
    for alias in itertools.product(range(-aliases,aliases+1),repeat=3):
        p=p0+2*np.pi*np.array(alias)/h
        w,_=window(p,h*scale);c,_=coefficients(p,R,mass)
        weighted=scale**3*w[:,None]**2*c/L**3
        for at in range(0,len(unique),64):
            phase=np.exp(1j*(unique[at:at+64]*h*scale)@p.T)
            coeff[at:at+64]+=phase@weighted
    blocks=np.einsum('...a,aij->...ij',coeff[inverse].reshape(len(cells),len(cells),4),GAMMA)
    blocks+=np.eye(len(cells))[:,:,None,None]*PB
    return blocks.transpose(0,2,1,3).reshape(4*len(cells),4*len(cells))
