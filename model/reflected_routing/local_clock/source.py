"""Lift the transverse cut to a three-bit tag. The receiver occupies one chart address."""
import numpy as np

from ..kernel import multipliers, packet, require, step


def launch(size,center,tags):
    modes,_,small,_ = packet(size,center,2)
    periodic = np.zeros((len(modes[0]),size,8),complex)
    periodic[:,modes[1]%size] = small
    periodic = np.fft.ifft(periodic,axis=1)*np.sqrt(size)
    extent = size*tags
    lifted = np.zeros((len(modes[0]),extent,8),complex)
    y = np.arange(-size//2,size//2)
    lifted[:,y%extent] = periodic[:,y%size]
    source = np.fft.fft(lifted,axis=1)/np.sqrt(extent)
    angles = np.zeros((*source.shape[:2],4))
    angles[...,0] = 2*np.pi*modes[0][:,None]/size
    angles[...,1] = 2*np.pi*np.arange(extent)/extent
    require(abs(np.vdot(source,source)-1)<1e-11,'winding lift lost initial weight')
    return modes[0],angles,source


def exact_powers(initial,angles):
    shift = multipliers(angles)
    equal = np.ones(8)/np.sqrt(8)
    alpha = np.cos(angles).mean(axis=-1)
    energy = np.arccos(alpha)
    nontrivial = abs(alpha)<1
    w = np.zeros_like(shift)
    w[nontrivial] = ((shift[...,np.arange(8)^1]*equal-alpha[...,None]*equal)[nontrivial]
                     /np.sin(energy[nontrivial,None]))
    a = initial@equal
    b = (w.conj()*initial).sum(axis=-1)
    remaining = initial-a[...,None]*equal-b[...,None]*w
    odd = -(shift*remaining)[...,np.arange(8)^1]
    def power(t):
        c,s = np.cos(t*energy),np.sin(t*energy)
        return ((a*c-b*s)[...,None]*equal+(a*s+b*c)[...,None]*w
                +(odd if t%2 else remaining))
    return power


def aperture(state,xmodes,size,x,tags):
    spatial = np.einsum('i,ijk->jk',np.exp(2j*np.pi*x*xmodes/size),state)
    # At y=0, tag w has Fourier weight exp(2 pi i j w/tags).
    classes = spatial.reshape(size,tags,8).sum(axis=0)
    amplitudes = np.fft.ifft(classes,axis=0)*np.sqrt(tags)/size
    # z and the private address are also fixed at the receiver, each costing 1/L.
    return (abs(amplitudes)**2).sum(axis=1)/size**2


def observe(size,center,x,times,spec,analytic=False):
    tags = spec['tag_modulus']
    require(size>0 and size%2==0 and tags in (1,8) and 0<=x<size and len(times)>0 and
            all(isinstance(t,(int,np.integer)) and 0<=t<6*size for t in times) and
            all(b>a for a,b in zip(times,times[1:])),'invalid aperture or observation schedule')
    shots = spec['shots_per_four_volume']*size**4
    modes,angles,initial = launch(size,center,tags)
    shift = multipliers(angles)
    power = exact_powers(initial,angles) if analytic else None
    state = initial
    previous = 0
    records = []
    for tick in times:
        if analytic:
            state = power(tick)
        else:
            for _ in range(tick-previous):
                state = step(state,shift)
        previous = tick
        p = aperture(state,modes,size,x,tags)
        require(np.all(np.isfinite(p)) and p.sum()<=1 and tick<6*size,'out of winding domain')
        edges = np.floor(shots*p.cumsum()+0.5).astype(np.int64)
        records.append(np.diff(np.r_[np.int64(0),edges,np.int64(shots)]))
    return np.asarray(records)
