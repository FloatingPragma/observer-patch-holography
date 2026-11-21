"""Complete finite counts, the dimer clock and a frequency reader blind to source levels."""
import numpy as np
from scipy.linalg import hankel
from scipy.optimize import least_squares
from scipy.stats import binom

from model.candidate_laws.scalar_channel.carrier import require


def counts(rng, law, shots):
    """Digitize the cumulative law, then sample an exhaustive binary tree of partitions."""
    p = np.asarray(law,float)
    require(p.ndim >= 1 and np.isfinite(p).all() and p.min() >= -2e-12
            and np.max(abs(p.sum(axis=-1)-1)) < 1e-10, 'incomplete pointer law')
    size = p.shape[-1]
    rows = np.maximum(p.reshape(-1,size),0.)
    edges = np.rint(np.cumsum(rows/rows.sum(axis=1)[:,None],axis=1)*2**30).astype(np.int64)
    require(np.all(edges[:,-1]==2**30),'invalid digitized pointer')
    width = 1 << (size-1).bit_length()
    masses = np.zeros((len(rows),width),np.int64)
    masses[:,:size] = np.diff(edges,prepend=0,axis=1)
    totals = np.full((len(rows),1),shots,np.int64)
    while totals.shape[1] < width:
        halves = masses.reshape(len(rows),2*totals.shape[1],-1).sum(axis=2)
        denominator = halves[:,::2]+halves[:,1::2]
        probability = halves[:,::2]/np.maximum(denominator,1)
        q = (rng.integers(2**52,size=totals.shape)+.5)/2**52
        left = binom.ppf(q,totals,probability).astype(np.int64)
        totals = np.stack([left,totals-left],axis=-1).reshape(len(rows),-1)
    out = totals[:,:size].reshape(p.shape).astype(np.int32)
    require(np.all(out>=0) and np.all(out.sum(axis=-1)==shots),'counter lost a trial')
    return out


def reference(rng, times, shots):
    values,U = np.linalg.eigh([[1.,-1.],[-1.,1.]])
    wave = (np.exp(-1j*np.outer(times,values))*U[0]) @ U.T
    cos = abs(wave[:,0])**2
    sin = abs((wave[:,0]-1j*wave[:,1])/np.sqrt(2))**2
    p = np.column_stack([cos,sin])
    return counts(rng,np.stack([1-p,p],axis=-1),shots)


def read_clock(data, shots):
    require(data.ndim==3 and data.shape[1:]==(2,2) and data.dtype.kind in 'iu' and
            np.all(data>=0) and np.all(data.sum(axis=-1)==shots),'invalid clock counts')
    z,y = (2.*data[:,:,1]/shots-1).T
    require(np.max(abs(np.hypot(z,y)-1)) < .02,'clock contrast lost')
    phase = np.unwrap(np.angle(z+1j*y))/2
    require(np.all((np.diff(phase)>0)&(np.diff(phase)<np.pi/4)),'unresolved clock')
    return phase-phase[0]


def path_counts(rng, amplitude, shots):
    p = (1+np.column_stack([amplitude.real,amplitude.imag]))/2
    return counts(rng,np.stack([1-p,p],axis=-1),shots)


def energies(car, raw, shots, law, home=(0,), active=None):
    """Decode bond occupations and contact readings into energies in dimer units."""
    i,j = np.array(car.seams).T
    edges = len(i)
    parked = law == 'parked_neutral'
    size = car.n if parked else edges
    require(raw.shape[-2:]==(2*edges+size+car.n,3) and raw.dtype.kind in 'iu'
            and np.all(raw>=0) and np.all(raw.sum(axis=-1)==shots),'invalid local energy record')
    # Only the neutral oscillator can return occupation two. The charged bond modes
    # are fermionic and the contact setting measures a projector.
    require(not np.any(raw[..., :2*edges, 2]) and not np.any(raw[..., -car.n:, 2]),
            'double occupation in a binary energy detector')
    mean = (raw[...,1]+2.*raw[...,2])/shots
    charged = 2*mean[...,:2*edges].sum(axis=-1)
    field = mean[...,2*edges:2*edges+size]*(car.degree if parked else 2)
    if parked:
        remote = ~np.isin(np.arange(car.n)//12,home)
    else:
        remote = (i//12==j//12) & ~np.isin(i//12,home)
    if active is None:
        active = np.full(car.n,law!='dephased_contact')
    contact = mean[...,-car.n:] @ np.where(active,2,1)
    neutral = field.sum(axis=-1)
    return dict(charged_energy=charged.tolist(),neutral_energy=neutral.tolist(),
                remote_neutral_energy=field[...,remote].sum(axis=-1).tolist(),
                contact_energy=contact.tolist(),total_energy=(charged+neutral+contact).tolist())


def lines(data, ticks, spec):
    shots = spec['shots']
    require(data.shape==(len(ticks),2,2) and data.dtype.kind in 'iu' and np.all(data>=0)
            and np.all(data.sum(axis=-1)==shots),'invalid interferometer counts')
    z = 2.*data[:,:,1]/shots-1
    record = z[:,0]+1j*z[:,1]
    k = spec['fit_samples']
    y = record[:k]
    t = ticks[:k]
    order = 32
    block = hankel(y[:order],y[order-1:])
    U,s,Vh = np.linalg.svd(block[:,:-1],full_matrices=False)
    rank = int(np.sum(s > .018*s[0]))
    require(0 < rank < order,'unresolved pencil rank')
    roots = np.linalg.eigvals((U[:,:rank].conj().T @ block[:,1:] @ Vh[:rank].conj().T)/s[:rank])
    initial = np.sort(-np.angle(roots)/np.mean(np.diff(t)))
    weights = np.linalg.lstsq(np.exp(-1j*np.outer(t,initial)),y,rcond=None)[0].real

    def error(x):
        r = np.exp(-1j*np.outer(t,x[:rank])) @ x[rank:]-y
        return np.r_[r.real,r.imag]

    fit = least_squares(error,np.r_[initial,weights],max_nfev=100)
    freq, weight = fit.x[:rank],fit.x[rank:]
    tail = np.exp(-1j*np.outer(ticks[k:],freq)) @ weight-record[k:]
    require(fit.success and np.all(weight>0) and abs(sum(weight)-1)<.01,'invalid line fit')
    return dict(frequencies=freq.tolist(),weights=weight.tolist(),heldout_rms=float(np.sqrt(np.mean(abs(tail)**2))))
