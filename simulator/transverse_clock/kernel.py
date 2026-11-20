"""Eight-port reflection and flip-flop routing in Bernhard Mueller's slot convention."""
import math

import numpy as np

from simulator import law

FREE = law.LAWS["oph"]
PHASE_CONTROL = law.LAWS["control-phase-i"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def multipliers(angles):
    return np.exp(-1j * angles[..., law.SLOT_AXIS] * law.SLOT_SIGN)


def reflect(state, phase=FREE.equal_phase):
    require(state.ndim>=1 and state.shape[-1]==8,'the routing coin has eight slots')
    return law.coin(state, phase)


def step(state, shift, phase=FREE.equal_phase):
    """Fourier transform of law.step, checked on each enrolled wave number."""
    out = shift*reflect(state,phase)
    return out[..., law.PARTNER]


def real_step(state, phase=FREE.equal_phase):
    """Literal simulator update for the enrolled free and changed-coin runs."""
    require(state.ndim == 5 and all(state.shape[:-1]), 'four nonempty address coordinates required')
    laws = {FREE.equal_phase: FREE, PHASE_CONTROL.equal_phase: PHASE_CONTROL}
    return law.step(state, laws[phase])


def packet_angles(size, modes):
    angles = np.zeros((len(modes[0]),len(modes[1]),4))
    angles[...,0] = 2*np.pi*modes[0][:,None]/size
    angles[...,1] = 2*np.pi*modes[1][None,:]/size
    return angles


def packet(size, center, width):
    offsets = np.arange(-width,width+1)
    modes = [c+offsets for c in center]
    angles = packet_angles(size, modes)
    alpha = np.cos(angles).mean(axis=-1)
    require(np.max(abs(alpha))<1, 'preparation hits a singular band endpoint')
    energy = np.arccos(alpha)
    shift = multipliers(angles)
    equal = np.ones(8)/np.sqrt(8)
    partner = (shift[...,np.arange(8)^1]*equal-alpha[...,None]*equal)/np.sin(energy)[...,None]
    weights = np.array([math.comb(2*width,j) for j in range(len(offsets))],dtype=float)
    weights /= np.linalg.norm(weights)
    initial = (equal+1j*partner)*np.outer(weights,weights)[...,None]/np.sqrt(2)
    require(abs(np.vdot(initial,initial)-1)<1e-12,'unnormalized source')
    return modes,shift,initial,energy


def marginal(state, modes, size, axis):
    # Inactive transverse Fourier coordinates are orthogonal. All eight slots are counted.
    coefficients = np.moveaxis(state,axis,0)
    grid = np.zeros((size,*coefficients.shape[1:]),complex)
    grid[modes[axis]%size] = coefficients
    return size*np.sum(abs(np.fft.ifft(grid,axis=0))**2,axis=(1,2))


def counters(probabilities, shots):
    p = np.asarray(probabilities)
    require(isinstance(shots,(int,np.integer)) and not isinstance(shots,(bool,np.bool_))
            and 0<shots<=np.iinfo(np.int32).max,'exposures must fit a positive detector count')
    require(p.ndim==1 and p.size>0 and p.dtype.kind in 'fiu' and
            np.isfinite(p).all() and p.min()>=0 and abs(p.sum()-1)<2e-10,
            'incomplete detector law')
    edges = np.floor(shots*np.cumsum(p)+0.5).astype(np.int64)
    result = np.diff(edges,prepend=0)
    require(result.min()>=0 and result.sum()==shots,'lost a terminal outcome')
    return result.astype(np.int32)


def tilt_rate(tilt):
    # Stable form of acosh(mean(cosh(theta_j))).
    theta = np.asarray(tilt,dtype=float)
    require(theta.shape==(4,) and np.isfinite(theta).all(),'four real tilts required')
    return 2*math.asinh(float(np.linalg.norm(np.sinh(theta/2)))/2)


def passive_log_bound(time, distance, tilt):
    """Log upper bound on the transfer norm between supporting half-spaces."""
    theta = np.asarray(tilt,dtype=float)
    rate = tilt_rate(theta)
    require(time>=0 and distance>=0,'negative separation or time')
    d = math.exp(float(max(abs(theta))))
    return min(0., math.log(d)+math.log1p((1+d)*time*(time+1))
               +time*rate-distance*float(np.linalg.norm(theta)))
