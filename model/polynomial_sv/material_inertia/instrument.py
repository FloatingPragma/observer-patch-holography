"""Clock calibration, complete two-port counts and a line reader without Hamiltonian access."""
import math

import numpy as np


def counts(probability, shots):
    """Count all N uniformly spaced midpoint thresholds, including complementary outcomes."""
    p = np.asarray(probability)
    if not np.isfinite(p).all() or np.any(p < -1e-10) or np.any(p > 1+1e-10):
        raise ValueError('invalid pointer probability')
    bright = np.floor(np.clip(p,0,1)*shots+.5).astype(np.int32)
    return np.stack((shots-bright,bright), axis=-1)


def quadratures(amplitude, shots):
    return counts((1+np.stack([amplitude.real, amplitude.imag], axis=-1))/2, shots)


def signal(readings, length, shots):
    if (readings.shape != (length,2,2) or readings.dtype != np.int32 or
            np.any(readings < 0) or not np.all(readings.sum(axis=-1) == shots)):
        raise ValueError('incomplete interferometer menu')
    q = 2*readings[:,:,1]/shots - 1
    return q[:,0]+1j*q[:,1]


def clock(spec):
    ticks = np.array(spec['clock_ticks'])
    # Occupation and a balanced read in the one-load sector of an isolated seam.
    from scipy.linalg import expm
    H = .5*np.array([[1.,-1.],[-1.,1.]])
    arm = np.array([expm(-1j*t*H/spec['tick_divisor'])[:,0] for t in ticks])
    p = np.stack([abs(arm[:,0])**2, abs((arm[:,0]-1j*arm[:,1])/np.sqrt(2))**2], axis=-1)
    return counts(p,spec['shots'])


def clock_rate(data, spec):
    ticks = np.array(spec['clock_ticks'])
    z = signal(data,len(ticks),spec['shots'])
    angle = np.unwrap(np.angle(z))
    if (np.min(abs(z)) < .99 or not np.all(np.diff(angle)>0) or
            np.max(np.diff(angle)) >= np.pi/2):
        raise ValueError('invalid or stopped material clock')
    rate = float(ticks @ angle / (ticks @ ticks))
    if np.max(abs(rate*ticks-angle)) > .001:
        raise ValueError('unstable material clock')
    return rate


def tick_menu(spec):
    dyadic = 2**np.arange(spec['dyadic_orders'], dtype=np.int64)
    return np.r_[dyadic, 3*dyadic[-1]//4]


def read_line(data, rate, spec):
    ticks = tick_menu(spec)
    z = signal(data,len(ticks),spec['shots'])
    floor = spec['dominant_floor']
    # A source audit verifies the dominant weight separately. A poor return can disqualify it here.
    if min(abs(z)) < 2*floor-1-math.sqrt(2)/spec['shots']:
        return None
    angles = -np.angle(z)
    phase = angles[0]
    for wrapped in angles[1:-1]:
        expected = 2*phase
        phase = wrapped + 2*np.pi*round((expected-wrapped)/(2*np.pi))
        if abs(phase-expected) > spec['phase_gate']:
            return None
    frequency = phase/(ticks[-2]*rate)
    heldout = np.angle(np.exp(1j*(angles[-1]-frequency*ticks[-1]*rate)))
    if abs(heldout) > spec['phase_gate']:
        return None
    angular_bound = math.asin((1-floor+math.sqrt(2)/spec['shots'])/floor)
    return float(frequency), float(angular_bound/(ticks[-2]*rate))


def quotient(numerator, denominator):
    """Positive observable ratio with worst-case intervals, including an unresolved denominator."""
    a,da = numerator
    b,db = denominator
    if b <= db or a <= da:
        return None
    return [a/b, (a-da)/(b+db), (a+da)/(b-db)]


def band_ratios(lines):
    if any(v is None for v in lines):
        return None
    one, moving_one, cluster, moving_cluster = lines
    def shift(a,b):
        return b[0]-a[0], a[1]+b[1]
    return {'rest':quotient(cluster,one), 'inertia':quotient(shift(one,moving_one),shift(cluster,moving_cluster))}
