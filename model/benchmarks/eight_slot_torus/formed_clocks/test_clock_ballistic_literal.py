"""Independent real-space birth, signed Fourier amplitudes and finite receiver taps."""
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from . import clock_ballistic


MU = np.array([128, 129])*math.pi/4096
SLOTS = np.arange(8)
S = np.full(8, 1/math.sqrt(8))
OFFSETS = np.arange(-12,13)
EXTENT = 16385
V = 2/7
MATURE_TIMES = 7*np.array([448,560,672,784,896,1008,1120])


def _route(wave):
    """Literal reflected routing on the (1,1,L,1) carrier."""
    out = np.empty_like(wave)
    for slot in SLOTS:
        delta = (1 if slot == 4 else -1) if slot in (4,5) else 0
        out[:,slot^1] = np.roll(wave[:,slot],delta)
    return out


def _native_birth():
    coin = 2*np.outer(S,S)-np.eye(8)
    feeds = []
    for site,slot in ((-1,4),(1,5)):
        wave = np.zeros((9,8))
        wave[site] = coin[:,slot]
        feeds.append(_route(wave @ coin.T))
    incoming = feeds[0][:,5]*feeds[1][:,4]
    capture = math.sin(math.pi/4)*incoming
    assert abs(abs(capture[0])**2-.5) < 2e-14
    assert np.max(abs(capture[1:])) < 2e-14
    return abs(capture[0])**2


def _apertures():
    # Bracket each actual turning point without using the quadratic root formula.
    saddles = []
    for g in MU:
        def velocity(k):
            E = math.acos((3+math.cos(k))*math.cos(g)/4)
            return math.cos(g)*math.sin(k)/(4*math.sin(E))
        peak = minimize_scalar(lambda k:-velocity(k),bounds=(0,math.pi),method="bounded").x
        saddles.append([brentq(lambda k:velocity(k)-V,*interval,xtol=1e-14)
                        for interval in ((0,peak),(peak,math.pi))])
    saddles = np.array(saddles)
    result = []
    for target in (0,1):
        suppressed = np.r_[saddles[:,1-target],-saddles.ravel()]
        # Write H(z)=sum conjugate(q_j)*z**(j+12). The required notches
        # make H a multiple of the six-root polynomial. Its shifted
        # coefficient columns span every physically allowed tap vector.
        polynomial = np.polynomial.polynomial.polyfromroots(np.exp(1j*suppressed))
        basis = np.zeros((25,19),complex)
        for j in range(19):
            basis[j:j+7,j] = polynomial.conj()
        wanted = np.exp(1j*OFFSETS*saddles[:,target].mean())
        q = basis @ np.linalg.lstsq(basis,wanted,rcond=None)[0]
        q /= np.linalg.norm(q)
        np.testing.assert_allclose(q.conj() @ np.exp(1j*np.outer(OFFSETS,suppressed)),0,
                                   rtol=0,atol=1e-11)
        result.append(q)
    delta = np.eye(25)[12]
    return result+[delta,result[0]], saddles


def _born_fields():
    # No wrap can occur before age 129 on this 263-site chart. Embed the
    # two output fields, with their actual relative phase, in the long ring.
    sample = []
    for g in MU:
        wave = np.zeros((263,8),complex)
        wave[0] = S
        rail = []
        for age in range(130):
            if age >= 128:
                large = np.zeros(EXTENT,complex)
                sites = np.arange(-129,130)
                large[sites%EXTENT] = wave[sites%263] @ S
                rail.append(np.fft.fft(large))
            if age < 129:
                bright = wave.mean(axis=1,keepdims=True)
                coin = np.exp(-1j*g)*bright-np.exp(1j*g)*(wave-bright)
                theta = math.pi*(1-math.cos(math.pi*(age+1)/128))/4 if age<128 else math.pi/2
                phase = np.exp(-2j*theta)
                wave = ((1+phase)*coin+(1-phase)*_route(coin))/2
        sample.append(rail)
    return np.array(sample)


def test_every_ballistic_original_start_word():
    q, roots = _apertures()
    for r,g in enumerate(MU):
        np.testing.assert_allclose(clock_ballistic.saddles(g),roots[r],atol=3e-13,rtol=0)
    for i,name in enumerate(("low","high")):
        np.testing.assert_allclose(clock_ballistic.aperture(name),q[i],atol=2e-10,rtol=0)
    birth = _native_birth()
    raw = _born_fields()
    k = 2*np.pi*np.fft.fftfreq(EXTENT)
    energy = np.arccos(np.cos(MU)[:,None]*(3+np.cos(k))/4)
    expected = np.empty((4,len(MATURE_TIMES),6))
    for j,n in enumerate(MATURE_TIMES):
        # Cayley-Hamilton on the bright two-dimensional plane, applied to
        # the projected fields from the two consecutive *literal* steps.
        transform = (np.sin(n*energy)*raw[:,1]-np.sin((n-1)*energy)*raw[:,0])/np.sin(energy)
        z = 2*n//7
        local = np.fft.ifft(transform,axis=1)[:,z+OFFSETS]
        for i,ap in enumerate(q):
            a,b = math.sqrt(birth/2)*(local @ ap.conj())
            p = abs(a+np.array([1,1j,-1,-1j])*b)**2/4
            if i == 3:
                p = np.full(4,(abs(a)**2+abs(b)**2)/4)
            expected[i,j] = [*p,birth-p.sum(),1-birth]
    np.testing.assert_allclose(clock_ballistic.source(),expected,rtol=0,atol=2e-12)
