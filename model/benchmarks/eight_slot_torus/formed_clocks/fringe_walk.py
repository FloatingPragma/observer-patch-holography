"""Marked and coherently supplied pairs under the periodic native opening."""
from __future__ import annotations

import numpy as np

from .. import formation, law
from .clock_comparison import optical_counter

MEMBER = law.LAWS["torus"]
MU = MEMBER.rail_contact_phases[0] / 2
PERIOD = formation.AGE_PERIOD
OPEN = MEMBER.maturation[2]


def require(test, reason):
    if not test:
        raise ValueError(reason)


def appointment(size, window):
    """A local two-service optical return with a supplied nonwrapping counter."""
    if window == "ballistic":
        target, denominator = size, 1
    elif window in ("diffusive", "diffusive_double"):
        target, denominator = size**2, 32 if window == "diffusive" else 16
    else:
        raise ValueError("unknown appointment")
    returns = (target + denominator) // (2 * denominator)
    # Small cubes include a birth reading before the first return.
    return 2 * returns


def optical_reference():
    """Check the native two-cycle once. Larger counts use a wider classical register."""
    returns, powers = optical_counter((2,1,1,2), (0,0,0,0), (1,2))
    require(returns == {1:2, 2:4} and np.allclose(powers, (0,1), atol=2e-13, rtol=0),
            "reference does not trigger on even returns")


def angles(shape):
    axes = [2*np.pi*np.fft.fftfreq(n) for n in shape]
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)


class BrightWalk:
    """Orthonormal basis S and (FS-aS)/sqrt(1-a²); the second vector is unused at rest."""
    def __init__(self, wavevectors, mu=MU):
        k = np.asarray(wavevectors, dtype=float)
        require(k.shape[-1] == 3 and np.all(np.isfinite(k)), "invalid wave vectors")
        self.k = k
        # Use sin(k/2) to retain the small difference 1-a at long wavelength.
        deficit = (np.sin(k/2)**2).sum(axis=-1)/2
        self.a = 1-deficit
        self.b = np.sqrt(deficit*(2-deficit))
        self.routing = np.empty((*self.a.shape,2,2))
        self.routing[...,0,0] = self.a
        self.routing[...,0,1] = self.routing[...,1,0] = self.b
        self.routing[...,1,1] = -self.a
        self.coin = np.array([np.exp(-1j*mu), -np.exp(1j*mu)])
        self.mature = self.routing*self.coin
        self.identity = np.broadcast_to(np.eye(2), self.mature.shape)
        self.prefix = {0:self.identity}

    def opening(self, ages):
        if ages in self.prefix:
            return self.prefix[ages]
        start = max(j for j in self.prefix if j < ages)
        matrix = self.prefix[start]
        for j in range(start, ages):
            theta = law.hop_angle(j, MEMBER.maturation)
            stay = np.exp(-1j*theta)*np.cos(theta)
            hop = 1j*np.exp(-1j*theta)*np.sin(theta)
            matrix = ((stay*self.identity+hop*self.routing)*self.coin) @ matrix
        self.prefix[ages] = matrix
        return matrix

    def block(self, n):
        if n < OPEN:
            return self.opening(n)
        return np.linalg.matrix_power(self.mature, n-OPEN) @ self.opening(OPEN)

    def wave(self, n):
        require(isinstance(n, int) and n >= 0, "invalid number of services")
        periods, tail = divmod(n, PERIOD)
        start = self.block(tail)
        if periods:
            start = start @ np.linalg.matrix_power(self.block(PERIOD), periods)
        u, v = start[...,0,0], start[...,1,0]
        # F on the equal slot ray carries the incoming partner's translation.
        phases = np.zeros((*self.k.shape[:-1],8))
        phases[...,0:6:2] = self.k
        phases[...,1:6:2] = -self.k
        basis = np.zeros_like(phases, dtype=complex)
        np.divide(np.expm1(1j*phases)+(1-self.a)[...,None], self.b[...,None],
                  where=self.b[...,None] != 0, out=basis)
        wave = (u[...,None]+v[...,None]*basis)*law.EQUAL
        require(np.max(abs((abs(wave)**2).sum(axis=-1)-1))<3e-8, "lost bright-state weight")
        return wave


def contrasts(shape, order, boosts, windows):
    """One local detector. The full grid is used only to predict its marked-birth clicks."""
    grid = angles(shape)
    walk = BrightWalk(grid)
    grating = sum(2*np.pi*m*x/n for x,m,n in zip(np.indices(shape), order, shape, strict=True))
    analyzer = np.exp(-1j*grating)
    result = {}
    for window in windows:
        ticks = appointment(max(shape),window)
        wave = walk.wave(ticks)
        position = np.fft.ifftn(wave, axes=(0,1,2))
        density = (abs(position)**2).sum(axis=-1)
        require(abs(density.sum()-1)<3e-8,"incomplete spatial distribution")
        marked = np.sum(density*analyzer)
        for boost in boosts:
            low = tuple(boost*m % n for m,n in zip(order,shape,strict=True))
            high = tuple((boost+1)*m % n for m,n in zip(order,shape,strict=True))
            result[window,boost] = marked, np.vdot(wave[low], wave[high])
    return result


def jet_product(left, right):
    """Taylor coefficients through k² for two matrix functions."""
    return np.stack([sum(left[i]@right[j-i] for i in range(j+1)) for j in range(3)])


def jet_power(matrix, n):
    out = np.zeros_like(matrix)
    out[0] = np.eye(2)
    while n:
        if n%2:
            out = jet_product(out,matrix)
        matrix = jet_product(matrix,matrix)
        n //= 2
    return out


def cycle_jet(mu):
    """Full-period coefficients through order two in an axial wave number."""
    jet = np.array([[[1.,0.],[0.,-1.]], [[0.,.5],[.5,0.]], [[-.125,0.],[0.,.125]]])
    coin = np.array([np.exp(-1j*mu), -np.exp(1j*mu)])
    opening = jet_power(jet*coin,0)
    for age in range(OPEN):
        theta = law.hop_angle(age,MEMBER.maturation)
        gate = 1j*np.exp(-1j*theta)*np.sin(theta)*jet
        gate[0] += np.exp(-1j*theta)*np.cos(theta)*np.eye(2)
        opening = jet_product(gate*coin,opening)
    return jet_product(jet_power(jet*coin,PERIOD-OPEN), opening)


def envelope_coefficient():
    """Derivative of the isolated bright eigenphase of one *complete* age period."""
    cycle = cycle_jet(MU)
    even, other = np.diag(cycle[0])
    require(abs(even-1)<1e-10 and abs(other+1j)<1e-10, "unexpected age monodromy")
    second = cycle[2,0,0]+cycle[1,0,1]*cycle[1,1,0]/(even-other)
    require(abs((second/even).real)<1e-7, "nonunitary eigenphase derivative")
    return float(-(second/even).imag/PERIOD)


def finite_envelope(k, ticks, curvature):
    """Bound the vector error of the scalar parabolic forecast at one finite mode."""
    w=BrightWalk(k)
    periods,tail=divmod(ticks,PERIOD)
    cycle=w.block(PERIOD)
    eig,vec=np.linalg.eig(cycle)
    rest=np.exp(-1j*MU*PERIOD)
    line=int(np.argmin(abs(eig-rest)))
    v=vec[:,line]
    projected_leak=np.linalg.norm(v[1:])/np.linalg.norm(v)
    # Phase-distance bounds telescope for unitary powers. Prefix error includes
    # its age-dependent slot rotation. The cap of two holds for unit vectors.
    freq=curvature*np.dot(k,k)
    local=w.block(tail)[:,0]-np.array([np.exp(-1j*(MU+freq)*tail),0.])
    bound=(np.linalg.norm(local)+2*projected_leak+
           periods*abs(eig[line]-rest*np.exp(-1j*PERIOD*freq)))
    return min(2.,float(bound))+3e-8
