"""Open neutral material chain and all finite time-filter exits."""

import numpy as np
from scipy.linalg import schur

from ...row_20_electric import experiment as prior
from ...row_20_electric.schedule import DELTA, rotate_pairs


def end_words(c):
    a, b = 0, len(c) - 2
    joint = c[a, a+1]*c[b, b+1] - c[a, b]*c[a+1, b+1] + c[a, b+1]*c[a+1, b]
    return np.array([(1 + s*c[a, a+1] + t*c[b, b+1] + s*t*joint)/4
                     for s in (1, -1) for t in (1, -1)])


def aperture(angle, size):
    """Four independent coherent uniform labels; one local operation at their summed time."""
    x = (np.asarray(angle) + np.pi) % (2*np.pi) - np.pi
    ratio = np.sinc(size*x/(2*np.pi)) / np.sinc(x/(2*np.pi))
    return (np.exp(1j * (size-1)*x/2) * ratio)**4


def leakage(delta, size):
    if not 0 < delta < np.pi:
        raise ValueError('unresolved filter gap or alias')
    return min(1., 1/(size*np.sin(delta/2)))**8


class Guide:
    def __init__(self, n, mode, source_duty=1.):
        j, h = prior.couplings(8)
        duty = (2/65) * (.5 if mode == 'half' else 1.)
        self.period = DELTA*(1+duty)
        self.n = n
        r = np.eye(2*n)
        rotate_pairs(r, np.arange(0, 2*n, 2), -2*h*DELTA)
        links = np.arange(1, 2*n-1, 2)
        if mode == 'cut':
            links = np.delete(links, n//2-1)
        rotate_pairs(r, links[links != 1], 2*j*duty*DELTA)
        if 1 in links:
            rotate_pairs(r, np.array([1]), 2*j*duty*DELTA*source_duty)
        tri, basis = schur(r, output='complex')
        np.testing.assert_allclose(tri, np.diag(tri.diagonal()), atol=2e-13, rtol=0)
        angles = np.angle(tri.diagonal())
        chosen = np.flatnonzero(angles < -1e-12)
        chosen = chosen[np.argsort(-angles[chosen])]
        if len(chosen) != n:
            raise ArithmeticError('ambiguous normal-mode vacuum')
        self.r = r
        self.phi = -angles[chosen]
        self.modes = np.sqrt(2)*basis[:, chosen]
        self.c = -(self.modes @ self.modes.conj().T).imag
        np.testing.assert_allclose(self.c @ self.c, -np.eye(2*n), atol=2e-12, rtol=0)
        self.source = self.modes[0].conj()
        self.receiver = self.modes[-1]
        self.cs = self.c.copy()
        self.cs[0, 1:] *= -1
        self.cs[1:, 0] *= -1

    def end_read(self, t, sham):
        if sham:
            return end_words(self.c)
        r = np.linalg.matrix_power(self.r, t)
        return end_words(r @ self.cs @ r.T)

    def absorb(self, energy, start, sham):
        theta = energy*self.period
        create = self.receiver.conj()*aperture(theta+self.phi, 2048)
        if sham:
            weights = abs(create)**2
            vacuum_tag = float(weights @ abs(aperture(self.phi, 4096))**2)
            admitted = float(weights.sum())
            vacuum = 0.
        else:
            incoming = self.source*np.exp(-1j*start*self.phi)
            annihilate = self.receiver*aperture(theta-self.phi, 2048)
            vacuum = float(abs(annihilate @ incoming)**2)
            pair = create[:, None]*incoming[None, :] - incoming[:, None]*create[None, :]
            weights = np.triu(abs(pair)**2, 1)
            residual = abs(aperture(self.phi[:, None]+self.phi[None, :], 4096))**2
            vacuum_tag = vacuum + float((weights*residual).sum())
            admitted = vacuum + float(weights.sum())
        return np.array([vacuum_tag, admitted-vacuum_tag, 1-admitted]), vacuum

    def bright_bound(self, target, threshold):
        cutoff = (threshold-target)*self.period
        if self.phi[-1]*2 >= np.pi:
            raise ArithmeticError('two-particle residual aliases the vacuum')
        return leakage(cutoff, 2048) + leakage(float(self.phi[:2].sum()), 4096)
