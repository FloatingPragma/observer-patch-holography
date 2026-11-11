"""Register ensemble: a three-dimensional Z_N link field with weight exp(kappa sum_p Re W_p), numpy heat bath.

Links n[mu, x] in 0..N-1 stand for omega^n, omega = exp(2 pi i / N). The plaquette is
p = n_mu(x) + n_nu(x + mu) - n_mu(x + nu) - n_nu(x) modulo N and Re W_p = cos(2 pi p / N).
N = 3 is the colour centre (flux classes 2 and 4), N = 2 the weak centre (flux class 3).
Links of one direction with equal parity of the transverse coordinate sum share no plaquette and are drawn together
from their exact conditional distribution. The recorded tables were produced with a numba loop over the same sets.
"""
import numpy as np

PLANES = ((0, 1), (0, 2), (1, 2))


class Register:
    def __init__(self, L, N, kappa, seed=0, start='cold'):
        self.L, self.N, self.kappa = int(L), int(N), float(kappa)
        self.rng = np.random.default_rng(seed)
        shape = (3, self.L, self.L, self.L)
        self.n = np.zeros(shape, np.int64) if start == 'cold' else self.rng.integers(0, self.N, shape)
        self.cos = np.cos(2.0 * np.pi * np.arange(self.N) / self.N)
        g = np.indices(shape[1:])
        self.masks = [[(g[(mu + 1) % 3] + g[(mu + 2) % 3]) % 2 == par for par in (0, 1)] for mu in range(3)]

    def staples(self, mu):
        n, N = self.n, self.N
        out = []
        for nu in range(3):
            if nu == mu:
                continue
            out.append((np.roll(n[nu], -1, mu) - np.roll(n[mu], -1, nu) - n[nu]) % N)
            back_nu = np.roll(n[nu], 1, nu)
            out.append((-np.roll(back_nu, -1, mu) - np.roll(n[mu], 1, nu) + back_nu) % N)
        return np.stack(out)

    def sweep(self):
        for mu in range(3):
            for par in (0, 1):
                m = self.masks[mu][par]
                s = self.staples(mu)[:, m]
                w = np.stack([np.exp(self.kappa * self.cos[(v + s) % self.N].sum(0)) for v in range(self.N)])
                c = np.cumsum(w, axis=0)
                r = self.rng.random(c.shape[1]) * c[-1]
                self.n[mu][m] = (r[None, :] >= c).sum(0).clip(0, self.N - 1)

    def plaquettes(self):
        n, N = self.n, self.N
        return np.stack([(n[mu] + np.roll(n[nu], -1, mu) - np.roll(n[mu], -1, nu) - n[nu]) % N for mu, nu in PLANES])

    def mean_plaquette(self):
        return float(self.cos[self.plaquettes()].mean())

    def flux_share(self):
        """Share of plaquettes with nontrivial centre flux."""
        return float((self.plaquettes() != 0).mean())

    def gauge_transform(self, g):
        """n_mu(x) -> n_mu(x) + g(x) - g(x + mu) modulo N."""
        for mu in range(3):
            self.n[mu] = (self.n[mu] + g - np.roll(g, -1, mu)) % self.N


def one_plaquette_average(N, kappa):
    """Exact <Re W> of a single plaquette with weight exp(kappa Re W) over Z_N."""
    c = np.cos(2.0 * np.pi * np.arange(N) / N)
    w = np.exp(kappa * c)
    return float((w * c).sum() / w.sum())


def dual_critical_coupling(N):
    """Critical coupling from the duality with the Ising (N = 2) and three-state Potts (N = 3) models.

    Imported spin-model couplings: K_c = 0.22165463 (Ising), K_c = 0.550565 (Potts, Janke and Villanova 1997).
    """
    if N == 2:
        return float(-0.5 * np.log(np.tanh(0.22165463)))
    if N == 3:
        return float(np.log(1.0 + 3.0 / (np.exp(0.550565) - 1.0)) * 2.0 / 3.0)
    raise ValueError('N = 2 or 3')
