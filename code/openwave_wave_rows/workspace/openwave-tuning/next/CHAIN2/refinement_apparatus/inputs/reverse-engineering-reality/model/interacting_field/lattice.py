"""Local two-component CAR field with all spatial modes retained.

The supplied open cube has two fermion modes at each site. H=dGamma(h)+U D,
where h(k)=-c sum sigma_a sin(k_a) and D=sum n_up n_down. Execution is on the
EXACT N=2 invariant sector, not a Gaussian covariance or occupation truncation.
"""
from __future__ import annotations
from itertools import product
from math import isfinite
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

SIGMA = (np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex),
         np.diag([1., -1.]).astype(complex))


def real(value, name, *, nonnegative=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float)):
        raise ValueError(f"finite real {name} required")
    if not isfinite(value) or (nonnegative and value < 0):
        raise ValueError(f"finite nonnegative {name} required")
    return float(value)


class TwoFermions:
    """Exterior-square basis and arbitrary number-preserving local CAR gates."""
    def __init__(self, modes):
        if type(modes) is not int or not 2 <= modes <= 700:
            raise ValueError("2..700 modes required")
        self.modes = modes
        self.first, self.second = np.triu_indices(modes, 1)
        self.dimension = len(self.first)
        self._indices = {}

    def index(self, p, q):
        if type(p) is not int or type(q) is not int or not 0 <= p < q < self.modes:
            raise ValueError("two ordered distinct modes required")
        return p*(2*self.modes-p-1)//2 + q-p-1

    def basis_state(self, p, q):
        out = np.zeros(self.dimension, complex)
        out[self.index(p, q)] = 1
        return out

    def state(self, value):
        if np.ma.isMaskedArray(value):
            raise ValueError("masked state forbidden")
        a = np.asarray(value)
        if a.ndim not in (1, 2) or a.shape[0] != self.dimension or a.dtype.kind not in 'ficu' or not np.isfinite(a).all():
            raise ValueError("finite state columns on the complete N=2 sector required")
        return np.asarray(a, complex)

    def pair_indices(self, p, q):
        self.index(p, q)
        if (p, q) not in self._indices:
            r = np.array([r for r in range(self.modes) if r != p and r != q])
            def one(a):
                lo, hi = np.minimum(a, r), np.maximum(a, r)
                idx = lo*(2*self.modes-lo-1)//2+hi-lo-1
                return idx, np.where(a < r, 1., -1.)
            ip, sp = one(p); iq, sq = one(q)
            self._indices[p, q] = (ip, iq, sp*sq, self.index(p, q))
        return self._indices[p, q]

    def apply_pair_inplace(self, state, p, q, gate):
        """Lift all four local Fock inputs, including vacuum and double occupancy.

        Gate parity/number conservation is a premise checked on its original
        matrix. CAR orientation signs implement the exact routed local map.
        This does not use Wick factorization of the evolving non-Gaussian state.
        """
        gate = np.asarray(gate, complex)
        if gate.shape != (4, 4) or not np.isfinite(gate).all():
            raise ValueError("finite complete two-mode gate required")
        forbidden = gate.copy(); forbidden[0, 0] = forbidden[3, 3] = 0
        forbidden[1:3, 1:3] = 0
        if np.linalg.norm(forbidden, 2) > 2e-10:
            raise ValueError("gate must preserve every local number sector")
        ip, iq, signs, both = self.pair_indices(p, q)
        a, b = state[ip].copy(), state[iq].copy()
        sign = signs if state.ndim == 1 else signs[:, None]
        occupied = state[both].copy()
        # Vacuum is retained, including its actual source-executed scalar.
        if gate[0, 0] != 1:
            state *= gate[0, 0]
        state[ip] = gate[1, 1]*a + sign*gate[1, 2]*b
        state[iq] = sign*gate[2, 1]*a + gate[2, 2]*b
        state[both] = gate[3, 3]*occupied

    def apply_pair(self, state, p, q, gate):
        out = self.state(state).copy()
        self.apply_pair_inplace(out, p, q, gate)
        return out

    def one_body(self, state):
        """gamma=A A*, Tr gamma=2 ||psi||^2; no state normalization."""
        state = self.state(state)
        if state.ndim != 1:
            raise ValueError("one pure state required")
        a = np.zeros((self.modes, self.modes), complex)
        a[self.first, self.second] = state
        a[self.second, self.first] = -state
        return a@a.conj().T

    def lift_one_body(self, h):
        h = sparse.coo_matrix(h)
        if h.shape != (self.modes, self.modes):
            raise ValueError("one-body shape mismatch")
        rows, cols, vals = [], [], []
        for p, q, z in zip(h.row, h.col, h.data, strict=True):
            if p == q:
                idx = np.flatnonzero((self.first == p) | (self.second == p))
                rows.append(idx); cols.append(idx); vals.append(np.full(len(idx), z))
            else:
                lo, hi = sorted((int(p), int(q)))
                ip, iq, signs, _ = self.pair_indices(lo, hi)
                rows.append(ip if p == lo else iq)
                cols.append(iq if p == lo else ip)
                vals.append(z*signs)
        if not vals:
            return sparse.csr_matrix((self.dimension,self.dimension),dtype=complex)
        return sparse.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                                 shape=(self.dimension, self.dimension)).tocsr()


class WeylLattice(TwoFermions):
    def __init__(self, radius, speed=0.125, interaction=1.0):
        if type(radius) is not int or not 1 <= radius <= 3:
            raise ValueError("retained execution radius must be 1..3")
        self.radius = radius
        self.speed = real(speed, 'coordinate speed', nonnegative=True)
        self.interaction = real(interaction, 'interaction', nonnegative=True)
        self.sites = tuple(product(range(-radius, radius+1), repeat=3))
        self.site_index = {x: i for i, x in enumerate(self.sites)}
        self.volume = len(self.sites)
        super().__init__(2*self.volume)
        groups = [[] for _ in range(6)]
        rows, cols, values = [], [], []
        for i, x in enumerate(self.sites):
            for axis in range(3):
                y = list(x); y[axis] += 1
                j = self.site_index.get(tuple(y))
                if j is None:
                    continue
                for a in range(2):
                    for b in range(2):
                        z = 1j*self.speed*SIGMA[axis][a, b]/2
                        if z != 0:
                            p, q = 2*i+a, 2*j+b
                            groups[2*axis+x[axis] % 2].append((p, q, complex(z)))
                            rows.extend((p, q)); cols.extend((q, p)); values.extend((z, z.conjugate()))
        self.groups = tuple(tuple(g) for g in groups)
        self.single_h = sparse.coo_matrix((values, (rows, cols)), shape=(self.modes, self.modes)).tocsr()
        self.contact_indices = np.array([self.index(2*i, 2*i+1) for i in range(self.volume)])
        diag = np.zeros(self.dimension); diag[self.contact_indices] = self.interaction
        self.hamiltonian = self.lift_one_body(self.single_h)+sparse.diags(diag)
        self.origin = 2*self.site_index[(0, 0, 0)]
        self.neighbor = 2*self.site_index[(1, 0, 0)]

    def exact_evolve(self, state, duration):
        duration = real(duration, 'model time')
        return expm_multiply(-1j*duration*self.hamiltonian, self.state(state))


def node_table():
    return [{'pi_coordinates': list(bits), 'jacobian_chirality': -(-1)**sum(bits)}
            for bits in product((0, 1), repeat=3)]


def symbol(momentum, speed=.125):
    if not isinstance(momentum, (list, tuple)) or len(momentum) != 3:
        raise ValueError("three coordinate momenta required")
    c = real(speed, 'speed', nonnegative=True)
    return -c*sum((np.sin(real(k, 'momentum'))*s for k, s in zip(momentum, SIGMA, strict=True)),
                  np.zeros((2, 2), complex))
