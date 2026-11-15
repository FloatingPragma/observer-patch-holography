#!/usr/bin/env python3
"""n3dlib: the declared cubic routing of twelve-port carriers and the candidate H_m = D + m Gamma on it.

Own construction from the definition in
  reverse-engineering-reality/model/native_spatial_refinement/PROOF.md  (section "Construction and complete records"):
ports (0, +-1, +-phi) and cyclic permutations in the order of `raw_frame`; internal seams = pairs at squared
distance four; in each coordinate-zero group the two ports with positive product of the other two coordinates are
joined to the port of the positive neighbouring cell obtained by flipping the cyclic last nonzero coordinate.
Nothing here imports or writes into the repository.  numpy / scipy only.

Index conventions
  cell (x, y, z) of a torus (Lx, Ly, Lz)      ->  c = (x * Ly + y) * Lz + z
  port i of cell c                            ->  12 c + i
  seam s of cell c (s = 0..29 internal, 30..35 external; external seam 30 + 2 a + j is link j of axis a,
                    from port p of cell c to port q of cell c + e_a)        ->  36 c + s
  position (ports, then transit seams)        ->  ports 0..P-1, transit P..P+M-1
  component of a cell (Bloch index)           ->  0..11 ports, 12..47 transit seams
Seam orientation: first endpoint (amplitude +1 of D) = lower port index for internal seams, the port p of cell c for
an external seam.  The orientation is a gauge (ASSUMPTIONS.md D2).
"""
from __future__ import annotations

import itertools
import math

import numpy as np
import scipy.sparse as sp

PHI = (1 + 5 ** 0.5) / 2
C3SQ = 1.0 / 8.0                       # proved second-order coefficient of the lowest port band (NSR)

# ---------------------------------------------------------------------------------------------------------
# exact arithmetic in Z[phi]: (a, b) = a + b phi, phi^2 = phi + 1
# ---------------------------------------------------------------------------------------------------------
Z0, Z1, ZP = (0, 0), (1, 0), (0, 1)


def zadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def zneg(x):
    return (-x[0], -x[1])


def zmul(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)


def zsign(x):
    """sign of a + b phi; a + b phi = 0 only for a = b = 0 (phi irrational)"""
    if x == (0, 0):
        return 0
    return 1 if x[0] + x[1] * PHI > 0 else -1


def frame():
    """The twelve ports in the order of the repository's raw_frame, exact in Z[phi]."""
    m1, mp = zneg(Z1), zneg(ZP)
    return ((m1, ZP, Z0), (Z1, ZP, Z0), (m1, mp, Z0), (Z1, mp, Z0),
            (Z0, m1, ZP), (Z0, Z1, ZP), (Z0, m1, mp), (Z0, Z1, mp),
            (ZP, Z0, m1), (ZP, Z0, Z1), (mp, Z0, m1), (mp, Z0, Z1))


def frame_float():
    return np.array([[c[0] + c[1] * PHI for c in p] for p in frame()])


def internal_seams():
    """Thirty pairs (i < j) at squared distance exactly four."""
    v = frame()
    out = []
    for i in range(12):
        for j in range(i + 1, 12):
            d = Z0
            for a, b in zip(v[i], v[j]):
                diff = zadd(a, zneg(b))
                d = zadd(d, zmul(diff, diff))
            if d == (4, 0):
                out.append((i, j))
    assert len(out) == 30
    return tuple(out)


def links(control=False):
    """External links per axis: ((p, q), (p', q')), port p of cell c joined to port q of cell c + e_axis.

    control=True: the repository's alternative axis-zero partition (for axis 0 the positive ports are those with
    positive coordinate one, and that coordinate is flipped across the link)."""
    v = frame()
    out = []
    for axis in range(3):
        row = []
        for i, p in enumerate(v):
            if p[axis] != Z0:
                continue
            positive = zsign(zmul(p[(axis + 1) % 3], p[(axis + 2) % 3])) > 0
            flip = (axis + 2) % 3
            if control and axis == 0:
                positive = zsign(p[1]) > 0
                flip = 1
            if positive:
                q = list(p)
                q[flip] = zneg(q[flip])
                row.append((i, v.index(tuple(q))))
        assert len(row) == 2
        out.append(tuple(row))
    assert sorted(x for row in out for e in row for x in e) == list(range(12))
    return tuple(out)


INTERNAL = internal_seams()


def cell_seams(control=False):
    """36 seams of a cell: (first port, second port, axis) with axis = -1 for internal seams; the second port of
    an external seam lies in the cell c + e_axis."""
    out = [(a, b, -1) for a, b in INTERNAL]
    for axis, row in enumerate(links(control)):
        for p, q in row:
            out.append((p, q, axis))
    return tuple(out)


# ---------------------------------------------------------------------------------------------------------
# the torus graph
# ---------------------------------------------------------------------------------------------------------
class Torus:
    def __init__(self, shape, control=False):
        if isinstance(shape, int):
            shape = (shape,) * 3
        self.shape = tuple(int(s) for s in shape)
        assert all(s >= 2 for s in self.shape), "periodic size at least two per axis"
        self.control = control
        lx, ly, lz = self.shape
        self.cells = lx * ly * lz
        self.P = 12 * self.cells
        self.M = 36 * self.cells
        cs = cell_seams(control)
        coords = np.array(list(itertools.product(range(lx), range(ly), range(lz))))
        self.coords = coords
        cid = lambda c: (c[:, 0] % lx * ly + c[:, 1] % ly) * lz + c[:, 2] % lz
        first = np.zeros((self.cells, 36), dtype=np.int64)
        second = np.zeros((self.cells, 36), dtype=np.int64)
        base = 12 * np.arange(self.cells)
        for s, (a, b, axis) in enumerate(cs):
            first[:, s] = base + a
            if axis < 0:
                second[:, s] = base + b
            else:
                nxt = coords.copy()
                nxt[:, axis] += 1
                second[:, s] = 12 * cid(nxt) + b
        self.first = first.reshape(-1)          # seam index 36 c + s
        self.second = second.reshape(-1)
        self.seam_axis = np.tile(np.array([c[2] for c in cs]), self.cells)

    @property
    def seams(self):
        return np.stack([self.first, self.second], axis=1)

    def incidence(self, phases=None):
        """P x M sparse incidence B (amplitude +1 at the first endpoint, -1 at the second).
        phases: optional (M, 2) array theta of half-hop phases, B[port, e] = sign * exp(+i theta[e, side])."""
        m = self.M
        rows = np.concatenate([self.first, self.second])
        cols = np.concatenate([np.arange(m), np.arange(m)])
        if phases is None:
            vals = np.concatenate([np.ones(m), -np.ones(m)])
        else:
            vals = np.concatenate([np.exp(1j * phases[:, 0]), -np.exp(1j * phases[:, 1])])
        return sp.csr_matrix((vals, (rows, cols)), shape=(self.P, m))

    def hamiltonian(self, mass=0.0, phases=None):
        """H_m = [[m, B], [B^dagger, -m]] on ports (+) transit seams, sparse."""
        b = self.incidence(phases)
        return sp.bmat([[mass * sp.identity(self.P), b], [b.conj().T, -mass * sp.identity(self.M)]], format="csr")

    def laplacian(self):
        b = self.incidence()
        return (b @ b.T).tocsr()

    def connected(self):
        from scipy.sparse.csgraph import connected_components
        b = self.incidence()
        adj = abs(b) @ abs(b).T
        n, _ = connected_components(adj, directed=False)
        return n == 1

    def degrees(self):
        return np.bincount(np.concatenate([self.first, self.second]), minlength=self.P)

    def position_of_components(self):
        """Declared position apparatus A (see SOLUTION.md): offsets of the 48 components of a cell."""
        return apparatus_A(self.control)


def apparatus_A(control=False):
    """Every port and every internal transit seam of cell n at the point n; the transit position of an external
    seam of axis a at n + e_a / 2 (midpoint of the two cell centres).  Returns 48 x 3 offsets."""
    off = np.zeros((48, 3))
    for s, (_, _, axis) in enumerate(cell_seams(control)):
        if axis >= 0:
            off[12 + s, axis] = 0.5
    return off


def apparatus_B(control=False, radius=0.25):
    """Alternative apparatus: port i at n + radius * (unit icosahedron vertex i); a transit seam at the midpoint of
    its two ports (for an external seam the second port lies in the next cell)."""
    f = frame_float()
    f = radius * f / np.linalg.norm(f[0])
    off = np.zeros((48, 3))
    off[:12] = f
    for s, (a, b, axis) in enumerate(cell_seams(control)):
        shift = np.zeros(3)
        if axis >= 0:
            shift[axis] = 1.0
        off[12 + s] = 0.5 * (f[a] + f[b] + shift)
    return off


# ---------------------------------------------------------------------------------------------------------
# Bloch objects.  Convention: amplitude on component i of cell n is s_i exp(i k.n)
# ---------------------------------------------------------------------------------------------------------
def bloch_incidence(k, control=False):
    """B(k): (..., 12, 36).  Internal seam (a, b): +1, -1.  External seam (p, q, axis): B[p] = 1,
    B[q] = -exp(-i k_axis), so that (B^dagger u)_e = u_p - exp(i k_axis) u_q."""
    k = np.asarray(k, float)
    out = np.zeros(k.shape[:-1] + (12, 36), complex)
    for s, (a, b, axis) in enumerate(cell_seams(control)):
        out[..., a, s] = 1.0
        if axis < 0:
            out[..., b, s] = -1.0
        else:
            out[..., b, s] = -np.exp(-1j * k[..., axis])
    return out


def port_symbol(k, control=False):
    """H(k) = 6 I - A_ico - sum_a sum_{(p,q) in E_a} (e^{i k_a} |p><q| + e^{-i k_a} |q><p|), shape (..., 12, 12)."""
    k = np.asarray(k, float)
    out = np.zeros(k.shape[:-1] + (12, 12), complex)
    out[..., np.arange(12), np.arange(12)] = 6.0
    for a, b in INTERNAL:
        out[..., a, b] -= 1.0
        out[..., b, a] -= 1.0
    for axis, row in enumerate(links(control)):
        ph = np.exp(1j * k[..., axis])
        for p, q in row:
            out[..., p, q] -= ph
            out[..., q, p] -= np.conj(ph)
    return out


def hm_symbol(k, mass, control=False):
    b = bloch_incidence(k, control)
    out = np.zeros(b.shape[:-2] + (48, 48), complex)
    out[..., :12, 12:] = b
    out[..., 12:, :12] = np.conj(np.swapaxes(b, -1, -2))
    idx = np.arange(48)
    out[..., idx, idx] = np.where(idx < 12, mass, -mass)
    return out


def lowest_band(k, control=False):
    """Lowest port band: (lambda0, u0) with the gauge sum_i u0_i real positive (overlap with the uniform vector)."""
    w, v = np.linalg.eigh(port_symbol(k, control))
    lam = w[..., 0]
    u = v[..., :, 0]
    s = u.sum(-1)
    u = u * (np.conj(s) / np.abs(s))[..., None]
    return lam, u


def band_velocity(k, u, control=False):
    """Hellmann-Feynman gradient of lambda0: d lambda0 / d k_a = 2 sum_{(p,q) in E_a} Im(e^{i k_a} conj(u_p) u_q)."""
    k = np.asarray(k, float)
    out = np.zeros(k.shape)
    for axis, row in enumerate(links(control)):
        ph = np.exp(1j * k[..., axis])
        for p, q in row:
            out[..., axis] += 2 * np.imag(ph * np.conj(u[..., p]) * u[..., q])
    return out


def positive_spinor(k, mass, control=False):
    """Positive-energy spinor of the lowest port band: (E, lambda0, u0, s) with s = (P, Q) of 48 components,
    P = sqrt((E + m) / 2E) u0, Q = B^dagger u0 / sqrt(2E (E + m))."""
    lam, u = lowest_band(k, control)
    lam = np.maximum(lam, 0.0)
    e = np.sqrt(mass * mass + lam)
    b = bloch_incidence(k, control)
    bu = np.einsum("...pe,...p->...e", np.conj(b), u)
    s = np.zeros(lam.shape + (48,), complex)
    with np.errstate(divide="ignore", invalid="ignore"):
        s[..., :12] = np.sqrt((e + mass) / (2 * e))[..., None] * u
        s[..., 12:] = bu / np.sqrt(2 * e * (e + mass))[..., None]
    return e, lam, u, s


# ---------------------------------------------------------------------------------------------------------
# packets on the torus (momentum representation, exact for the torus)
# ---------------------------------------------------------------------------------------------------------
def kgrid(shape):
    ks = [2 * np.pi * np.fft.fftfreq(n) for n in shape]
    return np.stack(np.meshgrid(*ks, indexing="ij"), axis=-1)


class Packet:
    """Single-branch positive-energy packet of the lowest port band.

    Momentum amplitude f(k) = g(k) exp(-|k - k0|^2 / (4 sigma^2)), g(k) = 1 - exp(-|k|^2 / (2 kappa^2)) (vanishes
    at k = 0; kappa=None gives g = 1), k - k0 reduced to (-pi, pi].  The state is
    Psi(n, i) = N^{-1/2} sum_k f(k) s_i(k) exp(i k.(n - n0)), normalised to one."""

    def __init__(self, shape, k0, sigma, mass, n0, kappa=None, control=False, cut=1e-13):
        if isinstance(shape, int):
            shape = (shape,) * 3
        self.shape = tuple(shape)
        self.control = control
        self.mass = float(mass)
        self.k0 = np.asarray(k0, float)
        self.sigma = float(sigma)
        self.kappa = kappa
        self.n0 = np.asarray(n0, float)
        kg = kgrid(self.shape)
        dk = (kg - self.k0 + np.pi) % (2 * np.pi) - np.pi
        f = np.exp(-(dk ** 2).sum(-1) / (4 * sigma ** 2))
        if kappa is not None:
            f = f * (1 - np.exp(-(kg ** 2).sum(-1) / (2 * kappa ** 2)))
        mask = f > cut * f.max()
        if self.mass == 0.0:
            mask &= (kg ** 2).sum(-1) > 0          # the positive branch does not exist at k = 0 for m = 0
        self.mask = mask
        self.k = kg[mask]                          # (n_k, 3)
        f = f[mask]
        self.f = f / np.sqrt((f ** 2).sum())
        self.energy, self.lam, self.u, self.s = positive_spinor(self.k, self.mass, control)
        self.w = self.f ** 2                       # momentum weights, sum to one
        self.phase0 = np.exp(-1j * (self.k @ self.n0))

    # ---- predictions from the prepared momentum distribution (no fitted parameter) ----
    def velocity_lattice(self):
        return band_velocity(self.k, self.u, self.control) / (2 * self.energy)[:, None]

    def velocity_continuum(self, coeff=None):
        c = np.full(3, C3SQ) if coeff is None else np.asarray(coeff, float)
        e = np.sqrt(self.mass ** 2 + (c * self.k ** 2).sum(-1))
        return c * self.k / e[:, None]

    def moments_of(self, v):
        mean = (self.w[:, None] * v).sum(0)
        second = np.einsum("k,ka,kb->ab", self.w, v, v)
        return mean, second - np.outer(mean, mean)

    def central_ray(self, coeff=None):
        c = np.full(3, C3SQ) if coeff is None else np.asarray(coeff, float)
        return c * self.k0 / math.sqrt(self.mass ** 2 + (c * self.k0 ** 2).sum())

    # ---- states ----
    def khat(self, t=0.0):
        """Momentum amplitudes at time t, shape (n_k, 48)."""
        return (self.f * self.phase0 * np.exp(-1j * self.energy * t))[:, None] * self.s

    def realspace(self, t=0.0, amplitudes=None):
        """Psi(n, i) at time t on the torus, shape (Lx, Ly, Lz, 48)."""
        amp = self.khat(t) if amplitudes is None else amplitudes
        full = np.zeros(self.shape + (48,), complex)
        full[self.mask] = amp
        n = np.prod(self.shape)
        return np.fft.ifftn(full, axes=(0, 1, 2)) * math.sqrt(n)


def measure(psi, offsets, shell=2, weights_only=False):
    """Centroid, covariance tensor, total probability and edge-shell probability of a state psi(n, i) under the
    position apparatus `offsets` (48 x 3): component i of cell n sits at n + offsets[i].  Coordinates are the
    arithmetic ones in [0, L) (valid while the edge-shell probability is negligible)."""
    shape = psi.shape[:3]
    nc = psi.shape[-1]                                     # 48, or 12 when only the ports are passed
    p = np.abs(psi) ** 2                                   # (Lx, Ly, Lz, nc)
    total = p.sum()
    axes = [np.arange(n, dtype=float) for n in shape]
    pc = p.sum(-1)                                         # cell marginal
    comp = p.reshape(-1, nc).sum(0)                        # component populations
    # first moments
    m1 = np.array([(pc.sum(axis=tuple(j for j in range(3) if j != a)) * axes[a]).sum() for a in range(3)])
    cen = (m1 + comp @ offsets) / total
    # second moments: sum p(n,i) (n + o_i)_a (n + o_i)_b
    second = np.zeros((3, 3))
    for a in range(3):
        for b in range(a, 3):
            if a == b:
                nn = (pc.sum(axis=tuple(j for j in range(3) if j != a)) * axes[a] ** 2).sum()
            else:
                marg = pc.sum(axis=3 - a - b)              # remaining two axes, in increasing order (a < b)
                nn = (marg * np.outer(axes[a], axes[b])).sum()
            second[a, b] = second[b, a] = nn
    # cross terms n_a o_b: need sum_n n_a p(n, i) for each component
    na = np.zeros((3, nc))
    for a in range(3):
        marg = p.sum(axis=tuple(j for j in range(3) if j != a))     # (L_a, nc)
        na[a] = axes[a] @ marg
    cross = na @ offsets                                    # (3, 3): sum n_a o_b p
    second = second + cross + cross.T + (offsets.T * comp) @ offsets
    cov = second / total - np.outer(cen, cen)
    edge = np.ones(shape, bool)
    edge[shell:shape[0] - shell, shell:shape[1] - shell, shell:shape[2] - shell] = False
    return {"centroid": cen, "cov": cov, "norm": float(total), "edge": float(pc[edge].sum() / total),
            "port_weight": float(comp[:12].sum() / total)}


# ---------------------------------------------------------------------------------------------------------
# real-space propagation (Chebyshev expansion of exp(-i H t)); used as the second method
# ---------------------------------------------------------------------------------------------------------
def chebyshev_propagate(h, psi, t, bound):
    """exp(-i h t) psi for a Hermitian sparse h with spectrum inside [-bound, bound]."""
    from scipy.special import jv
    a = bound * t
    order = int(abs(a) + 12 * abs(a) ** (1 / 3) + 30)
    hs = h / bound
    t0 = psi.astype(complex)
    t1 = hs @ t0
    out = jv(0, a) * t0 + 2 * (-1j) * jv(1, a) * t1
    for n in range(2, order + 1):
        t0, t1 = t1, 2 * (hs @ t1) - t0
        out += 2 * (-1j) ** n * jv(n, a) * t1
    return out


# ---------------------------------------------------------------------------------------------------------
# native commit dissipator on bare positions (zero-holonomy sector), exact element-wise form
# ---------------------------------------------------------------------------------------------------------
def kraus_bare(n_ports, seams, seam):
    """The ten local-coarse C/E Kraus operators of one seam on bare positions (ports, then transit seams):
    list of (label, source indices, target indices); every operator is 1/2 times a partial permutation
    (fair coin, equal C/E mixture).  A transit position is 'elsewhere' for every seam."""
    n = n_ports + len(seams)
    a, b = seams[seam]
    idx = np.arange(n)
    not_a, not_b = idx[idx != a], idx[idx != b]
    els = idx[(idx != a) & (idx != b)]
    A, Bp = np.array([a]), np.array([b])
    return [(("C", 0, 1), A, Bp), (("C", 0, 0), not_a, not_a), (("C", 1, -1), Bp, A), (("C", 1, 0), not_b, not_b),
            (("E", 0, 1), A, Bp), (("E", 0, -1), Bp, Bp), (("E", 0, 0), els, els),
            (("E", 1, 1), A, A), (("E", 1, -1), Bp, A), (("E", 1, 0), els, els)]
