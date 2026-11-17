"""X1f stage 2 library: the four-component Dirac source and a response engine for any Bloch generator.

Source (the other team's `common_dynamics/spin/dirac_proposal.py::bloch(p, h, m)` at k = h p, mu = h m; read, not
imported): one step
    U(k) = S(k) exp(-i mu beta),   S(k) = exp(-i k_z alpha_z) exp(-i k_y alpha_y) exp(-i k_x alpha_x),
    alpha_j = sigma_j (x) X,  beta = 1 (x) Z.
State (X1c's): the modes with principal quasi-energy in (-pi, 0), i.e. the negative branch of the strictly local
Hermitian generator
    G = sin H = i (U - U^dagger)/2 = cos(mu) G_hop + sin(mu) G_rest,
    G_hop = i (S - S^dagger)/2,   G_rest = (S beta + beta S^dagger)/2.
Handles (declared): cell fields nu_h, nu_r that scale the two terms, each hop of the generator between two cells
by the mean of their factors:  G' = (1/2){1 + nu_h, cos(mu) G_hop} + (1/2){1 + nu_r, sin(mu) G_rest}.
Uniform nu_h = nu_r is an exact null direction (G' = (1 + nu) G).
Variant "without the reflection points": G + w G_gap, G_gap = beta sin k_x sin k_y sin k_z (scaled with nu_h).
"""
import itertools
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import x1f  # noqa: E402   (stage-1 library: entropy, modular, CubicBloch for the scalar control)

EPS = 1e-13
X2 = np.array([[0, 1], [1, 0]], complex)
Y2 = np.array([[0, -1j], [1j, 0]], complex)
Z2 = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)
ALPHA = (np.kron(X2, X2), np.kron(Y2, X2), np.kron(Z2, X2))        # their pauli('XX'), pauli('XY'), pauli('XZ')
BETA = np.kron(I2, Z2)                                             # their pauli('ZI')
I4 = np.eye(4, dtype=complex)


def dag(a):
    return np.conj(np.swapaxes(a, -1, -2))


def dirac_S(k):
    k = np.atleast_2d(np.asarray(k, float))
    S = np.broadcast_to(I4, (len(k), 4, 4)).copy()
    for a in range(3):
        f = np.cos(k[:, a])[:, None, None] * I4 - 1j * np.sin(k[:, a])[:, None, None] * ALPHA[a]
        S = f @ S
    return S


def dirac_U(k, mu):
    return dirac_S(k) @ (math.cos(mu) * I4 - 1j * math.sin(mu) * BETA)


def dirac_parts(k, mu, w=0.0):
    """(c G_hop [+ w G_gap], s G_rest): the two terms of the generator with their coefficients."""
    k = np.atleast_2d(np.asarray(k, float))
    S = dirac_S(k)
    hop = math.cos(mu) * 0.5j * (S - dag(S))
    if w:
        hop = hop + w * np.prod(np.sin(k), axis=1)[:, None, None] * BETA
    rest = math.sin(mu) * 0.5 * (S @ BETA + BETA @ dag(S))
    return hop, rest


def speed(mu):
    return 1.0 if mu == 0 else math.sqrt(mu / math.tan(mu))


def ball_cells(r2max):
    r = int(math.floor(math.sqrt(r2max) + 1e-9))
    return np.array([c for c in itertools.product(range(-r, r + 1), repeat=3)
                     if c[0] ** 2 + c[1] ** 2 + c[2] ** 2 <= r2max], dtype=int)


def cut_faces(cells):
    s = set(map(tuple, cells))
    n = 0
    for c in s:
        for a in range(3):
            for d in (1, -1):
                t = list(c)
                t[a] += d
                n += tuple(t) not in s
    return n


class Lat:
    """Filled negative branch of a Bloch generator on an L^3 block with a twisted momentum grid."""

    def __init__(self, L, G_of_k, twist=(0.5, 0.5, 0.5)):
        self.L = L
        self.twist = np.asarray(twist, float)
        n = np.arange(L)
        self.grid = np.stack(np.meshgrid(n, n, n, indexing="ij"), -1).reshape(-1, 3)
        self.k = 2 * np.pi * (self.grid + self.twist) / L
        self.Nk = L ** 3
        self.G = G_of_k(self.k)
        self.nb = self.G.shape[-1]
        self.E, self.U = np.linalg.eigh(self.G)
        self.occ = self.E < 0
        Pk = (self.U * self.occ[:, None, :]) @ dag(self.U)
        self.Preal = np.fft.ifftn(Pk.reshape(L, L, L, self.nb, self.nb), axes=(0, 1, 2))

    def idx_shift(self, shift):
        g = (self.grid + np.asarray(shift)) % self.L
        return g[:, 0] * self.L * self.L + g[:, 1] * self.L + g[:, 2]

    def P_delta(self, D):
        """P(x - y) for unwrapped integer differences D (m, 3)."""
        D = np.asarray(D)
        i = D % self.L
        ph = np.exp(2j * np.pi * (D @ self.twist) / self.L)
        return ph[:, None, None] * self.Preal[i[:, 0], i[:, 1], i[:, 2]]

    def region_C(self, cells):
        nc, nb = len(cells), self.nb
        D = (cells[:, None, :] - cells[None, :, :]).reshape(-1, 3)
        C = self.P_delta(D).reshape(nc, nc, nb, nb).transpose(0, 2, 1, 3).reshape(nc * nb, nc * nb)
        return 0.5 * (C + C.conj().T)

    def modular(self, cells, eps=EPS):
        C = self.region_C(cells)
        S, k, c = x1f.modular(C, eps)
        nc, nb = len(cells), self.nb
        T = np.real(np.einsum("ij,ji->i", k, np.eye(len(C)) - 2 * C)).reshape(nc, nb)
        return S, k.reshape(nc, nb, nc, nb), c, T


def kappa(L, twist, cells, kpad4, shift=(0, 0, 0)):
    """kappa(k)[b, a] = sum_{j,i} exp(-i k.x_j) kpad4[j, b, i, a] exp(+i (k + q).x_i) for every k of the twisted grid,
    by accumulation on the lattice of differences and one FFT (needs |x_j - x_i| < L/2 in every component)."""
    cells = np.asarray(cells)
    nc, nb = kpad4.shape[0], kpad4.shape[1]
    q = 2 * np.pi * np.asarray(shift) / L
    D = (cells[:, None, :] - cells[None, :, :])                    # x_j - x_i
    assert np.abs(D).max() < L / 2, "region too large for the block"
    ph = np.exp(-2j * np.pi * (D @ np.asarray(twist)) / L) * np.exp(1j * (cells @ q))[None, :]
    vals = (kpad4.transpose(0, 2, 1, 3) * ph[:, :, None, None]).reshape(nc * nc, nb, nb)
    Di = (D % L).reshape(-1, 3)
    M = np.zeros((L, L, L, nb, nb), dtype=complex)
    np.add.at(M, (Di[:, 0], Di[:, 1], Di[:, 2]), vals)
    return np.fft.fftn(M, axes=(0, 1, 2)).reshape(L ** 3, nb, nb)


def response(L, twist, E, U, occ, idx2, cells, kpad4, shift, vertices, chunk=4096, masks=None):
    """First-order change of S(region at the origin) under dG with <k + q| dG |k> = V(k) for each V in vertices
    (arrays (Nk, nb, nb)); q = 2 pi shift / L, idx2 = index of k + q.  Returns one complex number per vertex.
    masks: optional list of boolean arrays over k; then returns (len(vertices), len(masks)) partial sums over the
    momenta kept by each mask (a decomposition of the same response by the momentum of the virtual transition)."""
    Nk = len(E)
    kap = kappa(L, twist, cells, kpad4, shift)
    if masks is not None:
        A = dag(U) @ kap @ U[idx2]
        f1, f2 = occ.astype(float), occ[idx2].astype(float)
        df = f2[:, :, None] - f1[:, None, :]
        dE = E[idx2][:, :, None] - E[:, None, :]
        F = np.zeros_like(dE)
        m = df != 0
        F[m] = df[m] / dE[m]
        Mk = U @ (A * np.swapaxes(F, 1, 2)) @ dag(U[idx2])
        per_k = np.array([np.einsum("kab,kba->k", Mk, V) for V in vertices])
        return np.array([[per_k[i][mk].sum() for mk in masks] for i in range(len(vertices))]) / Nk
    out = np.zeros(len(vertices), dtype=complex)
    for s0 in range(0, Nk, chunk):
        sl = slice(s0, min(Nk, s0 + chunk))
        U1, E1, f1 = U[sl], E[sl], occ[sl].astype(float)
        U2, E2, f2 = U[idx2[sl]], E[idx2[sl]], occ[idx2[sl]].astype(float)
        A = dag(U1) @ kap[sl] @ U2                                 # rows at k, columns at k + q
        df = f2[:, :, None] - f1[:, None, :]
        dE = E2[:, :, None] - E1[:, None, :]
        F = np.zeros_like(dE)
        m = df != 0
        F[m] = df[m] / dE[m]
        Mk = U1 @ (A * np.swapaxes(F, 1, 2)) @ dag(U2)
        for i, V in enumerate(vertices):
            out[i] += np.einsum("kab,kba->", Mk, V[sl])
    return out / Nk


def susceptibility(E, U, occ, idx2, Va, Vb):
    """chi_ab(q) = -(1/Nk) sum conj(W_a) W_b F: static second-order energy per cell, E2 = -(1/4) chi for cos(q.x)."""
    U2, E2, f2 = U[idx2], E[idx2], occ[idx2].astype(float)
    f1 = occ.astype(float)
    df = f2[:, :, None] - f1[:, None, :]
    dE = E2[:, :, None] - E[:, None, :]
    F = np.zeros_like(dE)
    m = df != 0
    F[m] = df[m] / dE[m]
    Wa = dag(U2) @ Va @ U
    Wb = dag(U2) @ Vb @ U
    return -complex((np.conj(Wa) * Wb * F).sum()) / len(E)


class DiracSource:
    def __init__(self, L, mu, w=0.0, twist=(0.5, 0.5, 0.5)):
        self.L, self.mu, self.w = L, mu, w
        self.lat = Lat(L, lambda k: sum(dirac_parts(k, mu, w)), twist)
        self.hop, self.rest = dirac_parts(self.lat.k, mu, w)
        self.gap_min = float(np.abs(self.lat.E).min())

    def vertices(self, shift):
        """V_h, V_r for nu = exp(i q.x): (1/2)[G_part(k) + G_part(k + q)]."""
        i2 = self.lat.idx_shift(shift)
        return i2, 0.5 * (self.hop + self.hop[i2]), 0.5 * (self.rest + self.rest[i2])

    def respond(self, cells, kpad4, shift):
        """(Z_h, Z_r): response of S(ball at the origin) to nu_h, nu_r = exp(i q.x)."""
        lat = self.lat
        i2, Vh, Vr = self.vertices(shift)
        return response(lat.L, lat.twist, lat.E, lat.U, lat.occ, i2, cells, kpad4, shift, [Vh, Vr])

    def tau(self, cells, kpad4):
        """tau_n(k) = <n, k| sum_x k_x |n, k> and the Bloch matrix of the translation sum."""
        lat = self.lat
        kap0 = kappa(lat.L, lat.twist, cells, kpad4)
        return np.real(np.einsum("kbn,kba,kan->kn", np.conj(lat.U), kap0, lat.U)), kap0


# ---- planar reduction: a slab of layers in a ring, transverse momentum as a parameter ---------------------
def slab_kernel(G_parts_of_k, kperp, Nx, h, tw=0.5, eps=EPS):
    """One transverse momentum.  Region: layers x = 0 .. h-1 of a ring of Nx layers (twist tw along x).
    Returns S(kperp) and, for each part of the generator, the kernel K_part[x'] = dS/d(nu_part on layer x')."""
    n = np.arange(Nx)
    kx = 2 * np.pi * (n + tw) / Nx
    ks = np.stack([kx, np.full(Nx, kperp[0]), np.full(Nx, kperp[1])], axis=1)
    parts = G_parts_of_k(ks)
    G = sum(parts)
    nb = G.shape[-1]
    E, U = np.linalg.eigh(G)
    occ = E < 0
    Pk = (U * occ[:, None, :]) @ dag(U)
    d = np.arange(-(h - 1), h)
    ph = np.exp(1j * np.outer(d, kx)) / Nx
    Pd = (ph @ Pk.reshape(Nx, -1)).reshape(-1, nb, nb)             # P(x - y) for x - y = d
    idx = (np.arange(h)[:, None] - np.arange(h)[None, :]) + (h - 1)
    C = Pd[idx].transpose(0, 2, 1, 3).reshape(h * nb, h * nb)
    C = 0.5 * (C + C.conj().T)
    S, k, c = x1f.modular(C, eps)
    k4 = k.reshape(h, nb, h, nb)
    f1 = occ.astype(float)
    Z = np.zeros((len(parts), Nx), dtype=complex)
    xs = np.arange(h)
    for m in range(Nx):
        q = 2 * np.pi * m / Nx
        i2 = (n + m) % Nx
        # kappa(k)[b, a] = sum_{j,i} exp(-i k x_j) k4[j, b, i, a] exp(i (k + q) x_i)
        e1 = np.exp(-1j * np.outer(kx, xs))
        e2 = np.exp(1j * np.outer(kx + q, xs))
        kap = np.einsum("kj,jbia,ki->kba", e1, k4, e2)
        A = dag(U) @ kap @ U[i2]
        df = f1[i2][:, :, None] - f1[:, None, :]
        dE = E[i2][:, :, None] - E[:, None, :]
        F = np.zeros_like(dE)
        mk = df != 0
        F[mk] = df[mk] / dE[mk]
        Mk = U @ (A * np.swapaxes(F, 1, 2)) @ dag(U[i2])
        for p, Gp in enumerate(parts):
            V = 0.5 * (Gp + Gp[i2])
            Z[p, m] = np.einsum("kab,kba->", Mk, V) / Nx
    K = np.fft.fft(Z, axis=1).real / Nx                            # K[x'] = (1/Nx) sum_q Z(q) exp(-i q x')
    return S, K, float(np.abs(E).min())
