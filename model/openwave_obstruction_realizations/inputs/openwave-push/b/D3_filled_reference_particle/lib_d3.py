"""D3 library: one particle above the free filled reference of the autonomous v5 source, order e^2 (numpy only).

Source: autonomous v5 (44c00889) in its CAR lift (df34a04c), declared infinite-lattice member of D1 (G = n in the
reference, matter at root ports). Reference: the free lower-band Slater filling of K_A(k) = mu beta + sum_a alpha_a
sin k_a with all eight species, tensored with the coexact vacuum. The source owner admits it as a comparison
reference only (6c3dedc4). It is not the interacting ground state. Every function below is about that reference.

Reused, byte-identical copies in consumed/: D1's lib_d1.py with lib21.py and geometry.json (field modes, route
currents, vertex beta_a(k,j) = i sqrt(a/(2 omega)) u^dagger J_a), the row-23 candidate's rc_lib.py (Clifford
matrices, Coulomb kernel v(k) = A [L(k)^-1]_00 through the frozen symbols).

Cell clock: H_cell = h H_phys, mu = h m. Many-fermion generator, expanded in the route amplitude:

    H = sum omega c^dagger c + sum_p psi_p^dagger K_A(p) psi_p                      (free)
      + sum_{p,k,j} c_{kj}^dagger psi_{p-k}^dagger M(p;k,j) psi_p + h.c.            (order e, one quantum)
      - sum_p psi_p^dagger [sum_a (|beta_a|^2/2) alpha_a sin p_a] psi_p             (order e^2, two-quantum part of the hop)
      + C(n) = (a/2) (G-n)^T L^+ (G-n)                                              (order e^2, instantaneous)

    M(p;k,j) = sum_a alpha_a i beta_a(k,j) exp(i k_a/2) cos(p_a - k_a/2)   (D1's one-quantum vertex with z_a = 1)

Order-e^2 self-energy of the state psi^dagger_{K,+,s} |reference> (4 x 4 matrix in spin space, energy E):

    S_C(K)    = (1/2) sum_q v(q) K_A(K-q)/E(K-q)              exchange with the filling plus the particle's own term
    S_sg(K)   = - sum_a (|beta_a|^2/2) alpha_a sin K_a
    S_+(K,E)  = sum_{k,j} M(K)^dagger P_+(K-k) M(K) / (E - E(K-k) - omega)
    S_Z(K,E)  = sum_{k,j} N P_-(K+k) N^dagger / (E + E(K+k) + omega),  N = M(K+k;k,j)   (pair fluctuation blocked by the particle)

Signed carrier without filling (control, D1): S_Z is replaced by S_L = sum M^dagger P_-(K-k) M / (E + E(K-k) - omega)
and S_C by the constant E_self.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'consumed', 'D1_dressed_charge'))
sys.path.insert(0, os.path.join(HERE, 'consumed', 'ROW_23'))
import lib_d1 as d1  # noqa: E402  (D1, byte-identical copy)
import rc_lib as rc  # noqa: E402  (row-23 candidate, byte-identical copy)

A_EL = d1.A_EL
ALPHA = np.array(rc.ALPHA)
BETA = rc.BETA
EYE4 = np.eye(4, dtype=complex)
HS = (1 / 8, 1 / 16, 1 / 32)
_gx, _gw = np.polynomial.legendre.leggauss(80)
_gx, _gw = (_gx + 1) / 2, _gw / 2
# integral of 1/r^2 over the unit cube centred at the origin = 12 int_0^1 int_0^1 du dv / (1 + u^2 + v^2) = 7.674...
# (the first identity run used the constant of the cube [0,1]^3, 2.38, in one informational line; kept as *_first_run_1622)
CUBE = float(12 * np.sum(_gw[:, None] * _gw[None, :] / (1 + _gx[:, None] ** 2 + _gx[None, :] ** 2)))
DIRS = {'x': (1, 0, 0), 'y': (0, 1, 0), 'z': (0, 0, 1), 'xy': (1, 1, 0), 'xyz': (1, 1, 1)}


def unit(name):
    v = np.array(DIRS[name], float)
    return v / np.linalg.norm(v)


def complex_and_paths():
    geom = d1.geometry()
    return d1.lib21.FieldComplex(geom), d1.route_paths(geom), geom


# ---------------------------------------------------------------- quadrature points
def points(L, levels=0, M=11):
    """Momenta and weights: grid L^3 without k = 0; the k = 0 cell is refined `levels` times by an M^3 grid
    (M odd) whose own central cell is left out again. Weights sum to 1 - (L M^levels)^-3."""
    n = np.array(list(np.ndindex(L, L, L)), float)
    k = 2 * np.pi * n / L
    k = (k + np.pi) % (2 * np.pi) - np.pi
    ks = [k[1:]]
    ws = [np.full(len(k) - 1, 1.0 / L ** 3)]
    side = 2 * np.pi / L
    wt = 1.0 / L ** 3
    for _ in range(levels):
        j = np.array(list(np.ndindex(M, M, M)), float) - (M - 1) / 2
        j = j[np.any(j != 0, axis=1)]
        side = side / M
        wt = wt / M ** 3
        ks.append(j * side)
        ws.append(np.full(len(j), wt))
    return np.concatenate(ks), np.concatenate(ws), side


class Field:
    """Per momentum: 24 coexact frequencies W (nk,24), vertex beta (nk,24,3), Coulomb kernel v (nk,), weights."""

    def __init__(self, fc, paths, k, wt, symbols=None):
        self.k = k
        self.wt = wt
        nk = len(k)
        self.W = np.zeros((nk, 24))
        G = np.zeros((nk, 24, 3), complex)
        for i, kk in enumerate(k):
            w, vec = fc.modes(kk)
            assert len(w) == 24, (kk, len(w))
            self.W[i] = w
            for a in range(3):
                G[i, :, a] = vec.conj().T @ fc.current(paths[a], kk)
        self.beta = 1j * np.sqrt(A_EL / (2 * self.W))[:, :, None] * G
        self.v = rc.v_bare(k, symbols) if symbols is not None else np.zeros(nk)
        self.b2 = np.einsum('n,nja->a', wt, abs(self.beta) ** 2)        # |beta_a|^2 of this quadrature

    def off(self):
        """Coupling switched off."""
        f = object.__new__(Field)
        f.k, f.wt, f.W = self.k, self.wt, self.W
        f.beta = np.zeros_like(self.beta)
        f.v = np.zeros_like(self.v)
        f.b2 = np.zeros(3)
        return f

    def subset(self, n):
        """First n momenta (the levels are stored in order: a coarser nesting is a prefix)."""
        f = object.__new__(Field)
        f.k, f.wt, f.W, f.beta, f.v = self.k[:n], self.wt[:n], self.W[:n], self.beta[:n], self.v[:n]
        f.b2 = np.einsum('n,nja->a', f.wt, abs(f.beta) ** 2)
        return f


# ---------------------------------------------------------------- free matter
def KA(q, mu, z=(1.0, 1.0, 1.0)):
    q = np.asarray(q, float)
    out = np.broadcast_to(mu * BETA, q.shape[:-1] + (4, 4)).copy()
    for a in range(3):
        out = out + (z[a] * np.sin(q[..., a]))[..., None, None] * ALPHA[a]
    return out


def EQ(q, mu, z=(1.0, 1.0, 1.0)):
    q = np.asarray(q, float)
    return np.sqrt(mu * mu + np.sum((np.asarray(z) * np.sin(q)) ** 2, axis=-1))


def upper_basis(K, mu, z=(1.0, 1.0, 1.0)):
    """Orthonormal basis (4 x 2) of the upper band at K, from the two beta = +1 columns."""
    P = (EYE4 + KA(K, mu, z) / EQ(K, mu, z)) / 2
    q, _ = np.linalg.qr(P[:, [0, 2]])
    return q


# ---------------------------------------------------------------- self-energy
def sandwich(Cs, P):
    """For each C (nk,3,3) in Cs: sum_n sum_{a,c} C[n,a,c] alpha_a P[n] alpha_c (4 x 4)."""
    AP = np.einsum('axy,nyz->naxz', ALPHA, P)
    return [np.einsum('ncxz,czw->xw', np.einsum('nac,naxz->ncxz', C, AP), ALPHA) for C in Cs]


def sigma(fd, K, mu, E=None, filled=True, zk=(1.0, 1.0, 1.0), zv=(1.0, 1.0, 1.0), eta=0.0, sign_exchange=1, sign_pair=1):
    """Order-e^2 self-energy matrices at momentum K and energy E (default: the free upper-branch energy).
    filled=False: signed carrier without filling (lower-branch term with broadening eta, constant E_self).
    zk, zv: hop factors in the kinematics and in the vertex (1 at strict order e^2; D1's convention uses z_a)."""
    K = np.asarray(K, float)
    zk = np.asarray(zk, float)
    zv = np.asarray(zv, float)
    E0 = float(EQ(K, mu, zk))
    if E is None:
        E = E0
    k, wt, W, beta = fd.k, fd.wt, fd.W, fd.beta
    ph = np.exp(0.5j * k)
    out = {'E0': E0}
    # particle emits one quantum and stays in the upper branch
    m = 1j * beta * (ph * np.cos(K[None, :] - k / 2) * zv[None, :])[:, None, :]
    Ef = EQ(K - k, mu, zk)
    HKf = KA(K - k, mu, zk) / Ef[:, None, None]
    Pp = (EYE4 + HKf) / 2
    Pm = (EYE4 - HKf) / 2
    den = E - Ef[:, None] - W
    mh = np.swapaxes(np.conj(m), 1, 2)                                   # (nk,3,24)

    def cmat(f):                                                         # sum_j wt f conj(m_a) m_c -> (nk,3,3)
        return np.matmul(mh * (wt[:, None] * f)[:, None, :], m)

    C = cmat(1 / den)
    dC = cmat(-1 / den ** 2)
    out['S_plus'], out['dS_plus'] = sandwich([C, dC], Pp)
    out['plus_mismatch_max'] = float(den.max())
    out['plus_open_modes'] = int(np.sum(den >= 0))
    if filled:
        nv = 1j * beta * (ph * np.cos(K[None, :] + k / 2) * zv[None, :])[:, None, :]
        Ez = EQ(K + k, mu, zk)
        Pmz = (EYE4 - KA(K + k, mu, zk) / Ez[:, None, None]) / 2
        denz = E + Ez[:, None] + W
        nt = np.swapaxes(nv, 1, 2)                                       # (nk,3,24): n_a conj(n_c)
        D = np.matmul(nt * (wt[:, None] / denz)[:, None, :], np.conj(nv))
        dD = np.matmul(nt * (-wt[:, None] / denz ** 2)[:, None, :], np.conj(nv))
        sp, dsp = sandwich([D, dD], Pmz)
        out['S_pair'], out['dS_pair'] = sign_pair * sp, sign_pair * dsp
        sc = 0.5 * np.einsum('n,nxy->xy', wt * fd.v, HKf)
        if sign_exchange == 1:
            out['S_coulomb'] = sc
        else:                                   # wrong sign of the exchange term: E_self + sum v P_-
            out['S_coulomb'] = float(np.sum(wt * fd.v)) * EYE4 - sc
    else:
        denl = E + Ef[:, None] - W
        pv = denl / (denl ** 2 + eta ** 2) if eta > 0 else 1 / denl
        lor = (eta / np.pi) / (denl ** 2 + eta ** 2) if eta > 0 else np.zeros_like(denl)
        Cl = cmat(pv)
        Gl = cmat(2 * np.pi * lor)
        Wl = cmat(np.ones_like(denl))
        out['S_lower'], out['Gamma_lower'], out['W_lower'] = sandwich([Cl, Gl, Wl], Pm)
        out['W_upper'] = sandwich([Wl], Pp)[0]
        out['lower_mismatch_min'] = float(denl.min())
        out['lower_mismatch_max'] = float(denl.max())
        out['S_coulomb'] = 0.5 * float(np.sum(wt * fd.v)) * EYE4
    out['S_seagull'] = -np.einsum('a,axy->xy', 0.5 * fd.b2 * np.sin(K) * zk, ALPHA)
    return out


def block(S, U):
    return U.conj().T @ S @ U


def shift(fd, K, mu, **kw):
    """Strict order e^2: eigenvalues of the 2 x 2 upper-band block of the on-shell self-energy, parts, and Z."""
    s = sigma(fd, K, mu, **kw)
    U = upper_basis(K, mu, kw.get('zk', (1.0, 1.0, 1.0)))
    names = [n for n in ('S_coulomb', 'S_seagull', 'S_plus', 'S_pair', 'S_lower') if n in s]
    tot = sum(s[n] for n in names)
    B = block(tot, U)
    ev = np.linalg.eigvalsh((B + B.conj().T) / 2)
    parts = {n: float(np.trace(block(s[n], U)).real / 2) for n in names}
    dn = [n for n in ('dS_plus', 'dS_pair') if n in s]
    dE = float(sum(np.trace(block(s[n], U)).real for n in dn) / 2)
    return {'E0': s['E0'], 'mean': float(ev.mean()), 'split': float(ev[1] - ev[0]), 'parts': parts,
            'dSigma_dE': dE, 'Z_strict': 1 + dE, 'Z_resummed': 1 / (1 - dE), 'anti_hermitian': float(abs(B - B.conj().T).max()),
            'plus_mismatch_max': s['plus_mismatch_max'], 'plus_open_modes': s['plus_open_modes'], 'raw': s}


def curvature(fd, direction, mu, delta=None, **kw):
    """Second derivative along a unit direction at K = 0 of the mean order-e^2 shift (and of each part), by the
    fourth-order even stencil at 0, delta, 2 delta; also the spin splitting at delta."""
    n = np.asarray(direction, float)
    if delta is None:
        delta = 0.08 * mu
    s0 = shift(fd, np.zeros(3), mu, **kw)
    s1 = shift(fd, delta * n, mu, **kw)
    s2 = shift(fd, 2 * delta * n, mu, **kw)

    def d2(f0, f1, f2):
        return 2 * (16 * (f1 - f0) - (f2 - f0)) / (12 * delta ** 2)

    parts = {p: d2(s0['parts'][p], s1['parts'][p], s2['parts'][p]) for p in s0['parts']}
    return {'delta': delta, 'rest_shift': s0['mean'], 'rest_parts': s0['parts'], 'Z_strict_rest': s0['Z_strict'],
            'Z_resummed_rest': s0['Z_resummed'], 'curv_mean': d2(s0['mean'], s1['mean'], s2['mean']), 'curv_parts': parts,
            'split_at_delta': s1['split'], 'split_at_2delta': s2['split'], 'rest_split': s0['split']}


def speeds(mu, rest_shift, curv):
    """c^2 = (rest energy) x (curvature of E): strict order e^2 and the unexpanded product."""
    lin = 1 + rest_shift / mu + mu * curv
    full = (mu + rest_shift) * (1 / mu + curv)
    return {'c2_strict': lin, 'c_strict': 1 + 0.5 * (rest_shift / mu + mu * curv), 'c2_product': full,
            'c_product': math.sqrt(full) if full > 0 else float('nan'),
            'kinetic_mass_over_m_strict': 1 - mu * curv, 'kinetic_mass_over_m_resummed': 1 / (1 + mu * curv) if 1 + mu * curv > 0 else float('nan')}


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items() if k != 'raw'}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, complex):
        return [o.real, o.imag]
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


def dump(obj, name):
    with open(os.path.join(HERE, name), 'w') as fh:
        json.dump(jsonable(obj), fh, indent=1)
        fh.write('\n')
