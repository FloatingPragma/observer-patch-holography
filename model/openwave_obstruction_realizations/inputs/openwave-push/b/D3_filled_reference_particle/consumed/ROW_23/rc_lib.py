"""Row 23 library (producer P23, second Fable session).

Matter sector of source `common_v4`: the routed Dirac macro step
    W = S_z S_y S_x exp(-i mu beta),   S_a psi(x) = P_a psi(x-e_a) + Q_a psi(x+e_a),
with P_a = (1+alpha_a)/2, Q_a = 1-P_a, in the units one macro = one time
unit, one cell = one length unit (mu = h*m is the mass-coin angle per macro).
The static held charge acts through the frozen longitudinal phase
exp(-i A q^T L^+ q / 2), A = 1/2, of the LOCAL12 graph: on a matter unit
(q_charge = -1) at the root port of cell x it is the diagonal kick
exp(-i V(x)), V = -A sum_y G_00(x-y) Q(y).

Own implementation; the frozen producers are imported read-only in
`frozen_modules` for the identity gates only.
"""
from pathlib import Path
import hashlib
import importlib.util
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
META = Path('/Users/muellerberndt/Projects/oph-meta')
TUNING = META/'openwave-tuning'
A_ELECTRIC = 0.5
E2_BARE = 1.0/3.0

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1., -1.]).astype(complex)
I2 = np.eye(2, dtype=complex)
# frozen little-endian labels: 'XX','XY','XZ' -> kron(sigma_a, X); 'ZI' -> kron(I, Z)
ALPHA = (np.kron(SX, SX), np.kron(SY, SX), np.kron(SZ, SX))
BETA = np.kron(I2, SZ)
EYE4 = np.eye(4, dtype=complex)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- Bloch symbol
def wilson_weight(k):
    return np.sum(1.0-np.cos(np.asarray(k, float)), axis=-1)


def bloch(k, mu, r=0.0):
    """U(k) = e^{-i k_z a_z} e^{-i k_y a_y} e^{-i k_x a_x} e^{-i M beta}.

    r = 0 is the source.  r != 0 is the mutated law of control K2
    (mass-coin angle M = mu + r sum_a (1 - cos k_a)); it is not the source.
    """
    k = np.asarray(k, float)
    M = mu + r*wilson_weight(k)
    U = np.cos(M)[..., None, None]*EYE4 - 1j*np.sin(M)[..., None, None]*BETA
    for a in range(3):
        F = np.cos(k[..., a])[..., None, None]*EYE4 - 1j*np.sin(k[..., a])[..., None, None]*ALPHA[a]
        U = F @ U
    return U


def grid(L):
    k = 2*np.pi*np.arange(L)/L
    return np.stack(np.meshgrid(k, k, k, indexing='ij'), axis=-1)


# ---------------------------------------------------------------- real space
def step(state, mu):
    """One free macro on an array of shape (L,L,L,4,...)."""
    coin = np.cos(mu)*EYE4 - 1j*np.sin(mu)*BETA
    out = np.einsum('ab,xyzb...->xyza...', coin, state)
    for a in range(3):
        P = (EYE4+ALPHA[a])/2
        Q = EYE4-P
        out = (np.einsum('ab,xyzb...->xyza...', P, np.roll(out, 1, axis=a))
               + np.einsum('ab,xyzb...->xyza...', Q, np.roll(out, -1, axis=a)))
    return out


def step_matrix(L, mu, r=0.0):
    """Dense one-body macro step, index ((x*L+y)*L+z)*4+s."""
    N = 4*L**3
    ident = np.eye(N, dtype=complex).reshape(L, L, L, 4, N)
    if r != 0.0:
        w = wilson_weight(grid(L))
        f = np.fft.fftn(ident, axes=(0, 1, 2))
        sign = np.real(np.diag(BETA))
        phase = np.exp(-1j*r*w[..., None]*sign[None, None, None, :])
        ident = np.fft.ifftn(f*phase[..., None], axes=(0, 1, 2))
    return step(ident, mu).reshape(N, N)


def principal_negative_projector(U):
    """Projector on the eigenmodes with quasi-energy in (-pi,0) (U = e^{-i eps}).

    eps in (-pi,0)  <=>  -sin(eps) > 0  <=>  positive part of (U-U^dagger)/(2i).
    """
    Hs = (U-U.conj().T)/(2j)
    w, V = np.linalg.eigh(Hs)
    pos = w > 0
    Vp = V[:, pos]
    return Vp @ Vp.conj().T, int(pos.sum()), float(np.min(np.abs(w)))


def cell_density(P, L):
    return np.real(np.diag(P)).reshape(L, L, L, 4).sum(axis=-1)


def kick_matrix_diag(V):
    """Diagonal of exp(-i V(x)) on (cell, spinor)."""
    return np.repeat(np.exp(-1j*np.asarray(V, float).reshape(-1)), 4)


# ---------------------------------------------------------------- frozen static generator
def codex_symbols():
    path = TUNING/'next/CHAIN3/coulomb_response_connection/run.py'
    return load_module('codex_coulomb_response_run', path), path


def green_root(L, symbols):
    """Root-port Green function of the LOCAL12 graph Laplacian on the L^3 torus.

    Returns (g_tilde(k) real array, G_00(x) real array); zero mode removed
    (pseudo-inverse), as in the frozen static generator.
    """
    H = symbols(grid(L)).copy()
    H[0, 0, 0] += np.ones((12, 12))/12
    inv = np.linalg.inv(H)
    inv[0, 0, 0] -= np.ones((12, 12))/12
    gt = inv[..., 0, 0]
    assert np.max(abs(gt.imag)) < 1e-12
    gt = gt.real.copy()
    G = np.fft.ifftn(gt)
    assert np.max(abs(G.imag)) < 1e-12
    return gt, G.real.copy()


def v_bare(k, symbols):
    """v(k) = A [L(k)^-1]_00 for nonzero k (array of shape (...,3))."""
    H = symbols(np.asarray(k, float))
    inv = np.linalg.inv(H)
    return A_ELECTRIC*inv[..., 0, 0].real


# ---------------------------------------------------------------- one-loop function
_X, _W = np.polynomial.legendre.leggauss(400)
_X = (_X+1)/2
_W = _W/2
_XX = _X*(1-_X)


def uehling_I(t):
    """I(t) = int_0^1 x(1-x) ln(1 + x(1-x) t) dx."""
    t = np.asarray(t, float)
    return np.sum(_W*_XX*np.log1p(_XX*t[..., None]), axis=-1)


def uehling_slope(t):
    """2 t I'(t): logarithmic slope dI/dln q; tends to 1/3."""
    t = np.asarray(t, float)
    return 2*np.sum(_W*_XX*_XX*t[..., None]/(1+_XX*t[..., None]), axis=-1)


def half_slope_t():
    lo, hi = 1e-3, 1e4
    for _ in range(200):
        mid = np.sqrt(lo*hi)
        if uehling_slope(mid) < 1/6:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo*hi))
