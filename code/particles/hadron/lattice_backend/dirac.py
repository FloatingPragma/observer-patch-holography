"""Clover-improved Wilson Dirac operator, CG solver, and point propagators.

Hopping-parameter normalization:

    (D psi)(x) = C(x) psi(x)
                 - kappa sum_mu [ (1 - gamma_mu) U_mu(x) psi(x + mu)
                                 + (1 + gamma_mu) U_mu(x - mu)^dag psi(x - mu) ]

with the clover site matrix

    C(x) = 1 - kappa c_SW sum_{mu<nu} sigma_{mu nu} (x) F_hat_{mu nu}(x),
    F_hat_{mu nu} = (1/(8i)) (Q_{mu nu} - Q_{mu nu}^dag) - trace part,

where Q is the four-leaf clover sum. F_hat is Hermitian traceless, so C is
Hermitian and the operator keeps exact gamma5-hermiticity,
gamma5 D gamma5 = D^dag, which the tests verify numerically. The clover
normalization (unit weight over mu<nu) is a declared convention of this
diagnostic engine.

Fermion boundary conditions: periodic in space, antiperiodic in time,
realized by flipping the sign of the time-direction links on the wrap.
"""

from __future__ import annotations

import numpy as np

from .core import GAMMA, GAMMA5, IDSPIN, NDIM, plaquette_field, shift, sigma_munu

_HOP_PROJ = [(IDSPIN - GAMMA[mu], IDSPIN + GAMMA[mu]) for mu in range(NDIM)]
_SIGMA_MUNU = {(mu, nu): sigma_munu(mu, nu) for mu in range(NDIM) for nu in range(mu + 1, NDIM)}


def antiperiodic_time(u: np.ndarray) -> np.ndarray:
    """Copy of the gauge field with the fermion time boundary sign applied."""
    ubc = u.copy()
    ubc[3, -1] = -ubc[3, -1]
    return ubc


def _apply_spin(mat: np.ndarray, psi: np.ndarray) -> np.ndarray:
    return np.einsum("su,txyzuc->txyzsc", mat, psi)


def _apply_color(link: np.ndarray, psi: np.ndarray) -> np.ndarray:
    return np.einsum("txyzab,txyzsb->txyzsa", link, psi)


def _apply_color_dag(link: np.ndarray, psi: np.ndarray) -> np.ndarray:
    return np.einsum("txyzba,txyzsb->txyzsa", np.conj(link), psi)


def clover_field_strength(u: np.ndarray, mu: int, nu: int) -> np.ndarray:
    """Hermitian traceless clover average F_hat_{mu nu}(x)."""
    p1 = plaquette_field(u, mu, nu)
    # leaf in (nu, -mu): shift the (mu, nu) plaquette taken at x - mu, rotated
    umu, unu = u[mu], u[nu]
    umu_d = np.conj(np.swapaxes(umu, -1, -2))
    unu_d = np.conj(np.swapaxes(unu, -1, -2))
    # P_{nu,-mu}(x) = U_nu(x) U_mu(x+nu-mu)^dag U_nu(x-mu)^dag U_mu(x-mu)
    p2 = unu @ shift(shift(umu_d, nu, +1), mu, -1) @ shift(unu_d, mu, -1) @ shift(umu, mu, -1)
    # P_{-mu,-nu}(x) = U_mu(x-mu)^dag U_nu(x-mu-nu)^dag U_mu(x-mu-nu) U_nu(x-nu)
    p3 = (shift(umu_d, mu, -1) @ shift(shift(unu_d, mu, -1), nu, -1)
          @ shift(shift(umu, mu, -1), nu, -1) @ shift(unu, nu, -1))
    # P_{-nu,mu}(x) = U_nu(x-nu)^dag U_mu(x-nu) U_nu(x+mu-nu) U_mu(x)^dag
    p4 = shift(unu_d, nu, -1) @ shift(umu, nu, -1) @ shift(shift(unu, mu, +1), nu, -1) @ umu_d
    q = p1 + p2 + p3 + p4
    f = (q - np.conj(np.swapaxes(q, -1, -2))) / 8.0j
    tr = np.trace(f, axis1=-2, axis2=-1) / 3.0
    f = f - tr[..., None, None] * np.eye(3)
    return f


def clover_site_matrix(u: np.ndarray, kappa: float, c_sw: float) -> np.ndarray:
    """C(x) as a (site, 4, 3, 4, 3) tensor acting on (spin, color)."""
    shape = u.shape[1:5]
    c = np.zeros((*shape, 4, 3, 4, 3), dtype=complex)
    eye_sc = np.einsum("su,ab->saub", IDSPIN, np.eye(3, dtype=complex))
    c += eye_sc
    if c_sw != 0.0:
        for (mu, nu), sig in _SIGMA_MUNU.items():
            f = clover_field_strength(u, mu, nu)
            c -= kappa * c_sw * np.einsum("su,txyzab->txyzsaub", sig, f)
    return c


def _apply_site_matrix(c: np.ndarray, psi: np.ndarray) -> np.ndarray:
    return np.einsum("txyzsaub,txyzub->txyzsa", c, psi)


class WilsonClover:
    """Wilson-clover operator on a fixed gauge configuration."""

    def __init__(self, u: np.ndarray, kappa: float, c_sw: float = 1.0):
        self.kappa = float(kappa)
        self.c_sw = float(c_sw)
        self.ubc = antiperiodic_time(u)
        self.clover = clover_site_matrix(self.ubc, self.kappa, self.c_sw)

    def apply(self, psi: np.ndarray) -> np.ndarray:
        out = _apply_site_matrix(self.clover, psi)
        for mu in range(NDIM):
            proj_minus, proj_plus = _HOP_PROJ[mu]
            fwd = _apply_color(self.ubc[mu], shift(psi, mu, +1))
            out -= self.kappa * _apply_spin(proj_minus, fwd)
            bwd = _apply_color_dag(shift(self.ubc[mu], mu, -1), shift(psi, mu, -1))
            out -= self.kappa * _apply_spin(proj_plus, bwd)
        return out

    def apply_dag(self, psi: np.ndarray) -> np.ndarray:
        """D^dag psi = gamma5 D gamma5 psi, exact by gamma5-hermiticity."""
        return _apply_spin(GAMMA5, self.apply(_apply_spin(GAMMA5, psi)))

    def apply_normal(self, psi: np.ndarray) -> np.ndarray:
        return self.apply_dag(self.apply(psi))


def cg_normal(
    op: WilsonClover,
    rhs: np.ndarray,
    tol: float = 1e-9,
    max_iter: int = 4000,
) -> tuple[np.ndarray, int, float]:
    """Solve D x = rhs via CG on the normal equations D^dag D x = D^dag rhs.
    Returns (x, iterations, final relative residual of the normal system)."""
    b = op.apply_dag(rhs)
    x = np.zeros_like(rhs)
    r = b.copy()
    p = r.copy()
    rr = float(np.real(np.vdot(r, r)))
    b_norm = float(np.sqrt(np.real(np.vdot(b, b))))
    if b_norm == 0.0:
        return x, 0, 0.0
    it = 0
    for it in range(1, max_iter + 1):
        ap = op.apply_normal(p)
        alpha = rr / float(np.real(np.vdot(p, ap)))
        x += alpha * p
        r -= alpha * ap
        rr_new = float(np.real(np.vdot(r, r)))
        if np.sqrt(rr_new) <= tol * b_norm:
            rr = rr_new
            break
        p = r + (rr_new / rr) * p
        rr = rr_new
    return x, it, float(np.sqrt(rr) / b_norm)


def point_source(shape, spin: int, color: int, t0: int = 0) -> np.ndarray:
    src = np.zeros((*shape, 4, 3), dtype=complex)
    src[t0, 0, 0, 0, spin, color] = 1.0
    return src


def point_propagator(
    op: WilsonClover,
    shape,
    t0: int = 0,
    tol: float = 1e-9,
) -> tuple[np.ndarray, dict]:
    """Full point propagator S[t, x, y, z, s, c, s0, c0] from 12 CG solves."""
    prop = np.zeros((*shape, 4, 3, 4, 3), dtype=complex)
    iters, resids = [], []
    for s0 in range(4):
        for c0 in range(3):
            sol, it, res = cg_normal(op, point_source(shape, s0, c0, t0), tol=tol)
            prop[..., s0, c0] = sol
            iters.append(it)
            resids.append(res)
    info = {"cg_iterations": iters, "cg_relative_residuals": resids}
    return prop, info
