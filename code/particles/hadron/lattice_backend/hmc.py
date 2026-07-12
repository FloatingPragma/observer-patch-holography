"""Two-flavor pseudofermion HMC for the dynamical branch of the engine.

The demonstration Hamiltonian is

    H = sum_links Tr[pi^2] + S_g(U; beta) + phi^dag (D^dag D)^-1 phi

with pi Hermitian traceless momenta (pi = sum_a pi^a T^a, pi^a ~ N(0,1),
Tr[T^a T^b] = delta^{ab}/2), the Wilson plaquette gauge action, and the
unimproved Wilson operator (c_SW = 0) for the degenerate light doublet.
The clover term is used in valence measurements only; carrying its force in
the toy dynamical branch is a declared simplification of this diagnostic
engine, recorded in the export manifest.

Flow equations (verified by the finite-difference force test):

    dU/dt  = i pi U
    dpi/dt = -(F_gauge + F_fermion) / 2

where F is the Hermitian traceless matrix with dS/dtau = Tr[Q F] for the
link rotation U -> exp(i tau Q) U.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import NDIM, staple_sum, wilson_gauge_action
from .dirac import WilsonClover, cg_normal, _apply_spin, _apply_color
from .core import GAMMA, IDSPIN, shift


def _herm_traceless(m: np.ndarray) -> np.ndarray:
    h = 0.5 * (m + np.conj(np.swapaxes(m, -1, -2)))
    tr = np.trace(h, axis1=-2, axis2=-1) / 3.0
    return h - tr[..., None, None] * np.eye(3)


def _im_herm(m: np.ndarray) -> np.ndarray:
    """(M - M^dag)/(2i), traceless part."""
    a = (m - np.conj(np.swapaxes(m, -1, -2))) / 2.0j
    tr = np.trace(a, axis1=-2, axis2=-1) / 3.0
    return a - tr[..., None, None] * np.eye(3)


def sample_momenta(rng: np.random.Generator, shape) -> np.ndarray:
    """Gaussian su(3) momenta as Hermitian traceless matrices."""
    a = rng.normal(size=(NDIM, *shape, 3, 3)) + 1j * rng.normal(size=(NDIM, *shape, 3, 3))
    return _herm_traceless(a / np.sqrt(2.0))


def kinetic_energy(pi: np.ndarray) -> float:
    return float(np.real(np.einsum("...ij,...ji->...", pi, pi).sum()))


def exp_i_herm(q: np.ndarray) -> np.ndarray:
    """exp(i Q) for Hermitian Q via eigendecomposition."""
    w, v = np.linalg.eigh(q)
    return np.einsum("...ij,...j,...kj->...ik", v, np.exp(1j * w), np.conj(v))


def gauge_force(u: np.ndarray, beta: float) -> np.ndarray:
    """F with dS_g/dtau = Tr[Q F] under U_mu(x) -> exp(i tau Q) U_mu(x)."""
    f = np.empty_like(u)
    for mu in range(NDIM):
        f[mu] = (beta / 3.0) * _im_herm(u[mu] @ staple_sum(u, mu))
    return f


def pseudofermion_force(
    op: WilsonClover,
    u_bc: np.ndarray,
    x_sol: np.ndarray,
) -> np.ndarray:
    """Fermion force for S_pf = phi^dag (D^dag D)^-1 phi at X = (D^dag D)^-1 phi.

    dS = -2 Re[Y^dag dD X] with Y = D X. For each link, the forward hop gives
    -i kappa Tr[Q B1] and the backward hop +i kappa Tr[Q C2]; the force is
    F = -2 kappa Im_H(B1) + 2 kappa Im_H(C2).
    """
    kappa = op.kappa
    y_sol = op.apply(x_sol)
    f = np.zeros((NDIM, *x_sol.shape[:4], 3, 3), dtype=complex)
    for mu in range(NDIM):
        proj_minus = IDSPIN - GAMMA[mu]
        proj_plus = IDSPIN + GAMMA[mu]
        link = u_bc[mu]
        x_fwd = shift(x_sol, mu, +1)
        y_fwd = shift(y_sol, mu, +1)
        # B1_{ba} = sum_s [P- U X(x+mu)]_{s b} conj(Y(x))_{s a}
        b1 = np.einsum(
            "txyzsb,txyzsa->txyzba",
            _apply_spin(proj_minus, _apply_color(link, x_fwd)),
            np.conj(y_sol),
        )
        # C2_{cb} = sum_s X(x)_{s c} [conj(P+ Y(x+mu)) U^dag]_{s b}
        z2 = _apply_spin(proj_plus, y_fwd)
        w2 = np.einsum("txyzsa,txyzba->txyzsb", np.conj(z2), np.conj(link))
        c2 = np.einsum("txyzsc,txyzsb->txyzcb", x_sol, w2)
        f[mu] = -2.0 * kappa * _im_herm(b1) + 2.0 * kappa * _im_herm(c2)
    return f


@dataclass
class HMCResult:
    accepted: bool
    delta_h: float
    action_gauge: float
    action_pf: float
    cg_iterations: int


class TwoFlavorHMC:
    """Leapfrog HMC for the degenerate light doublet on the Wilson action."""

    def __init__(
        self,
        beta: float,
        kappa: float,
        n_steps: int = 20,
        trajectory_length: float = 1.0,
        cg_tol: float = 1e-10,
    ):
        self.beta = beta
        self.kappa = kappa
        self.n_steps = n_steps
        self.dt = trajectory_length / n_steps
        self.cg_tol = cg_tol

    def _pf_solution(self, u: np.ndarray, phi: np.ndarray) -> tuple[np.ndarray, WilsonClover, int]:
        """X = (D^dag D)^-1 phi by CG on the normal operator."""
        op = WilsonClover(u, self.kappa, c_sw=0.0)
        b = phi
        x = np.zeros_like(phi)
        r = b.copy()
        p = r.copy()
        rr = float(np.real(np.vdot(r, r)))
        b_norm = float(np.sqrt(rr))
        it = 0
        if b_norm > 0.0:
            for it in range(1, 20000):
                ap = op.apply_normal(p)
                alpha = rr / float(np.real(np.vdot(p, ap)))
                x += alpha * p
                r -= alpha * ap
                rr_new = float(np.real(np.vdot(r, r)))
                if np.sqrt(rr_new) <= self.cg_tol * b_norm:
                    break
                p = r + (rr_new / rr) * p
                rr = rr_new
        return x, op, it

    def _force(self, u: np.ndarray, phi: np.ndarray) -> tuple[np.ndarray, int]:
        x_sol, op, iters = self._pf_solution(u, phi)
        f = gauge_force(u, self.beta) + pseudofermion_force(op, op.ubc, x_sol)
        return f, iters

    def action_pf(self, u: np.ndarray, phi: np.ndarray) -> float:
        x_sol, _, _ = self._pf_solution(u, phi)
        return float(np.real(np.vdot(phi, x_sol)))

    def integrate(
        self,
        u: np.ndarray,
        pi: np.ndarray,
        phi: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, int]:
        """Leapfrog integration of one trajectory; returns (U', pi', cg_count).
        Time-reversible: integrating (U', -pi') returns to (U, -pi)."""
        u_new = u.copy()
        total_cg = 0
        f, it = self._force(u_new, phi)
        total_cg += it
        pi = pi - 0.5 * self.dt * (f / 2.0)
        for step in range(self.n_steps):
            u_new = exp_i_herm(self.dt * pi) @ u_new
            f, it = self._force(u_new, phi)
            total_cg += it
            weight = 0.5 if step == self.n_steps - 1 else 1.0
            pi = pi - weight * self.dt * (f / 2.0)
        return u_new, pi, total_cg

    def trajectory(
        self,
        rng: np.random.Generator,
        u: np.ndarray,
    ) -> tuple[np.ndarray, HMCResult]:
        shape = u.shape[1:5]
        pi = sample_momenta(rng, shape)
        # heatbath pseudofermion: phi = D^dag eta with Gaussian eta
        eta = (rng.normal(size=(*shape, 4, 3)) + 1j * rng.normal(size=(*shape, 4, 3))) / np.sqrt(2.0)
        op0 = WilsonClover(u, self.kappa, c_sw=0.0)
        phi = op0.apply_dag(eta)

        s_g0 = wilson_gauge_action(u, self.beta)
        s_pf0 = float(np.real(np.vdot(eta, eta)))
        h0 = kinetic_energy(pi) + s_g0 + s_pf0

        u_new, pi, total_cg = self.integrate(u, pi, phi)

        s_g1 = wilson_gauge_action(u_new, self.beta)
        s_pf1 = self.action_pf(u_new, phi)
        h1 = kinetic_energy(pi) + s_g1 + s_pf1
        delta_h = h1 - h0
        accept = (delta_h <= 0.0) or (rng.uniform() < np.exp(-delta_h))
        result = HMCResult(
            accepted=bool(accept),
            delta_h=float(delta_h),
            action_gauge=float(s_g1 if accept else s_g0),
            action_pf=float(s_pf1 if accept else s_pf0),
            cg_iterations=total_cg,
        )
        return (u_new if accept else u), result
