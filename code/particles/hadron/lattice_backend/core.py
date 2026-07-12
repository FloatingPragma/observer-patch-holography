"""SU(3) lattice gauge core for the real diagnostic hadron engine.

Conventions follow OPH_HADRON_BACKEND_REFERENCE_SPEC.md where it fixes them:
isotropic Wilson plaquette action, SU(3) fundamental representation, periodic
gauge boundary conditions in all four directions, Euclidean Dirac matrices
with gamma_4 Hermitian and gamma_5 = gamma_1 gamma_2 gamma_3 gamma_4.

Field layout:
    U[mu, t, x, y, z]      -> complex128 (3, 3) link matrices
    psi[t, x, y, z, s, c]  -> spin s in 0..3, color c in 0..2

Site shifts use numpy.roll on the axis of direction mu (axis mu + 1 for gauge
fields, axis mu for spinors), which realizes the periodic lattice.
"""

from __future__ import annotations

import numpy as np

NDIM = 4
NCOL = 3
NSPIN = 4

# Euclidean Dirac-Pauli basis: gamma_k = [[0, -i sigma_k], [i sigma_k, 0]],
# gamma_4 = diag(1, 1, -1, -1). All four are Hermitian.
_SIGMA = [
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
]


def _gamma_spatial(k: int) -> np.ndarray:
    g = np.zeros((4, 4), dtype=complex)
    g[:2, 2:] = -1j * _SIGMA[k]
    g[2:, :2] = 1j * _SIGMA[k]
    return g


GAMMA = np.array([
    _gamma_spatial(0),
    _gamma_spatial(1),
    _gamma_spatial(2),
    np.diag([1.0, 1.0, -1.0, -1.0]).astype(complex),
])
GAMMA5 = GAMMA[0] @ GAMMA[1] @ GAMMA[2] @ GAMMA[3]
IDSPIN = np.eye(NSPIN, dtype=complex)
# charge conjugation C = gamma_4 gamma_2: C gamma_mu C^-1 = -gamma_mu^T
CCONJ = GAMMA[3] @ GAMMA[1]


def sigma_munu(mu: int, nu: int) -> np.ndarray:
    """sigma_{mu nu} = (i/2) [gamma_mu, gamma_nu]."""
    return 0.5j * (GAMMA[mu] @ GAMMA[nu] - GAMMA[nu] @ GAMMA[mu])


def cold_start(shape: tuple[int, int, int, int]) -> np.ndarray:
    unit = np.zeros((NDIM, *shape, NCOL, NCOL), dtype=complex)
    unit[..., np.arange(NCOL), np.arange(NCOL)] = 1.0
    return unit


def project_su3(m: np.ndarray) -> np.ndarray:
    """Project matrices onto SU(3) by polar decomposition and phase removal."""
    u, _, vh = np.linalg.svd(m)
    w = u @ vh
    det = np.linalg.det(w)
    return w * (det ** (-1.0 / 3.0))[..., None, None]


def random_su3_near_identity(rng: np.random.Generator, shape, eps: float) -> np.ndarray:
    """Random SU(3) elements exp(i eps H) with H Hermitian traceless."""
    a = rng.normal(size=(*shape, NCOL, NCOL)) + 1j * rng.normal(size=(*shape, NCOL, NCOL))
    h = 0.5 * (a + np.conj(np.swapaxes(a, -1, -2)))
    tr = np.trace(h, axis1=-2, axis2=-1) / NCOL
    h = h - tr[..., None, None] * np.eye(NCOL)
    w, v = np.linalg.eigh(h)
    phase = np.exp(1j * eps * w)
    return np.einsum("...ij,...j,...kj->...ik", v, phase, np.conj(v))


def _mdag(m: np.ndarray) -> np.ndarray:
    return np.conj(np.swapaxes(m, -1, -2))


def shift(field: np.ndarray, mu: int, sign: int, site_axis0: int = 0) -> np.ndarray:
    """Field at x + sign*mu_hat. site_axis0 is the axis index of the t axis."""
    return np.roll(field, -sign, axis=site_axis0 + mu)


def plaquette_field(u: np.ndarray, mu: int, nu: int) -> np.ndarray:
    """U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag at every site."""
    umu = u[mu]
    unu = u[nu]
    unu_xpmu = shift(unu, mu, +1)
    umu_xpnu = shift(umu, nu, +1)
    return umu @ unu_xpmu @ _mdag(umu_xpnu) @ _mdag(unu)


def average_plaquette(u: np.ndarray) -> float:
    """<(1/3) Re Tr U_p> over all sites and the 6 plaquette orientations."""
    total = 0.0
    count = 0
    for mu in range(NDIM):
        for nu in range(mu + 1, NDIM):
            p = plaquette_field(u, mu, nu)
            total += float(np.real(np.trace(p, axis1=-2, axis2=-1)).sum())
            count += p[..., 0, 0].size
    return total / (NCOL * count)


def wilson_gauge_action(u: np.ndarray, beta: float) -> float:
    """S_g = beta sum_p (1 - (1/3) Re Tr U_p), per the reference spec."""
    volume = u[0][..., 0, 0].size
    n_plaq = 6 * volume
    return beta * n_plaq * (1.0 - average_plaquette(u))


def staple_sum(u: np.ndarray, mu: int) -> np.ndarray:
    """Sum of the six staples around U_mu(x); S_g contains
    -(beta/3) Re Tr [U_mu(x) staple^dag] up to the constant term."""
    total = np.zeros_like(u[mu])
    for nu in range(NDIM):
        if nu == mu:
            continue
        unu = u[nu]
        umu = u[mu]
        # forward staple: U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag
        total += shift(unu, mu, +1) @ _mdag(shift(umu, nu, +1)) @ _mdag(unu)
        # backward staple: U_nu(x+mu-nu)^dag U_mu(x-nu)^dag U_nu(x-nu)
        unu_back = shift(unu, nu, -1)
        umu_back = shift(umu, nu, -1)
        total += _mdag(shift(unu_back, mu, +1)) @ _mdag(umu_back) @ unu_back
    return total


def gauge_transform(u: np.ndarray, g: np.ndarray) -> np.ndarray:
    """U_mu(x) -> g(x) U_mu(x) g(x+mu)^dag."""
    out = np.empty_like(u)
    for mu in range(NDIM):
        out[mu] = g @ u[mu] @ _mdag(shift(g, mu, +1, site_axis0=0))
    return out


def gauge_transform_spinor(psi: np.ndarray, g: np.ndarray) -> np.ndarray:
    """psi(x) -> g(x) psi(x)."""
    return np.einsum("txyzab,txyzsb->txyzsa", g, psi)


# ---------------------------------------------------------------------------
# Quenched updates: Cabibbo-Marinari heatbath + overrelaxation on the three
# SU(2) subgroups, applied checkerboard-wise per direction so every link in a
# sub-step sees a fixed staple.
# ---------------------------------------------------------------------------

_SU2_SUBGROUPS = ((0, 1), (0, 2), (1, 2))


def _parity_mask(shape: tuple[int, int, int, int], parity: int) -> np.ndarray:
    t, x, y, z = np.indices(shape)
    return ((t + x + y + z) % 2) == parity


def _su2_extract(w: np.ndarray, i: int, j: int) -> np.ndarray:
    """Least-squares SU(2)-proportional part of the (i, j) 2x2 submatrix,
    returned as coefficients (a0, a1, a2, a3) with w ~ a0 + i a_k sigma_k."""
    a0 = 0.5 * np.real(w[..., i, i] + w[..., j, j])
    a1 = 0.5 * np.imag(w[..., i, j] + w[..., j, i])
    a2 = 0.5 * np.real(w[..., i, j] - w[..., j, i])
    a3 = 0.5 * np.imag(w[..., i, i] - w[..., j, j])
    return np.stack([a0, a1, a2, a3], axis=-1)


def _su2_to_matrix(a: np.ndarray, i: int, j: int, ncol: int = NCOL) -> np.ndarray:
    """Embed SU(2) coefficients into an SU(3) matrix on rows/cols (i, j)."""
    out = np.zeros((*a.shape[:-1], ncol, ncol), dtype=complex)
    for d in range(ncol):
        out[..., d, d] = 1.0
    out[..., i, i] = a[..., 0] + 1j * a[..., 3]
    out[..., i, j] = a[..., 2] + 1j * a[..., 1]
    out[..., j, i] = -a[..., 2] + 1j * a[..., 1]
    out[..., j, j] = a[..., 0] - 1j * a[..., 3]
    return out


def _kennedy_pendleton(rng: np.random.Generator, k: np.ndarray) -> np.ndarray:
    """Sample a0 with density ~ sqrt(1-a0^2) exp(k a0) for the SU(2) heatbath."""
    a0 = np.empty_like(k)
    todo = np.ones(k.shape, dtype=bool)
    # guard against tiny k: fall back to uniform a0 in [-1, 1]
    small = k < 1e-8
    a0[small] = rng.uniform(-1.0, 1.0, size=int(small.sum()))
    todo &= ~small
    while todo.any():
        n = int(todo.sum())
        r1 = rng.uniform(size=n)
        r2 = rng.uniform(size=n)
        r3 = rng.uniform(size=n)
        kk = k[todo]
        x = 1.0 + (np.log(r1) + (np.cos(2.0 * np.pi * r2) ** 2) * np.log(r3)) / kk
        accept = (rng.uniform(size=n) ** 2) <= 1.0 - 0.5 * (1.0 - x)
        # numerical guard: reject x outside [-1, 1]
        accept &= x >= -1.0
        idx = np.flatnonzero(todo)
        good = idx[accept]
        vals = x[accept]
        a0.flat[good] = vals
        todo.flat[good] = False
    return a0


def _su2_heatbath(rng: np.random.Generator, a: np.ndarray, beta_eff: float) -> np.ndarray:
    """Heatbath for SU(2) with distribution ~ exp(beta_eff/2 Tr(g v^dag)),
    where a are the coefficients of v = a0 + i a_k sigma_k (not normalized)."""
    norm = np.sqrt(np.maximum(np.sum(a * a, axis=-1), 1e-300))
    k = beta_eff * norm
    a0 = _kennedy_pendleton(rng, k)
    r = np.sqrt(np.maximum(1.0 - a0 * a0, 0.0))
    ct = rng.uniform(-1.0, 1.0, size=a0.shape)
    st = np.sqrt(np.maximum(1.0 - ct * ct, 0.0))
    ph = rng.uniform(0.0, 2.0 * np.pi, size=a0.shape)
    g = np.stack([a0, r * st * np.cos(ph), r * st * np.sin(ph), r * ct], axis=-1)
    # new link submatrix g v^-1_normalized^dag: compose g with v^dag / |v|
    vd = np.stack([a[..., 0], -a[..., 1], -a[..., 2], -a[..., 3]], axis=-1) / norm[..., None]
    return _su2_compose(g, vd)


def _su2_compose(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Quaternion product of SU(2) coefficient vectors."""
    p0, p1, p2, p3 = (p[..., d] for d in range(4))
    q0, q1, q2, q3 = (q[..., d] for d in range(4))
    return np.stack([
        p0 * q0 - p1 * q1 - p2 * q2 - p3 * q3,
        p0 * q1 + p1 * q0 - p2 * q3 + p3 * q2,
        p0 * q2 + p1 * q3 + p2 * q0 - p3 * q1,
        p0 * q3 - p1 * q2 + p2 * q1 + p3 * q0,
    ], axis=-1)


def _su2_or_reflect(a: np.ndarray) -> np.ndarray:
    """Overrelaxation reflection g -> v_norm^dag g^dag v_norm^dag maximizing
    the action-preserving move; here reduced to g_new = vd * g^dag * vd with
    vd the normalized conjugate of the staple projection."""
    norm = np.sqrt(np.maximum(np.sum(a * a, axis=-1), 1e-300))
    vd = np.stack([a[..., 0], -a[..., 1], -a[..., 2], -a[..., 3]], axis=-1) / norm[..., None]
    return _su2_compose(vd, vd)


def _apply_su2_update(
    rng: np.random.Generator,
    u_mu: np.ndarray,
    staple: np.ndarray,
    mask: np.ndarray,
    beta: float,
    overrelax: bool,
) -> None:
    """One Cabibbo-Marinari pass over the three SU(2) subgroups on masked sites."""
    for (i, j) in _SU2_SUBGROUPS:
        w = u_mu @ staple
        a = _su2_extract(w[mask], i, j)
        if overrelax:
            g = _su2_or_reflect(a)
        else:
            g = _su2_heatbath(rng, a, beta * 2.0 / NCOL)
        gm = _su2_to_matrix(g, i, j)
        u_mu[mask] = gm @ u_mu[mask]


def sweep(
    rng: np.random.Generator,
    u: np.ndarray,
    beta: float,
    n_or: int = 1,
) -> None:
    """One quenched update sweep: heatbath plus n_or overrelaxation passes,
    checkerboarded per direction. Updates u in place and reunitarizes."""
    shape = u.shape[1:5]
    masks = [_parity_mask(shape, p) for p in (0, 1)]
    for do_or in [False] + [True] * n_or:
        for mu in range(NDIM):
            for mask in masks:
                # Re Tr[U_mu(x) staple_sum] reproduces the six plaquettes that
                # contain the link, so the weight uses the staple sum directly.
                st = staple_sum(u, mu)
                _apply_su2_update(rng, u[mu], st, mask, beta, do_or)
    for mu in range(NDIM):
        u[mu] = project_su3(u[mu])
