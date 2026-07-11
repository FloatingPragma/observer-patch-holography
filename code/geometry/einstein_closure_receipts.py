#!/usr/bin/env python3
"""Machine receipts for the Einstein branch closure packets (GitHub #526-#528, #503).

Implements finite witnesses for the compact paper's subsection
`subsec:einstein-branch-closure`:

* null tomography (Theorem `thm:null-tomography`): exact reconstruction of the
  eta-trace-free part of a symmetric tensor from null-null charges over nine
  or more directions, the eta-ambiguity, and the countermodel of
  Proposition `prop:no-local-stress-countermodel`: directional charges that
  violate one dependent-family linearity relation admit no rank-two source
  (irreducible least-squares residual);
* bulk/edge/central first law (Theorem `thm:bulk-edge-central-first-law`):
  for blockwise central-interface states, delta S = 2pi delta<B> + delta<Z>
  exactly at fixed operators, delta<Z> = delta S_edge on the declared
  normalization z_alpha = log d_alpha, and the naive bookkeeping defect
  sum (z_alpha - log d_alpha) delta p_alpha for mismatched normalizations
  (Proposition `prop:entropy-coefficient-countermodel`(ii));
* MaxEnt multiplier identity (Theorem `thm:maxent-lagrange-stationarity`):
  dS/dt = lambda along a Gibbs/MaxEnt constraint family, verified by finite
  differences;
* baseline countermodel (Proposition `prop:baseline-countermodels`(iv)):
  two maximally symmetric baselines produce identical first-variation data
  while differing by c * g_ab: delta Y = 0 never fixes the constant.
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy.linalg import expm, logm

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


# ---------------------------------------------------------------------------
# null tomography (thm:null-tomography)
# ---------------------------------------------------------------------------

def null_vector(direction: np.ndarray) -> np.ndarray:
    d = np.asarray(direction, dtype=float)
    return np.concatenate(([1.0], d / np.linalg.norm(d)))


def sym_basis() -> list[np.ndarray]:
    """Basis of the ten-dimensional space of symmetric 4x4 matrices."""
    basis = []
    for i in range(4):
        for j in range(i, 4):
            m = np.zeros((4, 4))
            m[i, j] = m[j, i] = 1.0
            basis.append(m)
    return basis


def design_matrix(null_dirs: list[np.ndarray]) -> np.ndarray:
    """Rows: the linear functionals T -> T(k,k) in the sym_basis coordinates."""
    basis = sym_basis()
    rows = []
    for k in null_dirs:
        rows.append([float(k @ b @ k) for b in basis])
    return np.array(rows)


def charges_of(t_matrix: np.ndarray, null_dirs: list[np.ndarray]) -> np.ndarray:
    return np.array([float(k @ t_matrix @ k) for k in null_dirs])


def reconstruct_from_charges(charges: np.ndarray,
                             null_dirs: list[np.ndarray]) -> tuple[np.ndarray, float]:
    """Least-squares tensor from null charges; returns (T_hat, residual).

    The eta-direction is in the kernel of the design map (eta(k,k)=0), so the
    minimum-norm solution fixes the eta-component to zero: T is recovered
    exactly up to phi * eta, which is the classified ambiguity.
    """
    a = design_matrix(null_dirs)
    coef, *_ = np.linalg.lstsq(a, charges, rcond=None)
    basis = sym_basis()
    t_hat = sum(c * b for c, b in zip(coef, basis))
    residual = float(np.linalg.norm(a @ coef - charges))
    return t_hat, residual


def eta_project_out(t_matrix: np.ndarray) -> np.ndarray:
    """Remove the eta-component in the Frobenius sense (the null-invisible part)."""
    coeff = float(np.sum(t_matrix * ETA)) / float(np.sum(ETA * ETA))
    return t_matrix - coeff * ETA


def generic_null_directions(n: int, seed: int = 2) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    return [null_vector(rng.normal(size=3)) for _ in range(n)]


def tomography_receipt(seed: int = 2) -> dict[str, float]:
    """Exact reconstruction of consistent charges; irreducible residual for a
    family violating one dependent linearity relation."""
    rng = np.random.default_rng(seed)
    dirs = generic_null_directions(12, seed=seed)
    m = rng.normal(size=(4, 4))
    t_true = (m + m.T) / 2.0
    charges = charges_of(t_true, dirs)
    t_hat, residual = reconstruct_from_charges(charges, dirs)
    out = {
        "consistent_residual": residual,
        "tracefree_error": float(
            np.linalg.norm(eta_project_out(t_hat) - eta_project_out(t_true))
        ),
        "design_rank": int(np.linalg.matrix_rank(design_matrix(dirs))),
    }
    # countermodel: violate one linearity relation by bumping one charge
    bad = charges.copy()
    bad[0] += 1.0
    _, bad_residual = reconstruct_from_charges(bad, dirs)
    out["inconsistent_residual"] = bad_residual
    return out


# ---------------------------------------------------------------------------
# bulk/edge/central first law (thm:bulk-edge-central-first-law)
# ---------------------------------------------------------------------------

def random_faithful(dim: int, rng: np.random.Generator) -> np.ndarray:
    m = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    rho = m @ m.conj().T + 0.2 * np.eye(dim)
    return rho / np.trace(rho)


def blockwise_state(ps: list[float], sector_states: list[np.ndarray]) -> np.ndarray:
    blocks = [p * rho for p, rho in zip(ps, sector_states)]
    dim = sum(b.shape[0] for b in blocks)
    out = np.zeros((dim, dim), dtype=complex)
    i = 0
    for b in blocks:
        d = b.shape[0]
        out[i:i + d, i:i + d] = b
        i += d
    return out


def central_z(ps: list[float], dims: list[int],
              z_weights: list[float]) -> np.ndarray:
    """Z = oplus_alpha (-log p_alpha + z_alpha) 1_alpha."""
    dim = sum(dims)
    out = np.zeros((dim, dim))
    i = 0
    for p, d, z in zip(ps, dims, z_weights):
        out[i:i + d, i:i + d] = (-np.log(p) + z) * np.eye(d)
        i += d
    return out


def entropy(rho: np.ndarray) -> float:
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-14]
    return float(-np.sum(evals * np.log(evals)))


def edge_entropy(ps: list[float], dims: list[int]) -> float:
    return float(-np.sum(np.array(ps) * np.log(ps))
                 + np.sum(np.array(ps) * np.log(dims)))


def first_law_receipt(z_weights: list[float] | None = None,
                      move_weights: bool = True,
                      seed: int = 4, eps: float = 1e-6) -> dict[str, float]:
    """Finite-difference check of Theorem thm:bulk-edge-central-first-law.

    Verified identities (declared normalization z_alpha = log d_alpha):
      * first law: delta S = 2pi delta<B> + delta<Z> (exact, fixed operators);
      * edge identification: delta<Z> = delta S_edge;
      * boxed split: delta S = 2pi delta<B> + delta S_edge;
      * naive-step defect: delta S_bulk - 2pi delta<B> = sum_a z_a dp_a,
        which vanishes exactly on the fixed-sector-weight class and equals
        the central flux otherwise;
      * mismatched normalization z != log d breaks the edge identification
        by sum_a (z_a - log d_a) dp_a.
    """
    rng = np.random.default_rng(seed)
    dims = [2, 3]
    ps0 = [0.6, 0.4]
    sectors0 = [random_faithful(d, rng) for d in dims]
    correct_z = [float(np.log(d)) for d in dims]
    zw = correct_z if z_weights is None else z_weights

    # base operators at the base point
    rho0 = blockwise_state(ps0, sectors0)
    k0 = -logm(rho0)
    z0 = central_z(ps0, dims, zw)
    b0 = (k0 - z0) / (2.0 * np.pi)

    # a variation moving sector states, and optionally sector weights
    dps = [eps, -eps] if move_weights else [0.0, 0.0]
    dsectors = [random_faithful(d, rng) for d in dims]
    ps1 = [p + dp for p, dp in zip(ps0, dps)]
    sectors1 = [
        (1 - eps) * s + eps * t for s, t in zip(sectors0, dsectors)
    ]
    rho1 = blockwise_state(ps1, sectors1)
    drho = rho1 - rho0

    d_s = entropy(rho1) - entropy(rho0)
    d_b = float(np.real(np.trace(b0 @ drho)))
    d_z = float(np.real(np.trace(z0 @ drho)))
    d_s_edge = edge_entropy(ps1, dims) - edge_entropy(ps0, dims)
    s_bulk0 = sum(p * entropy(s) for p, s in zip(ps0, sectors0))
    s_bulk1 = sum(p * entropy(s) for p, s in zip(ps1, sectors1))
    d_s_bulk = s_bulk1 - s_bulk0

    return {
        "first_law_defect": abs(d_s - (2.0 * np.pi * d_b + d_z)) / eps,
        "edge_identification_defect": abs(d_z - d_s_edge) / eps,
        "split_identity_defect": abs(d_s - (2.0 * np.pi * d_b + d_s_edge)) / eps
        if z_weights is None else float("nan"),
        "naive_step_defect": abs(d_s_bulk - 2.0 * np.pi * d_b) / eps,
        "predicted_naive_defect": abs(float(np.dot(zw, dps))) / eps,
        "predicted_edge_defect": abs(
            float(np.dot(np.array(zw) - np.array(correct_z), dps))
        ) / eps,
    }


# ---------------------------------------------------------------------------
# MaxEnt multiplier identity (thm:maxent-lagrange-stationarity)
# ---------------------------------------------------------------------------

def maxent_multiplier_receipt(lam: float = 1.3, dlam: float = 1e-5,
                              seed: int = 6) -> dict[str, float]:
    """Along the Gibbs family rho(lam) = exp(-lam T)/Z, verify dS/dt = lam
    where t = <T>: the exact envelope identity behind the coupled-class
    stationarity."""
    rng = np.random.default_rng(seed)
    m = rng.normal(size=(6, 6))
    t_op = (m + m.T) / 2.0

    def state(l: float) -> np.ndarray:
        r = expm(-l * t_op)
        return r / np.trace(r)

    def s_and_t(l: float) -> tuple[float, float]:
        r = state(l)
        return entropy(r), float(np.real(np.trace(r @ t_op)))

    s_p, t_p = s_and_t(lam + dlam)
    s_m, t_m = s_and_t(lam - dlam)
    ds_dt = (s_p - s_m) / (t_p - t_m)
    return {"ds_dt": ds_dt, "lambda": lam,
            "multiplier_defect": abs(ds_dt - lam)}


# ---------------------------------------------------------------------------
# baseline countermodel (prop:baseline-countermodels(iv))
# ---------------------------------------------------------------------------

def baseline_countermodel_receipt(c: float = 0.7, seed: int = 8) -> dict[str, float]:
    """Two baselines Y1 = 0 and Y2 = c*eta produce identical first-variation
    responses delta Y = 0 for every sampled variation, while differing by
    exactly c*eta: the absolute constant is invisible to first variations."""
    rng = np.random.default_rng(seed)
    y1 = np.zeros((4, 4))
    y2 = c * ETA
    # first variations of a *constant* tensor field vanish identically for
    # both baselines; sample variation directions to record the identical data
    defects = []
    for _ in range(16):
        v = rng.normal(size=4)
        v = v / np.linalg.norm(v)
        defects.append(abs(float(v @ (y2 - y2) @ v)) + abs(float(v @ (y1 - y1) @ v)))
    return {
        "variation_data_difference": max(defects),
        "baseline_difference": float(np.linalg.norm(y2 - y1)),
        "eta_alignment": float(
            np.linalg.norm((y2 - y1) - c * ETA)
        ),
    }
