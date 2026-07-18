"""Source-side radial-lift reference implementation for OPH issue #330.

The module implements the mathematical receipts proved in
``RADIAL_LIFT_THEOREMS_330.tex``:

* exact shell/window projection from a homogeneous isotropic 3-D curvature field;
* exact Mellin integral for the thin-shell power-law branch;
* source-derived amplitude conversion (never an angular-spectrum fit);
* a rigorous finite-window stability bound in a weighted Bessel Hilbert space;
* finite-basis SVD/null-space and prior-selected continuation diagnostics;
* the explicit dilation-intertwiner residual needed to derive the radial power law;
* fail-closed receipt construction and the TT/TE/EE promotion firewall.

This module does not run a Boltzmann solver and contains no observational target.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import null_space
from scipy.special import digamma, gammaln, spherical_jn


class RadialLiftInputError(ValueError):
    """Raised when an input violates the theorem contract."""


@dataclass(frozen=True)
class PrimordialAmplitude:
    """Exact source-side thin-shell amplitude conversion."""

    A_q: float
    theta: float
    A_zeta: float
    Z_q: float
    R_star: float
    k_pivot: float
    kR_pivot: float
    conversion_factor: float


@dataclass(frozen=True)
class WindowBound:
    """Certified finite-window deviation from the source shell."""

    ell: int
    theta: float
    I_ell_theta: float
    J_ell_theta: float
    eta: float
    shell_norm: float
    absolute_cl_bound: float
    relative_to_shell_cl_bound: float


@dataclass(frozen=True)
class NullSpaceReport:
    shape: tuple[int, int]
    rank: int
    nullity: int
    singular_values: list[float]
    effective_threshold: float
    condition_number_nonzero: float
    null_basis: list[list[float]]


@dataclass(frozen=True)
class PriorContinuation:
    """Minimum-Q exact continuation and its resolution operator."""

    p: list[float]
    residual: list[float]
    residual_norm: float
    objective: float
    resolution: list[list[float]]
    null_projector: list[list[float]]
    effective_rank: int


@dataclass(frozen=True)
class DilationReceipt:
    """Finite safe-band test of the physical dilation-intertwiner law."""

    theta: float
    scale_ratios: list[float]
    max_absolute_log_residual: float
    rms_log_residual: float
    passed: bool
    tolerance: float
    evaluated_pairs: int


# ---------------------------------------------------------------------------
# Exact shell formulas
# ---------------------------------------------------------------------------


def mellin_spherical_bessel_square(ell: int, theta: float) -> float:
    r"""Return ``I_l(theta)=∫ dln(x) x^{-theta} j_l(x)^2`` exactly.

    The formula is

    .. math::

        I_\ell(\theta)=\frac{\sqrt\pi}{4}
        \frac{\Gamma(1+\theta/2)}{\Gamma(3/2+\theta/2)}
        \frac{\Gamma(\ell-\theta/2)}
             {\Gamma(\ell+2+\theta/2)}.

    Absolute convergence holds for ``-2 < theta < 2*ell``.  The OPH retained
    band starts at ``ell=2``, so a common sufficient interval is ``(-2, 4)``.
    """

    ell = _integer_at_least(ell, 0, "ell")
    theta = _finite(theta, "theta")
    if not (-2.0 < theta < 2.0 * ell):
        raise RadialLiftInputError(
            f"Mellin integral requires -2 < theta < 2*ell; got ell={ell}, theta={theta}"
        )
    log_value = (
        0.5 * math.log(math.pi)
        - math.log(4.0)
        + gammaln(1.0 + 0.5 * theta)
        - gammaln(1.5 + 0.5 * theta)
        + gammaln(ell - 0.5 * theta)
        - gammaln(ell + 2.0 + 0.5 * theta)
    )
    value = math.exp(log_value)
    if not (math.isfinite(value) and value > 0.0):
        raise RadialLiftInputError("Mellin integral evaluated to a nonpositive/nonfinite value")
    return value


def derivative_mellin_norm(ell: int, theta: float) -> float:
    r"""Return ``J_l(theta)=∫ dln(x) x^{2-theta} j_l'(x)^2`` exactly.

    Integration by parts with the spherical-Bessel equation gives

    .. math::

        J_\ell(\theta)=I_\ell(\theta-2)
        -\left[\ell(\ell+1)-\frac{\theta(\theta+1)}2\right]I_\ell(\theta).

    The derivative norm is finite for ``0 < theta < 2*ell``.
    """

    ell = _integer_at_least(ell, 1, "ell")
    theta = _finite(theta, "theta")
    if not (0.0 < theta < 2.0 * ell):
        raise RadialLiftInputError("derivative Mellin norm requires 0 < theta < 2*ell")
    value = mellin_spherical_bessel_square(ell, theta - 2.0) - (
        ell * (ell + 1.0) - 0.5 * theta * (theta + 1.0)
    ) * mellin_spherical_bessel_square(ell, theta)
    # Roundoff can only produce a tiny negative number near a boundary.
    if value < 0.0 and abs(value) < 1e-13:
        value = 0.0
    if not (math.isfinite(value) and value >= 0.0):
        raise RadialLiftInputError("derived derivative norm is negative or nonfinite")
    return value


def screen_gamma_ratio_cl(ell: ArrayLike, A_q: float, theta: float) -> NDArray[np.float64] | float:
    """Return the exact gamma-ratio screen spectrum ``C_l^q``."""

    A_q = _positive(A_q, "A_q")
    theta = _finite(theta, "theta")
    ell_arr = np.asarray(ell, dtype=float)
    if np.any(~np.isfinite(ell_arr)) or np.any(ell_arr < 2.0):
        raise RadialLiftInputError("ell must be finite and at least 2")
    if not (-2.0 < theta < 4.0):
        raise RadialLiftInputError("the common retained-band interval is -2 < theta < 4")
    log_values = (
        math.log(A_q)
        + gammaln(ell_arr - 0.5 * theta)
        - gammaln(ell_arr + 2.0 + 0.5 * theta)
    )
    values = np.exp(log_values)
    if np.isscalar(ell):
        return float(values)
    return np.asarray(values, dtype=float)


def primordial_amplitude_from_screen(
    A_q: float,
    theta: float,
    *,
    Z_q: float,
    R_star: float,
    k_pivot: float,
) -> PrimordialAmplitude:
    r"""Convert source amplitude ``A_q`` to ``A_zeta`` without fitting.

    For

    .. math::

        \Delta_\zeta^2(k)=A_\zeta(k/k_*)^{-\theta},\qquad
        q(\hat n)=Z_q\,\zeta(R_*\hat n),

    the exact result is

    .. math::

        A_\zeta=\frac{A_q}{\pi^{3/2} Z_q^2(k_*R_*)^\theta}
        \frac{\Gamma(3/2+\theta/2)}{\Gamma(1+\theta/2)}.
    """

    A_q = _positive(A_q, "A_q")
    theta = _finite(theta, "theta")
    Z_q = _positive(Z_q, "Z_q")
    R_star = _positive(R_star, "R_star")
    k_pivot = _positive(k_pivot, "k_pivot")
    if not (-2.0 < theta < 4.0):
        raise RadialLiftInputError("thin-shell retained-band formula requires -2 < theta < 4")
    kR = k_pivot * R_star
    log_factor = (
        -1.5 * math.log(math.pi)
        - 2.0 * math.log(Z_q)
        - theta * math.log(kR)
        + gammaln(1.5 + 0.5 * theta)
        - gammaln(1.0 + 0.5 * theta)
    )
    factor = math.exp(log_factor)
    return PrimordialAmplitude(
        A_q=A_q,
        theta=theta,
        A_zeta=A_q * factor,
        Z_q=Z_q,
        R_star=R_star,
        k_pivot=k_pivot,
        kR_pivot=kR,
        conversion_factor=factor,
    )


def screen_amplitude_from_primordial(
    A_zeta: float,
    theta: float,
    *,
    Z_q: float,
    R_star: float,
    k_pivot: float,
) -> float:
    """Inverse of :func:`primordial_amplitude_from_screen`."""

    A_zeta = _positive(A_zeta, "A_zeta")
    theta = _finite(theta, "theta")
    Z_q = _positive(Z_q, "Z_q")
    R_star = _positive(R_star, "R_star")
    k_pivot = _positive(k_pivot, "k_pivot")
    if not (-2.0 < theta < 4.0):
        raise RadialLiftInputError("thin-shell retained-band formula requires -2 < theta < 4")
    log_factor = (
        1.5 * math.log(math.pi)
        + 2.0 * math.log(Z_q)
        + theta * math.log(k_pivot * R_star)
        + gammaln(1.0 + 0.5 * theta)
        - gammaln(1.5 + 0.5 * theta)
    )
    return A_zeta * math.exp(log_factor)


def thin_shell_powerlaw_cl(
    ell: ArrayLike,
    A_zeta: float,
    theta: float,
    *,
    Z_q: float,
    R_star: float,
    k_pivot: float,
) -> NDArray[np.float64] | float:
    """Exact thin-shell angular spectrum of the source power law."""

    A_q = screen_amplitude_from_primordial(
        A_zeta, theta, Z_q=Z_q, R_star=R_star, k_pivot=k_pivot
    )
    return screen_gamma_ratio_cl(ell, A_q, theta)


def primordial_amplitude_log_sensitivity(theta: float, kR_pivot: float) -> float:
    r"""Return ``∂_theta ln(A_zeta/A_q)`` at fixed ``Z_q`` and ``k_*R_*``."""

    theta = _finite(theta, "theta")
    kR_pivot = _positive(kR_pivot, "kR_pivot")
    if not (-2.0 < theta < 4.0):
        raise RadialLiftInputError("theta must lie in (-2, 4)")
    return float(
        -math.log(kR_pivot)
        + 0.5 * (digamma(1.5 + 0.5 * theta) - digamma(1.0 + 0.5 * theta))
    )


# ---------------------------------------------------------------------------
# Exact finite-window forward kernel and stability bound
# ---------------------------------------------------------------------------


def normalized_radial_window(radii: ArrayLike, weights: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Validate and normalize a positive discrete radial window."""

    r = np.asarray(radii, dtype=float)
    w = np.asarray(weights, dtype=float)
    if r.ndim != 1 or w.shape != r.shape or r.size == 0:
        raise RadialLiftInputError("radii and weights must be nonempty matching vectors")
    if np.any(~np.isfinite(r)) or np.any(r <= 0.0):
        raise RadialLiftInputError("all radii must be positive and finite")
    if np.any(~np.isfinite(w)) or np.any(w < 0.0):
        raise RadialLiftInputError("window weights must be finite and nonnegative")
    total = float(np.sum(w))
    if total <= 0.0:
        raise RadialLiftInputError("window weights must have positive total")
    return r, w / total


def window_transfer(
    ell: int,
    k: ArrayLike,
    radii: ArrayLike,
    weights: ArrayLike,
) -> NDArray[np.float64]:
    r"""Return ``Psi_l(k)=sum_i w_i j_l(k r_i)`` for a discrete window."""

    ell = _integer_at_least(ell, 0, "ell")
    k_arr = np.asarray(k, dtype=float)
    if k_arr.ndim != 1 or k_arr.size == 0 or np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0):
        raise RadialLiftInputError("k must be a nonempty positive finite vector")
    r, w = normalized_radial_window(radii, weights)
    values = spherical_jn(ell, np.outer(k_arr, r)) @ w
    return np.asarray(values, dtype=float)


def window_powerlaw_cl_quadrature(
    ell: int,
    A_zeta: float,
    theta: float,
    *,
    Z_q: float,
    k_pivot: float,
    k: ArrayLike,
    dlnk_weights: ArrayLike,
    radii: ArrayLike,
    radial_weights: ArrayLike,
) -> float:
    """Finite quadrature of the exact radial-window forward projection."""

    ell = _integer_at_least(ell, 0, "ell")
    A_zeta = _positive(A_zeta, "A_zeta")
    theta = _finite(theta, "theta")
    Z_q = _positive(Z_q, "Z_q")
    k_pivot = _positive(k_pivot, "k_pivot")
    k_arr = np.asarray(k, dtype=float)
    wk = np.asarray(dlnk_weights, dtype=float)
    if k_arr.ndim != 1 or wk.shape != k_arr.shape or k_arr.size == 0:
        raise RadialLiftInputError("k and dlnk_weights must be matching vectors")
    if np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0):
        raise RadialLiftInputError("k must be positive and finite")
    if np.any(~np.isfinite(wk)) or np.any(wk <= 0.0):
        raise RadialLiftInputError("dlnk_weights must be positive and finite")
    psi = window_transfer(ell, k_arr, radii, radial_weights)
    delta = A_zeta * (k_arr / k_pivot) ** (-theta)
    value = 4.0 * math.pi * Z_q**2 * float(np.dot(wk, delta * psi**2))
    return _positive(value, "projected C_l")


def finite_window_stability_bound(
    ell: int,
    theta: float,
    *,
    A_zeta: float,
    Z_q: float,
    k_pivot: float,
    R_star: float,
    radii: ArrayLike,
    radial_weights: ArrayLike,
) -> WindowBound:
    r"""Certify the finite-window deviation from a thin shell.

    In the Hilbert space with norm

    .. math::

        \|f\|_\theta^2=\int_0^\infty d\ln k\,k^{-\theta}|f(k)|^2,

    set ``f_r(k)=j_l(kr)`` and ``f_W=sum_i w_i f_{r_i}``.  The theorem proves

    .. math::

        \|f_W-f_{R_*}\|_\theta\le
        \eta=\frac{2\sqrt{J_\ell(\theta)}}{\theta}
        \sum_iw_i|r_i^{\theta/2}-R_*^{\theta/2}|,

    and then bounds the difference of squared norms.  This retains the Bessel
    ultraviolet decay and is integrable; it replaces the unsafe pointwise Taylor
    bound that can lose that decay after a supremum is taken.
    """

    ell = _integer_at_least(ell, 1, "ell")
    theta = _finite(theta, "theta")
    if not (0.0 < theta < 2.0 * ell):
        raise RadialLiftInputError("window stability theorem requires 0 < theta < 2*ell")
    A_zeta = _positive(A_zeta, "A_zeta")
    Z_q = _positive(Z_q, "Z_q")
    k_pivot = _positive(k_pivot, "k_pivot")
    R_star = _positive(R_star, "R_star")
    r, w = normalized_radial_window(radii, radial_weights)
    I = mellin_spherical_bessel_square(ell, theta)
    J = derivative_mellin_norm(ell, theta)
    a = 0.5 * theta
    eta = (2.0 * math.sqrt(J) / theta) * float(np.dot(w, np.abs(r**a - R_star**a)))
    shell_norm = R_star**a * math.sqrt(I)
    prefactor = 4.0 * math.pi * Z_q**2 * A_zeta * k_pivot**theta
    abs_bound = prefactor * eta * (2.0 * shell_norm + eta)
    shell_cl = prefactor * shell_norm**2
    rel_bound = abs_bound / shell_cl
    return WindowBound(
        ell=ell,
        theta=theta,
        I_ell_theta=I,
        J_ell_theta=J,
        eta=eta,
        shell_norm=shell_norm,
        absolute_cl_bound=abs_bound,
        relative_to_shell_cl_bound=rel_bound,
    )


# ---------------------------------------------------------------------------
# Finite radial operator, null-space and prior-selected continuation
# ---------------------------------------------------------------------------


def radial_projection_matrix(
    ell: ArrayLike,
    k: ArrayLike,
    dlnk_weights: ArrayLike,
    *,
    Z_q: float,
    radii: ArrayLike,
    radial_weights: ArrayLike,
) -> NDArray[np.float64]:
    r"""Build ``A_{ell,j}=4*pi*Z_q^2*w_j*|Psi_l(k_j)|^2``."""

    ell_arr = np.asarray(ell, dtype=int)
    k_arr = np.asarray(k, dtype=float)
    wk = np.asarray(dlnk_weights, dtype=float)
    Z_q = _positive(Z_q, "Z_q")
    if ell_arr.ndim != 1 or ell_arr.size == 0 or np.any(ell_arr < 0):
        raise RadialLiftInputError("ell must be a nonempty vector of nonnegative integers")
    if k_arr.ndim != 1 or wk.shape != k_arr.shape or k_arr.size == 0:
        raise RadialLiftInputError("k and dlnk_weights must be matching vectors")
    if np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0):
        raise RadialLiftInputError("k must be positive and finite")
    if np.any(~np.isfinite(wk)) or np.any(wk <= 0.0):
        raise RadialLiftInputError("dlnk_weights must be positive and finite")
    A = np.empty((ell_arr.size, k_arr.size), dtype=float)
    for i, lval in enumerate(ell_arr):
        psi = window_transfer(int(lval), k_arr, radii, radial_weights)
        A[i, :] = 4.0 * math.pi * Z_q**2 * wk * psi**2
    return A


def radial_null_space_report(matrix: ArrayLike, *, rtol: float = 1e-12) -> NullSpaceReport:
    """Return the complete finite-basis right-null report."""

    A = np.asarray(matrix, dtype=float)
    if A.ndim != 2 or A.size == 0 or np.any(~np.isfinite(A)):
        raise RadialLiftInputError("matrix must be a finite nonempty 2-D array")
    rtol = _positive(rtol, "rtol")
    singular = np.linalg.svd(A, compute_uv=False)
    threshold = rtol * (float(singular[0]) if singular.size else 0.0)
    rank = int(np.sum(singular > threshold))
    basis = null_space(A, rcond=rtol)
    retained = singular[singular > threshold]
    condition = float(retained[0] / retained[-1]) if retained.size else math.inf
    return NullSpaceReport(
        shape=(int(A.shape[0]), int(A.shape[1])),
        rank=rank,
        nullity=int(A.shape[1] - rank),
        singular_values=[float(x) for x in singular],
        effective_threshold=float(threshold),
        condition_number_nonzero=condition,
        null_basis=basis.T.tolist(),
    )


def minimum_prior_continuation(
    matrix: ArrayLike,
    screen_cl: ArrayLike,
    *,
    prior_center: ArrayLike,
    prior_precision: ArrayLike,
    rtol: float = 1e-12,
) -> PriorContinuation:
    r"""Return the exact minimum-prior continuation under ``A p=C``.

    It minimizes

    .. math::

        \frac12(p-p_0)^TQ(p-p_0)\quad\text{subject to}\quad Ap=C,

    and therefore returns

    .. math::

        p_*=p_0+Q^{-1}A^T(AQ^{-1}A^T)^+(C-Ap_0).

    This is a prior-selected representative, not a model-free inverse.  Its
    resolution and null projectors are returned explicitly.
    """

    A = np.asarray(matrix, dtype=float)
    C = np.asarray(screen_cl, dtype=float)
    p0 = np.asarray(prior_center, dtype=float)
    Q = np.asarray(prior_precision, dtype=float)
    if A.ndim != 2 or C.shape != (A.shape[0],) or p0.shape != (A.shape[1],):
        raise RadialLiftInputError("matrix, screen_cl and prior_center dimensions do not match")
    if Q.shape != (A.shape[1], A.shape[1]) or np.any(~np.isfinite(Q)):
        raise RadialLiftInputError("prior_precision must be a finite square matrix")
    if not np.allclose(Q, Q.T, rtol=0.0, atol=1e-12):
        raise RadialLiftInputError("prior_precision must be symmetric")
    q_eigs = np.linalg.eigvalsh(Q)
    if float(np.min(q_eigs)) <= 0.0:
        raise RadialLiftInputError("prior_precision must be positive definite")
    Qinv = np.linalg.inv(Q)
    G = A @ Qinv @ A.T
    Gplus = np.linalg.pinv(G, rcond=rtol, hermitian=True)
    resolution = Qinv @ A.T @ Gplus @ A
    p = p0 + Qinv @ A.T @ Gplus @ (C - A @ p0)
    residual = C - A @ p
    scale = max(float(np.linalg.norm(C)), 1.0)
    if float(np.linalg.norm(residual)) > 100.0 * rtol * scale:
        raise RadialLiftInputError("exact constraints are inconsistent at the requested tolerance")
    d = p - p0
    objective = 0.5 * float(d @ Q @ d)
    report = radial_null_space_report(A, rtol=rtol)
    I = np.eye(A.shape[1], dtype=float)
    return PriorContinuation(
        p=p.tolist(),
        residual=residual.tolist(),
        residual_norm=float(np.linalg.norm(residual)),
        objective=objective,
        resolution=resolution.tolist(),
        null_projector=(I - resolution).tolist(),
        effective_rank=report.rank,
    )


def source_powerlaw(k: ArrayLike, A_zeta: float, theta: float, k_pivot: float) -> NDArray[np.float64]:
    """Return the one-dimensional source-derived radial family."""

    k_arr = np.asarray(k, dtype=float)
    if k_arr.ndim != 1 or k_arr.size == 0 or np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0):
        raise RadialLiftInputError("k must be a nonempty positive finite vector")
    return _positive(A_zeta, "A_zeta") * (k_arr / _positive(k_pivot, "k_pivot")) ** (-_finite(theta, "theta"))


def forward_residual(matrix: ArrayLike, spectrum: ArrayLike, screen_cl: ArrayLike) -> Mapping[str, Any]:
    """Return a non-fitting forward residual receipt."""

    A = np.asarray(matrix, dtype=float)
    p = np.asarray(spectrum, dtype=float)
    C = np.asarray(screen_cl, dtype=float)
    if A.ndim != 2 or p.shape != (A.shape[1],) or C.shape != (A.shape[0],):
        raise RadialLiftInputError("matrix, spectrum and screen_cl dimensions do not match")
    predicted = A @ p
    residual = C - predicted
    denom = max(float(np.linalg.norm(C)), 1e-300)
    return {
        "predicted": predicted.tolist(),
        "residual": residual.tolist(),
        "absolute_l2_residual": float(np.linalg.norm(residual)),
        "relative_l2_residual": float(np.linalg.norm(residual) / denom),
    }


# ---------------------------------------------------------------------------
# Physical dilation-intertwiner receipt
# ---------------------------------------------------------------------------


def dilation_intertwiner_receipt(
    k: ArrayLike,
    delta_zeta_sq: ArrayLike,
    theta: float,
    *,
    scale_ratios: Sequence[float],
    tolerance: float,
) -> DilationReceipt:
    r"""Test ``Delta^2(b k)=b^{-theta} Delta^2(k)`` on a finite safe band.

    The theorem-level producer is the physical operator relation

    ``D_s^{-1} C_zeta D_s = exp(-theta*s) C_zeta``

    on the common ``dln(k)`` mode basis.  This finite function tests its diagonal
    consequence after interpolation; it does not substitute for the source DAG
    that proves the operator intertwiner.
    """

    k_arr = np.asarray(k, dtype=float)
    p = np.asarray(delta_zeta_sq, dtype=float)
    theta = _finite(theta, "theta")
    tolerance = _positive(tolerance, "tolerance")
    if k_arr.ndim != 1 or p.shape != k_arr.shape or k_arr.size < 3:
        raise RadialLiftInputError("k and delta_zeta_sq must be matching vectors with at least 3 points")
    if np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0) or np.any(np.diff(k_arr) <= 0.0):
        raise RadialLiftInputError("k must be strictly increasing, positive and finite")
    if np.any(~np.isfinite(p)) or np.any(p <= 0.0):
        raise RadialLiftInputError("delta_zeta_sq must be positive and finite")
    logk = np.log(k_arr)
    logp = np.log(p)
    residuals: list[float] = []
    ratios: list[float] = []
    for b_raw in scale_ratios:
        b = _positive(b_raw, "scale ratio")
        if math.isclose(b, 1.0):
            continue
        shifted = logk + math.log(b)
        mask = (shifted >= logk[0]) & (shifted <= logk[-1])
        if not np.any(mask):
            continue
        interpolated = np.interp(shifted[mask], logk, logp)
        res = interpolated - logp[mask] + theta * math.log(b)
        residuals.extend(float(x) for x in res)
        ratios.append(b)
    if not residuals:
        raise RadialLiftInputError("no scale-ratio pair lies inside the safe k band")
    arr = np.asarray(residuals, dtype=float)
    max_res = float(np.max(np.abs(arr)))
    rms = float(np.sqrt(np.mean(arr**2)))
    return DilationReceipt(
        theta=theta,
        scale_ratios=ratios,
        max_absolute_log_residual=max_res,
        rms_log_residual=rms,
        passed=bool(max_res <= tolerance),
        tolerance=tolerance,
        evaluated_pairs=int(arr.size),
    )


def approximate_dilation_shape_bound(
    k: ArrayLike,
    epsilon_log_slope: ArrayLike,
    *,
    k_pivot: float,
) -> NDArray[np.float64]:
    r"""Integrate a bound on ``|d ln Delta^2/d ln k + theta|``.

    Returns the cumulative upper bound

    .. math::

        \left|\ln\frac{\Delta^2(k)}{A_\zeta(k/k_*)^{-\theta}}\right|
        \le\left|\int_{\ln k_*}^{\ln k}\epsilon(e^t)dt\right|.
    """

    k_arr = np.asarray(k, dtype=float)
    eps = np.asarray(epsilon_log_slope, dtype=float)
    k_pivot = _positive(k_pivot, "k_pivot")
    if k_arr.ndim != 1 or eps.shape != k_arr.shape or k_arr.size < 2:
        raise RadialLiftInputError("k and epsilon_log_slope must be matching vectors")
    if np.any(~np.isfinite(k_arr)) or np.any(k_arr <= 0.0) or np.any(np.diff(k_arr) <= 0.0):
        raise RadialLiftInputError("k must be strictly increasing, positive and finite")
    if np.any(~np.isfinite(eps)) or np.any(eps < 0.0):
        raise RadialLiftInputError("epsilon_log_slope must be finite and nonnegative")
    logk = np.log(k_arr)
    # Include the pivot by interpolation in the cumulative trapezoid grid.
    if not (k_arr[0] <= k_pivot <= k_arr[-1]):
        raise RadialLiftInputError("k_pivot must lie inside the k grid")
    pivot_log = math.log(k_pivot)
    grid = np.unique(np.concatenate([logk, np.array([pivot_log])]))
    eps_grid = np.interp(grid, logk, eps)
    increments = 0.5 * (eps_grid[:-1] + eps_grid[1:]) * np.diff(grid)
    cumulative = np.concatenate([[0.0], np.cumsum(increments)])
    pivot_index = int(np.where(np.isclose(grid, pivot_log, rtol=0.0, atol=1e-14))[0][0])
    bound_grid = np.abs(cumulative - cumulative[pivot_index])
    return np.interp(logk, grid, bound_grid)


# ---------------------------------------------------------------------------
# Fail-closed receipts
# ---------------------------------------------------------------------------


RADIAL_RECEIPTS = {
    "SCR330_SOURCE_SHELL_EMBEDDING_RECEIPT",
    "SCR330_PHYSICAL_MODE_BASIS_RECEIPT",
    "SCR330_RADIAL_DILATION_INTERTWINER_RECEIPT",
    "SCR330_THIN_SHELL_MELLIN_LIFT_RECEIPT",
    "SCR330_FINITE_WINDOW_KERNEL_RECEIPT",
    "SCR330_RADIAL_NULL_REPORT",
    "SCR330_RADIAL_FORWARD_RESIDUAL_RECEIPT",
    "SCR330_RADIAL_TOMOGRAPHY_RECEIPT",
    "SCR330_RADIAL_PROMOTION_RECEIPT",
    "SCR330_TRANSFER_FIREWALL_RECEIPT",
}


def build_radial_receipt(
    *,
    receipt: str,
    passed: bool,
    claim_tier: str,
    source_dag: Mapping[str, Any],
    blockers: Iterable[str] = (),
    payload: Mapping[str, Any] | None = None,
    physical_tt_te_ee_claim: bool = False,
) -> dict[str, Any]:
    """Build a fail-closed radial-lift receipt."""

    if receipt not in RADIAL_RECEIPTS:
        raise RadialLiftInputError(f"unknown radial receipt: {receipt}")
    if claim_tier not in {"E0", "E1", "E2", "E3", "E4", "E5"}:
        raise RadialLiftInputError("claim_tier must be E0 through E5")
    blocker_set = set(str(x) for x in blockers)
    no_target_ancestor = not _dag_has_blacklisted_ancestor(source_dag)
    if passed and not no_target_ancestor:
        passed = False
        blocker_set.add("measurement_fit_or_likelihood_ancestor")
    if physical_tt_te_ee_claim and claim_tier != "E5":
        passed = False
        blocker_set.add("tt_te_ee_claim_before_E5")
    if receipt == "SCR330_RADIAL_PROMOTION_RECEIPT" and claim_tier != "E4":
        passed = False
        blocker_set.add("radial_primordial_promotion_requires_E4")
    canonical = json.dumps(source_dag, sort_keys=True, separators=(",", ":")).encode("utf-8")
    result: dict[str, Any] = {
        "schema_version": "scr330-radial-v2",
        "receipt": receipt,
        "passed": bool(passed),
        "claim_tier": claim_tier,
        "source_dag_hash": "sha256:" + hashlib.sha256(canonical).hexdigest(),
        "no_measurement_fit_likelihood_ancestor": bool(no_target_ancestor),
        "physical_tt_te_ee_claim": bool(physical_tt_te_ee_claim),
        "blockers": sorted(blocker_set),
    }
    if payload:
        result["payload"] = dict(payload)
    return result


def write_json(path: str | Path, value: Mapping[str, Any] | Sequence[Any]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def summary_payload(A_q: float, theta: float, *, Z_q: float, R_star: float, k_pivot: float) -> dict[str, Any]:
    """Return a compact source-side summary useful in simulator receipts."""

    amp = primordial_amplitude_from_screen(
        A_q, theta, Z_q=Z_q, R_star=R_star, k_pivot=k_pivot
    )
    return {"amplitude": asdict(amp), "A_zeta_over_A_q": amp.conversion_factor}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _dag_has_blacklisted_ancestor(dag: Mapping[str, Any]) -> bool:
    nodes = dag.get("nodes")
    if not isinstance(nodes, list):
        return True
    for node in nodes:
        if not isinstance(node, Mapping):
            return True
        if any(
            bool(node.get(key))
            for key in (
                "measurement",
                "fit",
                "likelihood",
                "posterior",
                "planck_shape",
                "target_calibrated_proxy",
                "metadata_unknown",
            )
        ):
            return True
    return False


def _finite(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RadialLiftInputError(f"{name} must be finite")
    return result


def _positive(value: float, name: str) -> float:
    result = _finite(value, name)
    if result <= 0.0:
        raise RadialLiftInputError(f"{name} must be positive")
    return result


def _integer_at_least(value: int, minimum: int, name: str) -> int:
    if isinstance(value, bool):
        raise RadialLiftInputError(f"{name} must be an integer")
    result = int(value)
    if result != value or result < minimum:
        raise RadialLiftInputError(f"{name} must be an integer at least {minimum}")
    return result


if __name__ == "__main__":
    theta = 1.630968209403959 / 48.0
    payload = summary_payload(1.0, theta, Z_q=1.0, R_star=1.0, k_pivot=1.0)
    print(json.dumps(payload, indent=2, sort_keys=True))
