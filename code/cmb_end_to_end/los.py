"""Line-of-sight integration: sources S(k,tau) -> C_ell (Proposition B5)."""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn


def tau_grid(bg, rec):
    """Dense around recombination + coarse late-time grid for ISW/reionization."""
    i = int(np.argmax(rec.g_grid))
    tau_star = bg.tau_grid[i]
    dense = np.arange(max(tau_star - 180.0, 5.0), tau_star + 260.0, 1.4)
    late = np.geomspace(dense[-1] + 5.0, bg.tau0 - 1.0, 150)
    return np.unique(np.concatenate([dense, late]))


def coarse_k_grid(kmax=0.32):
    lo = np.geomspace(6e-6, 0.008, 52, endpoint=False)
    hi = np.arange(0.008, kmax, 0.0014)
    return np.concatenate([lo, hi])


def cl_from_sources(bg, k_coarse, taus, S, ells, n_s, A_s, k_pivot=0.05,
                    dk_fine=1.2e-4):
    """S: array (n_k_coarse, n_tau). Returns D_ell [muK^2] at requested ells."""
    tau0 = bg.tau0
    # log spacing at low k (plateau/ISW region), linear through the acoustic range
    k_fine = np.concatenate([
        np.geomspace(k_coarse[0], 0.004, 600, endpoint=False),
        np.arange(0.004, k_coarse[-1], dk_fine),
    ])
    # spline sources in k at each tau (S is smooth in k on the acoustic scale)
    S_fine = np.empty((k_fine.size, taus.size))
    for j in range(taus.size):
        S_fine[:, j] = CubicSpline(k_coarse, S[:, j])(k_fine)

    chi = tau0 - taus                      # comoving distance to each source slice
    w_tau = np.gradient(taus)              # trapezoid-like weights
    prim = A_s * (k_fine / k_pivot) ** (n_s - 1.0)
    dlnk = np.gradient(np.log(k_fine))
    t_cmb_uk = bg.p.t_cmb * 1e6

    x_max = k_fine[-1] * chi.max() + 10.0
    x_grid = np.arange(0.0, x_max, 0.02)

    D = np.zeros(len(ells))
    for i, ell in enumerate(ells):
        jl_tab = spherical_jn(int(ell), x_grid)
        x = np.outer(k_fine, chi)
        jl = np.interp(x, x_grid, jl_tab)
        theta_l = (S_fine * jl) @ w_tau     # Theta_ell(k)
        cl = 4.0 * np.pi * np.sum(dlnk * prim * theta_l**2)
        D[i] = ell * (ell + 1) / (2 * np.pi) * cl * t_cmb_uk**2
    return D
