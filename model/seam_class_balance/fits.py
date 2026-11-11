"""Fits of the measured class balance: ln <A> = ln M0 + p ln S - c S, and the sector law of the untwisted link."""
import math

import numpy as np
from scipy import special

from .exact import S_U


def wls(X, y, e):
    W = 1.0 / np.asarray(e) ** 2
    A = X.T @ (X * W[:, None])
    cov = np.linalg.inv(A)
    par = cov @ (X.T @ (W * y))
    return par, cov, float((W * (y - X @ par) ** 2).sum()), len(y) - X.shape[1]


def fit_scan(S, y, e, p_fixed=1.0):
    """returns c, its error, p, ln M0, chi2, dof and the multiplicity at the chain coupling"""
    S, y, e = np.asarray(S, float), np.asarray(y, float), np.asarray(e, float)
    if p_fixed is None:
        X = np.stack([np.ones_like(S), np.log(S), -S], 1)
        par, cov, chi2, dof = wls(X, y, e)
        a, p, c = (float(v) for v in par)
        ic = 2
    else:
        X = np.stack([np.ones_like(S), -S], 1)
        par, cov, chi2, dof = wls(X, y - p_fixed * np.log(S), e)
        a, c = (float(v) for v in par)
        p = float(p_fixed)
        ic = 1
    sc = max(1.0, math.sqrt(chi2 / dof)) if dof > 0 else 1.0
    err = sc * np.sqrt(np.diag(cov))
    return {"c": c, "c_err": float(err[ic]), "p": p, "p_err": float(err[1]) if p_fixed is None else 0.0,
            "ln_M0": a, "ln_M0_err": float(err[0]), "chi2": chi2, "dof": dof, "error_scale": sc,
            "ln_M_at_S_U": a + p * math.log(S_U)}


def sector_law(p, pe, tau, nsig=8.0):
    """weight power w and heat time t from ln p_n = const + w ln n - t (n^2 - 1)/4 over the leading sectors"""
    p, pe = np.asarray(p, float), np.asarray(pe, float)
    n = np.arange(1, p.size + 1)
    ok = (p > nsig * pe) & (p > 1e-7)
    last = int(np.argmax(~ok)) if (~ok).any() else p.size
    use = np.arange(last)
    if use.size < 4:
        return None
    X = np.stack([np.ones(use.size), np.log(n[use]), -(n[use] ** 2 - 1) / 4.0], 1)
    par, cov, chi2, dof = wls(X, np.log(p[use]), pe[use] / p[use])
    sc = max(1.0, math.sqrt(chi2 / dof)) if dof > 0 else 1.0
    X2 = np.stack([np.ones(use.size), -(n[use] ** 2 - 1) / 4.0], 1)
    par2, cov2, chi22, dof2 = wls(X2, np.log(p[use]) - 2 * np.log(n[use]), pe[use] / p[use])
    sc2 = max(1.0, math.sqrt(chi22 / dof2)) if dof2 > 0 else 1.0
    return {"w": float(par[1]), "w_err": float(sc * math.sqrt(cov[1, 1])), "tau_over_t": float(tau / par[2]),
            "tau_over_t_w2": float(tau / par2[1]), "tau_over_t_w2_err": float(sc2 * math.sqrt(cov2[1, 1]) * tau / par2[1] ** 2),
            "sectors_used": int(use.size), "mean_casimir": float((p * (n * n - 1) / 4).sum()),
            "c_E_second_moment": float((p * (n * n - 1) / 4).sum() * 2 * tau / 3),
            "alternating_sum": float((p * (-1.0) ** (n - 1)).sum()), "sum": float(p.sum())}


def sector_weights(kind, par, nmax):
    """coefficients f_n of the plaquette weight k = sum_n f_n chi_n in the engine's normalisation, and tau"""
    n = np.arange(1, nmax + 1)
    if kind == 1:
        tau = par / 2
        lnnorm = 0.25 * tau + 0.5 * math.log(4 * math.pi / tau) - math.log(tau)
        return n * np.exp(-tau * (n * n - 1) / 4 - lnnorm), tau
    return 2 * n * special.iv(n, par) / par, 2.0 / par      # exp(beta x) = sum_n (2 n I_n(beta)/beta) chi_n


def estimate(nodes):
    """ln <A> from the block means of a staged run. nodes[k, b] = (mean D, mean exp(+a D), mean exp(-a D),
    mean exp(+2a D), mean exp(-2a D), ...), a = 1/(2N). Returns the half-step overlap estimator with its jackknife
    error, the full-step forward and backward estimators and the Simpson integral of <D>."""
    nodes = np.asarray(nodes, float)
    N = nodes.shape[0] - 1
    nb = nodes.shape[1]

    def jack(fun, b):
        full = fun(b.mean(0))
        tot = b.sum(0)
        v = np.array([fun((tot - b[i]) / (nb - 1)) for i in range(nb)])
        return full, math.sqrt((nb - 1) / nb * ((v - v.mean()) ** 2).sum())
    tot = var = 0.0
    for k in range(N + 1):
        if k == 0:
            v, e = jack(lambda m: math.log(m[1]), nodes[k])
        elif k == N:
            v, e = jack(lambda m: -math.log(m[2]), nodes[k])
        else:
            v, e = jack(lambda m: math.log(m[1]) - math.log(m[2]), nodes[k])
        tot += v; var += e * e
    m = nodes.mean(1)
    wts = np.ones(N + 1)
    if N % 2 == 0:
        wts[1:-1:2] = 4.0; wts[2:-1:2] = 2.0; wts /= 3.0 * N
    else:
        wts[0] = wts[-1] = 0.5; wts /= N
    return {"lnA": tot, "lnA_err": math.sqrt(var), "lnA_forward": float(np.log(m[:N, 3]).sum()),
            "lnA_backward": float(-np.log(m[1:, 4]).sum()), "lnA_TI": float((wts * m[:, 0]).sum())}
