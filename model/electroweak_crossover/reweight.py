"""Multi-histogram reweighting (Ferrenberg, Swendsen, Phys. Rev. Lett. 63 (1989) 1195) with block jackknife.

A run stores per measurement H = sum of hopping terms, S2 = sum R^2, S4 = sum R^4. At fixed beta_G the action
depends on the scanned parameters through  -beta_H H + beta_R (S4 - 2 S2),  so every series can be reweighted
to any nearby (beta_H, beta_R), that is to any nearby (x, y).
"""
import math

import numpy as np

from .lattice import condensate, lattice_couplings


def logsumexp(a, axis=None):
    m = np.max(a, axis=axis, keepdims=True)
    return np.squeeze(np.log(np.sum(np.exp(a - m), axis=axis, keepdims=True)) + m, axis=axis)


def tau_int(a):
    a = np.asarray(a, float); a = a - a.mean(); n = len(a); v = a.var()
    if v == 0 or n < 16:
        return 0.5
    t = 0.5
    for k in range(1, n // 4):
        t += np.dot(a[:-k], a[k:]) / ((n - k) * v)
        if k > 6 * t:
            break
    return max(t, 0.5)


class RunSet:
    """Runs at one lattice size and one beta_G. Each run: dict(beta_H, beta_R, series) with series[:, 0:3] = H, S2, S4."""

    def __init__(self, runs, n_side, beta_G, nblock=20):
        self.runs = []
        self.L = n_side; self.bG = beta_G
        for r in runs:
            s = np.asarray(r['series'], float)
            if len(s) < 4 * nblock:
                continue
            self.runs.append(dict(bH=r['beta_H'], bR=r['beta_R'], H=s[:, 0].copy(), S2=s[:, 1].copy(),
                                  Q=(s[:, 2] - 2 * s[:, 1]), n=len(s), tau=tau_int(s[:, 1])))
        self.V = self.L ** 3
        self.nblock = nblock
        self.bH0 = np.mean([r['bH'] for r in self.runs]); self.bR0 = np.mean([r['bR'] for r in self.runs])
        self.H0 = np.mean(np.concatenate([r['H'] for r in self.runs]))
        self.Q0 = np.mean(np.concatenate([r['Q'] for r in self.runs]))
        self._prep(None)

    def _prep(self, drop):
        """Concatenate data, optionally dropping block `drop` of every run."""
        Hs, Qs, S2s, ns = [], [], [], []
        for r in self.runs:
            n = r['n']; nb = self.nblock; bl = n // nb
            if drop is None:
                sl = np.ones(n, bool)
            else:
                sl = np.ones(n, bool); sl[drop * bl:(drop + 1) * bl] = False
            Hs.append(r['H'][sl] - self.H0); Qs.append(r['Q'][sl] - self.Q0); S2s.append(r['S2'][sl])
            ns.append(float(sl.sum()))          # plain counts: the pooled sample is the mixture of the runs
        self.H = np.concatenate(Hs); self.Q = np.concatenate(Qs); self.S2 = np.concatenate(S2s)
        self.offsets = np.concatenate([[0], np.cumsum([len(h) for h in Hs])])
        self.logN = np.log(np.array(ns))
        self.pb = np.array([[r['bH'] - self.bH0, r['bR'] - self.bR0] for r in self.runs])
        # minus reduced action of every config at every run parameter: K x N
        self.ms = self.pb[:, :1] * self.H[None, :] - self.pb[:, 1:] * self.Q[None, :]

    def solve(self, f0=None, tol=1e-10, itmax=200):
        """Free energies of the runs. Newton iteration on sum_n W_k(n) = 1 with
        W_k(n) = exp(f_k - s_k(n)) / sum_j N_j exp(f_j - s_j(n)); f_0 = 0 fixes the zero mode."""
        K = len(self.runs)
        ntot = self.ms.shape[1]
        logN = self.logN + math.log(ntot) - logsumexp(self.logN)      # weights rescaled to sum to ntot
        Nw = np.exp(logN)
        if f0 is None:
            # start: f_(k+1) - f_k = -1/2 (<D>_k + <D>_(k+1)), D = s_k - s_(k+1), along the runs ordered in beta_H
            order = np.argsort(self.pb[:, 0])
            f = np.zeros(K)
            for a, b in zip(order[:-1], order[1:]):
                D = self.ms[b] - self.ms[a]
                ia = slice(self.offsets[a], self.offsets[a + 1]); ib = slice(self.offsets[b], self.offsets[b + 1])
                f[b] = f[a] - 0.5 * (D[ia].mean() + D[ib].mean())
        else:
            f = f0.copy()
        f = f - f[0]
        if K > 1:
            for it in range(itmax):
                logden = logsumexp(logN[:, None] + f[:, None] + self.ms, axis=0)
                Wk = np.exp(f[:, None] + self.ms - logden[None, :])
                S = Wk.sum(1)
                g = S - 1.0
                J = np.diag(S) - (Wk @ Wk.T) * Nw[None, :]
                try:
                    d = np.linalg.solve(J[1:, 1:], -g[1:])
                except np.linalg.LinAlgError:
                    d = -g[1:] / S[1:]
                mx = np.max(np.abs(d))
                if mx > 2.0:
                    d *= 2.0 / mx
                f[1:] += d
                if mx < tol:
                    break
            self.resid = float(np.max(np.abs(g)))
        else:
            self.resid = 0.0
        self.f = f
        self.logNw = logN
        self.logden = logsumexp(logN[:, None] + f[:, None] + self.ms, axis=0)
        return f

    def at(self, bH, bR):
        """Mean and variance of the volume average of R^2 at (bH, bR), and the effective sample fraction."""
        lw = (bH - self.bH0) * self.H - (bR - self.bR0) * self.Q - self.logden
        lw -= lw.max()
        w = np.exp(lw); w /= w.sum()
        r2 = self.S2 / self.V
        m = np.dot(w, r2)
        var = np.dot(w, (r2 - m) ** 2)
        m3 = np.dot(w, (r2 - m) ** 3)
        m4 = np.dot(w, (r2 - m) ** 4)
        return m, var, 1.0 / (len(w) * np.dot(w, w)), m3, m4

    def hist(self, bH, bR, bins):
        lw = (bH - self.bH0) * self.H - (bR - self.bR0) * self.Q - self.logden
        lw -= lw.max()
        w = np.exp(lw); w /= w.sum()
        h, e = np.histogram(self.S2 / self.V, bins=bins, weights=w)
        return h, e

    def curve(self, params):
        """params: list of (bH, bR). Returns arrays mean, var, neff."""
        out = np.array([self.at(bH, bR) for bH, bR in params])
        return out[:, 0], out[:, 1], out[:, 2], out[:, 3], out[:, 4]

    def jackknife(self, fn):
        """fn(self) -> array. Returns (value, error) by delete-one-block jackknife over all runs at once."""
        self._prep(None); f_full = self.solve()
        full = np.asarray(fn(self), float)
        vals = []
        for b in range(self.nblock):
            self._prep(b); self.solve(f0=f_full)
            vals.append(np.asarray(fn(self), float))
        self._prep(None); self.solve(f0=f_full)
        vals = np.array(vals)
        err = np.sqrt((self.nblock - 1) * np.mean((vals - vals.mean(0)) ** 2, axis=0))
        return full, err


def peak(xs, ys):
    """Location and height of the maximum of a sampled curve, parabolic refinement."""
    i = int(np.argmax(ys))
    if i == 0 or i == len(ys) - 1:
        return xs[i], ys[i], False
    a, b, c = ys[i - 1], ys[i], ys[i + 1]
    den = a - 2 * b + c
    d = 0.5 * (a - c) / den if den != 0 else 0.0
    return xs[i] + d * (xs[1] - xs[0]), b - 0.25 * (a - c) * d, True


def half_width_sym(xs, ys, ipk=None):
    """Distance from the maximum to the point on the high side where the curve is half the maximum."""
    i = int(np.argmax(ys)) if ipk is None else ipk
    half = 0.5 * ys[i]
    for j in range(i, len(ys) - 1):
        if ys[j] >= half > ys[j + 1]:
            t = (ys[j] - half) / (ys[j] - ys[j + 1])
            return xs[j] + t * (xs[j + 1] - xs[j]) - xs[i]
    return float('nan')


def y_scan(rs, x, ygrid):
    """chi3 = L^3 beta_H^2 / beta_G Var(R2bar) and the continuum condensate along y at fixed x."""
    pars = [lattice_couplings(x, float(y), rs.bG) for y in ygrid]
    m, var, neff, m3, m4 = rs.curve(pars)
    bH = np.array([p[0] for p in pars])
    chi3 = rs.V * bH ** 2 / rs.bG * var
    cond = np.array([condensate(mm, rs.bG, b) for mm, b in zip(m, bH)])
    return m, chi3, cond, neff, m4 / np.maximum(var, 1e-300) ** 2
