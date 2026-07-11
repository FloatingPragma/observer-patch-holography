"""Per-k linear Einstein-Boltzmann integrator (conformal Newtonian gauge).

Implements Theorem B3 of cosmology/oph_boltzmann_transport_derivation.tex:
the Ma-Bertschinger (1995) system for photons (temperature + polarization),
massless neutrinos, baryons, CDM, and the two Einstein constraints, with
adiabatic initial conditions normalized to comoving curvature R = 1.

The OPH insertion points Gamma_rec and B_A are present as explicit zero
hooks (gated by #374); with both zero the system is standard.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

LMAX_G = 12   # photon temperature multipoles 0..LMAX_G
LMAX_P = 12   # photon polarization
LMAX_N = 12   # massless neutrinos


def state_layout():
    n = 0
    idx = {}
    for name, size in [("phi", 1), ("dc", 1), ("tc", 1), ("db", 1), ("tb", 1),
                       ("F", LMAX_G + 1), ("G", LMAX_P + 1), ("N", LMAX_N + 1)]:
        idx[name] = slice(n, n + size)
        n += size
    return idx, n


IDX, NDIM = state_layout()


class PerturbationSolver:
    def __init__(self, bg, rec):
        self.bg, self.rec = bg, rec
        p = bg.p
        self.H02 = p.H0_mpc ** 2
        self.om = dict(g=p.omg, n=p.omn, b=p.omb, c=p.omc, l=p.oml)
        self.Rnu = p.omn / (p.omg + p.omn)

    # ---- background helpers ----
    def _aH(self, a):
        return a * self.bg.H_of_a(np.asarray(a).reshape(1))[0]

    def _rho(self, a):
        """4 pi G a^2 rho_i in Mpc^-2: (3/2) H0^2 om_i a^2 * a^-3(or 4)."""
        f = 1.5 * self.H02
        return dict(
            g=f * self.om["g"] * a**-2, n=f * self.om["n"] * a**-2,
            b=f * self.om["b"] * a**-1, c=f * self.om["c"] * a**-1,
            l=f * self.om["l"] * a**2,
        )

    # ---- metric constraints ----
    def _psi_phidot(self, a, k, y, kdot):
        r = self._rho(a)
        aH = self._aH(a)
        F, G, N = y[IDX["F"]], y[IDX["G"]], y[IDX["N"]]
        sig_g, sig_n = F[2] / 2.0, N[2] / 2.0
        # anisotropic stress: k^2(phi - psi) = 12 pi G a^2 (rho+p) sigma
        shear = 3.0 * ((4.0 / 3.0) * r["g"] * sig_g + (4.0 / 3.0) * r["n"] * sig_n)
        phi = y[IDX["phi"]][0]
        psi = phi - shear / k**2
        # momentum constraint: phidot = -aH psi + 4 pi G a^2 (rho+p) theta / k^2
        th_g = 0.75 * k * F[1]
        th_n = 0.75 * k * N[1]
        mom = (r["c"] * y[IDX["tc"]][0] + r["b"] * y[IDX["tb"]][0]
               + (4.0 / 3.0) * r["g"] * th_g + (4.0 / 3.0) * r["n"] * th_n)
        phidot = -aH * psi + mom / k**2
        return psi, phidot

    # ---- RHS in ln a ----
    def rhs(self, lna, y, k):
        a = np.exp(lna)
        aH = self._aH(a)
        kdot = float(self.rec.kdot_of_tau(np.clip(self.bg.tau_of_a(lna), None, self.bg.tau0)))
        psi, phidot = self._psi_phidot(a, k, y, kdot)
        F, G, N = y[IDX["F"]], y[IDX["G"]], y[IDX["N"]]
        th_g = 0.75 * k * F[1]
        tb = y[IDX["tb"]][0]
        tc = y[IDX["tc"]][0]
        r = self._rho(a)
        Rgb = (4.0 * r["g"]) / (3.0 * r["b"])  # 4 rho_g / 3 rho_b
        Pi = F[2] + G[0] + G[2]

        dy = np.zeros_like(y)
        dy[IDX["phi"]] = phidot
        dy[IDX["dc"]] = -tc + 3.0 * phidot
        dy[IDX["tc"]] = -aH * tc + k**2 * psi
        dy[IDX["db"]] = -tb + 3.0 * phidot
        dy[IDX["tb"]] = (-aH * tb + k**2 * psi
                         + Rgb * kdot * (th_g - tb))

        dF = np.zeros(LMAX_G + 1)
        dF[0] = -k * F[1] + 4.0 * phidot
        dF[1] = (k / 3.0) * F[0] - (2.0 * k / 3.0) * F[2] + (4.0 * k / 3.0) * psi \
            + (4.0 / (3.0 * k)) * kdot * tb - kdot * F[1]
        dF[2] = (2.0 * k / 5.0) * F[1] - (3.0 * k / 5.0) * F[3] \
            - 0.9 * kdot * F[2] + 0.1 * kdot * (G[0] + G[2])
        for l in range(3, LMAX_G):
            dF[l] = k / (2 * l + 1.0) * (l * F[l - 1] - (l + 1) * F[l + 1]) - kdot * F[l]
        l = LMAX_G
        tau = float(self.bg.tau_of_a(lna))
        dF[l] = k * F[l - 1] - (l + 1) / tau * F[l] - kdot * F[l]
        dy[IDX["F"]] = dF

        dG = np.zeros(LMAX_P + 1)
        dG[0] = -k * G[1] + kdot * (-G[0] + 0.5 * Pi)
        for l in range(1, LMAX_P):
            src = 0.1 * kdot * Pi if l == 2 else 0.0
            dG[l] = k / (2 * l + 1.0) * (l * G[l - 1] - (l + 1) * G[l + 1]) - kdot * G[l] + src
        l = LMAX_P
        dG[l] = k * G[l - 1] - (l + 1) / tau * G[l] - kdot * G[l]
        dy[IDX["G"]] = dG

        dN = np.zeros(LMAX_N + 1)
        dN[0] = -k * N[1] + 4.0 * phidot
        dN[1] = (k / 3.0) * N[0] - (2.0 * k / 3.0) * N[2] + (4.0 * k / 3.0) * psi
        for l in range(2, LMAX_N):
            dN[l] = k / (2 * l + 1.0) * (l * N[l - 1] - (l + 1) * N[l + 1])
        l = LMAX_N
        dN[l] = k * N[l - 1] - (l + 1) / tau * N[l]
        dy[IDX["N"]] = dN

        # convert d/dtau -> d/dlna
        return dy / aH

    # ---- adiabatic initial conditions, R = 1 ----
    def initial(self, k):
        # start where k*tau small and kdot*tau large
        tau_i = min(0.05 / k, 1.0)
        lna_i = None
        # invert tau->a via background table
        a_i = float(self.bg.a_of_tau(tau_i))
        lna_i = np.log(a_i)
        Rnu = self.Rnu
        psi = 1.0 / (1.5 * (1.0 + 4.0 * Rnu / 15.0))  # normalizes R = 1
        phi = (1.0 + 2.0 * Rnu / 5.0) * psi
        y = np.zeros(NDIM)
        y[IDX["phi"]] = phi
        dg = -2.0 * psi
        y[IDX["dc"]] = 0.75 * dg
        y[IDX["db"]] = 0.75 * dg
        th = (k**2 * tau_i / 2.0) * psi
        y[IDX["tc"]] = th
        y[IDX["tb"]] = th
        F = np.zeros(LMAX_G + 1)
        F[0] = dg
        F[1] = (4.0 / (3.0 * k)) * th
        y[IDX["F"]] = F
        N = np.zeros(LMAX_N + 1)
        N[0] = dg
        N[1] = (4.0 / (3.0 * k)) * th
        N[2] = 2.0 * (k * tau_i) ** 2 / 30.0 * psi  # sigma_nu = (k tau)^2 psi/15
        y[IDX["N"]] = N
        return lna_i, y

    # ---- integrate one k, return sources on tau grid ----
    def sources(self, k, tau_out):
        lna_i, y0 = self.initial(k)
        tau_i = float(self.bg.tau_of_a(lna_i))
        mask = tau_out > tau_i * 1.05
        lna_out = np.log(np.asarray(self.bg.a_of_tau(tau_out[mask]), dtype=float))
        sol = solve_ivp(self.rhs, [lna_i, 0.0], y0, args=(k,), method="BDF",
                        t_eval=lna_out, rtol=3e-6, atol=1e-9,
                        dense_output=False)
        if not sol.success:
            raise RuntimeError(f"k={k}: {sol.message}")
        Y = sol.y
        a = np.exp(sol.t)
        tau = np.asarray(self.bg.tau_of_a(sol.t), dtype=float)
        g = self.rec.g_of_tau(tau)
        expmk = self.rec.expmk_of_tau(tau)
        kdot = self.rec.kdot_of_tau(tau)

        phi = Y[IDX["phi"], :][0]
        F0 = Y[IDX["F"], :][0]
        F1 = Y[IDX["F"], :][1]
        F2 = Y[IDX["F"], :][2]
        G0 = Y[IDX["G"], :][0]
        G2 = Y[IDX["G"], :][2]
        tb = Y[IDX["tb"], :][0]

        # psi, phidot pointwise
        psi = np.zeros_like(phi)
        phidot = np.zeros_like(phi)
        for j in range(Y.shape[1]):
            ps, pd = self._psi_phidot(a[j], k, Y[:, j], kdot[j])
            psi[j], phidot[j] = ps, pd
        psidot = np.gradient(psi, tau)

        Theta0 = F0 / 4.0
        Pi = (F2 + G0 + G2) / 4.0
        ub = tb / k

        S_mono = g * (Theta0 + psi + Pi / 4.0)
        S_isw = expmk * (phidot + psidot)
        S_dopp = np.gradient(g * ub, tau) / k
        S_pol = 0.75 / k**2 * np.gradient(np.gradient(g * Pi, tau), tau)
        S_full = np.zeros_like(tau_out)
        S_full[mask] = S_mono + S_isw + S_dopp + S_pol
        return S_full
