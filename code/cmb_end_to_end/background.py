"""Background cosmology + recombination for the OPH end-to-end CMB pipeline.

Implements the declared-import layer of cosmology/oph_boltzmann_transport_derivation.tex:
  - flat FRW background (BOLTZ import: H0, omega_b, omega_c from Planck 2018 baseline)
  - Saha + Peebles hydrogen recombination with helium Saha (BOLTZ-X1 declared import:
    atomic data only; no CMB-target information)
  - optical depth, visibility function, tanh reionization matched to tau_reio

Units: Mpc for lengths, conformal time in Mpc. Temperatures in K internally.
No CAMB/CLASS calls anywhere in this package.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq

# physical constants (SI, CODATA-class atomic imports; BOLTZ-X1/N1 class)
C_M_S = 2.99792458e8
G_SI = 6.67430e-11
K_B = 1.380649e-23
HBAR = 1.054571817e-34
M_E = 9.1093837015e-31
M_H = 1.6735575e-27
SIGMA_T_SI = 6.6524587321e-29  # Thomson cross-section, m^2 (alpha^2/m_e^2 import)
EPS0_H = 13.605693122994  # eV, hydrogen ionization
EV = 1.602176634e-19
MPC_M = 3.0856775814913673e22

# derived radiation density for T_cmb
def omega_gamma_h2(t_cmb: float) -> float:
    a_rad = 4.0 * 5.670374419e-8 / C_M_S  # radiation constant J m^-3 K^-4
    rho_g = a_rad * t_cmb**4 / C_M_S**2  # kg/m^3
    rho_crit_h2 = 3.0 * (100.0 * 1000.0 / MPC_M) ** 2 / (8.0 * np.pi * G_SI)
    return rho_g / rho_crit_h2


@dataclass
class Params:
    H0: float = 67.36                 # km/s/Mpc  (import)
    ombh2: float = 0.02237            # (import)
    omch2: float = 0.12               # (import)
    tau_reio: float = 0.0544          # (import)
    t_cmb: float = 2.7255             # K (import)
    n_eff: float = 3.044              # massless nu (mnu=0.06 ignored; ~<1% effect)
    y_he: float = 0.245               # helium mass fraction (import)
    A_s: float = 2.1e-9               # amplitude: GATED input (#330); import here
    n_s: float = 0.9649               # tilt: OPH candidate or baseline (declared per run)
    k_pivot: float = 0.05             # 1/Mpc

    h: float = field(init=False)
    omg: float = field(init=False)
    omn: float = field(init=False)
    omr: float = field(init=False)
    omb: float = field(init=False)
    omc: float = field(init=False)
    oml: float = field(init=False)

    def __post_init__(self):
        self.h = self.H0 / 100.0
        og_h2 = omega_gamma_h2(self.t_cmb)
        on_h2 = self.n_eff * (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * og_h2
        self.omg = og_h2 / self.h**2
        self.omn = on_h2 / self.h**2
        self.omr = self.omg + self.omn
        self.omb = self.ombh2 / self.h**2
        self.omc = self.omch2 / self.h**2
        self.oml = 1.0 - self.omr - self.omb - self.omc

    # H0 in 1/Mpc
    @property
    def H0_mpc(self) -> float:
        return self.H0 * 1000.0 / C_M_S / 1.0  # (km/s/Mpc)/(c) -> 1/Mpc


class Background:
    """Flat FRW background with conformal-time tables."""

    def __init__(self, p: Params, a_init: float = 1e-9, n_grid: int = 40000):
        self.p = p
        la = np.linspace(np.log(a_init), 0.0, n_grid)
        a = np.exp(la)
        H = self.H_of_a(a)                       # 1/Mpc (conformal H = a H_proper handled below)
        # conformal time: dtau = da/(a^2 H_proper); H_proper(a) in 1/Mpc
        integrand = 1.0 / (a**2 * H)
        tau = cumulative_trapezoid(integrand * a, la, initial=0.0)  # d(ln a) integration
        # add the analytic radiation-era offset from 0 to a_init
        tau0_off = a_init / (a_init**2 * self.H_of_a(np.array([a_init]))[0]) * 1.0
        tau = tau + tau0_off
        self.a_grid, self.tau_grid = a, tau
        self.tau_of_a = interp1d(la, tau, kind="cubic")
        self.a_of_tau = interp1d(tau, a, kind="cubic")
        self.tau0 = float(tau[-1])

    def H_of_a(self, a: np.ndarray) -> np.ndarray:
        """proper H(a) in 1/Mpc."""
        p = self.p
        return p.H0_mpc * np.sqrt(p.omr * a**-4 + (p.omb + p.omc) * a**-3 + p.oml)

    def aH(self, a):  # conformal Hubble in 1/Mpc
        return a * self.H_of_a(np.asarray(a))

    def rho_ratios(self, a):
        """Return (rho+p)-weighted densities in units of 3 H0^2 / (8 pi G) = 1."""
        p = self.p
        return dict(
            g=p.omg * a**-4, n=p.omn * a**-4,
            b=p.omb * a**-3, c=p.omc * a**-3, l=p.oml,
        )


class Recombination:
    """Saha + Peebles hydrogen with Saha helium; tanh reionization matched to tau."""

    def __init__(self, bg: Background, n_z: int = 20000):
        self.bg = bg
        p = bg.p
        z = np.linspace(2500.0, 0.0, n_z)
        a = 1.0 / (1.0 + z)
        T = p.t_cmb / a  # K
        # number densities (m^-3)
        rho_crit = 3.0 * (p.H0 * 1000.0 / MPC_M) ** 2 / (8 * np.pi * G_SI)
        n_b = p.omb * rho_crit / M_H * a**-3.0
        n_H = (1.0 - p.y_he) * n_b
        f_He = p.y_he / (4.0 * (1.0 - p.y_he))  # n_He/n_H

        # --- helium + hydrogen Saha to seed, then Peebles ODE for hydrogen ---
        def saha_xe(zi, Ti, nHi):
            # hydrogen Saha: x^2/(1-x) = S
            lam_T = np.sqrt(2 * np.pi * M_E * K_B * Ti) / (2 * np.pi * HBAR)
            S = lam_T**3 / nHi * np.exp(-EPS0_H * EV / (K_B * Ti))
            if S > 1e8:
                return 1.0
            return (-S + np.sqrt(S * S + 4 * S)) / 2.0

        # Peebles ODE in ln a
        LAM_2S = 8.227  # s^-1
        LYA = 121.5673e-9  # m

        def alpha_B(T4):  # m^3/s, Pequignot fit * 1.14 fudge (RECFAST class)
            t = T4 / 1.0e4
            return 1.14 * 4.309e-19 * t**-0.6166 / (1.0 + 0.6703 * t**0.5300)

        def beta_full(T):  # ionization coefficient, detailed balance with full eps0
            # inv_lam = 1/lambda_dB, so (m kT/2 pi hbar^2)^{3/2} = inv_lam**3
            inv_lam = np.sqrt(2 * np.pi * M_E * K_B * T) / (2 * np.pi * HBAR)
            return alpha_B(T) * inv_lam**3 * np.exp(-EPS0_H * EV / (K_B * T))

        def beta2(T):  # photoionization from n=2 (enters the C factor only)
            inv_lam = np.sqrt(2 * np.pi * M_E * K_B * T) / (2 * np.pi * HBAR)
            return alpha_B(T) * inv_lam**3 * np.exp(-EPS0_H * EV / (4.0 * K_B * T))

        def dxe_dlna(lna, xe):
            ai = np.exp(lna)
            Ti = p.t_cmb / ai
            nHi = (1.0 - p.y_he) * p.omb * rho_crit / M_H * ai**-3
            Hi = bg.H_of_a(np.array([ai]))[0] * C_M_S / MPC_M  # s^-1
            x = min(max(xe[0], 1e-10), 1.0 + 2 * f_He)
            xH = min(x, 1.0)
            # Peebles C factor (Ly-alpha escape vs 2s two-photon vs re-ionization from n=2)
            n1s = nHi * max(1.0 - xH, 1e-12)
            K = LYA**3 * n1s / (8.0 * np.pi * Hi)
            KL = K * LAM_2S
            Kb = K * beta2(Ti)
            C = (1.0 + KL) / (1.0 + KL + Kb)
            dx = C / Hi * (beta_full(Ti) * (1.0 - xH)
                           - alpha_B(Ti) * x * xH * nHi)
            return [dx]

        # seed with Saha until x_e < 0.985
        xe = np.zeros_like(z)
        i0 = 0
        for i, (zi, Ti, nHi) in enumerate(zip(z, T, n_H)):
            xs = saha_xe(zi, Ti, nHi)
            # add helium (singly ionized until z~2500, fully coupled; double He ignored below z 2500)
            xe[i] = xs + f_He * (1.0 if zi > 2200 else self._he_saha(Ti, nHi, f_He))
            if xs < 0.985:
                i0 = i
                break
        # Peebles from i0
        sol = solve_ivp(dxe_dlna, [np.log(a[i0]), 0.0], [xe[i0] - f_He * self._he_saha(T[i0], n_H[i0], f_He)],
                        t_eval=np.log(a[i0:]), method="LSODA", rtol=1e-8, atol=1e-12)
        xH_evo = np.clip(sol.y[0], 0.0, 1.0)
        xe[i0:] = xH_evo + f_He * np.array([self._he_saha(Ti, nHi, f_He) for Ti, nHi in zip(T[i0:], n_H[i0:])])

        # reionization: tanh in y=(1+z)^{3/2}, width dz=0.5, amplitude (1+f_He)
        self._z, self._xe_norec = z, xe.copy()
        self._nH_of_a = lambda ai: (1.0 - p.y_he) * p.omb * rho_crit / M_H * ai**-3.0
        self._f_He = f_He
        zre = self._solve_zre()
        xe_re = self._reio_profile(z, zre)
        self.xe = np.maximum(xe, xe_re)
        self.xe_of_z = interp1d(z, self.xe, kind="cubic", fill_value=(self.xe[-1], self.xe[0]), bounds_error=False)
        self.zre = zre

        # opacity kdot = a n_e sigma_T (conformal, 1/Mpc) on the background tau grid
        a_g = bg.a_grid
        z_g = 1.0 / a_g - 1.0
        xe_g = np.where(z_g > 2499.0, 1.0 + 2 * f_He, self.xe_of_z(np.clip(z_g, 0, 2499.9)))
        n_e = xe_g * self._nH_of_a(a_g)
        self.kdot_grid = n_e * SIGMA_T_SI * a_g * MPC_M  # 1/Mpc conformal
        # optical depth kappa(tau, tau0)
        kap = -cumulative_trapezoid(self.kdot_grid[::-1], bg.tau_grid[::-1], initial=0.0)[::-1]
        self.kappa_grid = kap
        self.expmk = np.exp(-kap)
        self.g_grid = self.kdot_grid * self.expmk
        self.kdot_of_tau = interp1d(bg.tau_grid, self.kdot_grid, kind="cubic")
        self.g_of_tau = interp1d(bg.tau_grid, self.g_grid, kind="cubic")
        self.expmk_of_tau = interp1d(bg.tau_grid, self.expmk, kind="cubic")

    @staticmethod
    def _he_saha(T, nH, f_He):
        # singly-ionized helium fraction via Saha (24.587 eV); returns x_He in [0,1]
        chi = 24.587387 * EV
        lam_T = np.sqrt(2 * np.pi * M_E * K_B * T) / (2 * np.pi * HBAR)
        S = 4.0 * lam_T**3 / nH * np.exp(-chi / (K_B * T))
        if S > 1e8:
            return 1.0
        return S / (1.0 + S)

    def _reio_profile(self, z, zre, dz=0.5):
        y = (1.0 + z) ** 1.5
        yre = (1.0 + zre) ** 1.5
        dy = 1.5 * np.sqrt(1.0 + zre) * dz
        return (1.0 + self._f_He) * 0.5 * (1.0 + np.tanh((yre - y) / dy))

    def _solve_zre(self):
        bg, p = self.bg, self.bg.p
        rho_crit = 3.0 * (p.H0 * 1000.0 / MPC_M) ** 2 / (8 * np.pi * G_SI)

        def tau_of_zre(zre):
            z = np.linspace(0.0, 50.0, 4000)
            a = 1.0 / (1 + z)
            xe = np.maximum(self._reio_profile(z, zre),
                            interp1d(self._z, self._xe_norec, bounds_error=False,
                                     fill_value=(self._xe_norec[-1], 1.0))(z))
            n_e = xe * self._nH_of_a(a)
            # dtau = n_e sigma_T c dt = n_e sigma_T dz c/((1+z) H)
            Hz = bg.H_of_a(a) * C_M_S / MPC_M
            integ = n_e * SIGMA_T_SI * C_M_S / ((1 + z) * Hz)
            return np.trapezoid(integ, z)

        return brentq(lambda zz: tau_of_zre(zz) - p.tau_reio, 3.0, 15.0, xtol=1e-3)
