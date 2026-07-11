"""Diagnostic: per-source-term D_ell decomposition (not part of the receipt)."""
import numpy as np
from multiprocessing import Pool
from scipy.integrate import solve_ivp

from background import Params, Background, Recombination
from perturbations import PerturbationSolver, IDX
from los import tau_grid, coarse_k_grid, cl_from_sources

_ps = _taus = _bg = _rec = None


def _init():
    global _ps, _taus, _bg, _rec
    p = Params()
    _bg = Background(p)
    _rec = Recombination(_bg)
    _ps = PerturbationSolver(_bg, _rec)
    _taus = tau_grid(_bg, _rec)


def _terms_one_k(k):
    ps, taus, bg, rec = _ps, _taus, _bg, _rec
    lna_i, y0 = ps.initial(k)
    tau_i = float(bg.tau_of_a(lna_i))
    mask = taus > tau_i * 1.05
    lna_out = np.log(np.asarray(bg.a_of_tau(taus[mask]), dtype=float))
    sol = solve_ivp(ps.rhs, [lna_i, 0.0], y0, args=(k,), method="BDF",
                    t_eval=lna_out, rtol=3e-6, atol=1e-9)
    Y = sol.y
    a = np.exp(sol.t)
    tau = np.asarray(bg.tau_of_a(sol.t), dtype=float)
    g = rec.g_of_tau(tau)
    expmk = rec.expmk_of_tau(tau)
    F0, F2 = Y[IDX["F"], :][0], Y[IDX["F"], :][2]
    G0, G2 = Y[IDX["G"], :][0], Y[IDX["G"], :][2]
    tb = Y[IDX["tb"], :][0]
    psi = np.zeros_like(F0)
    phidot = np.zeros_like(F0)
    for j in range(Y.shape[1]):
        psi[j], phidot[j] = ps._psi_phidot(a[j], k, Y[:, j], 0.0)
    psidot = np.gradient(psi, tau)
    Theta0 = F0 / 4.0
    Pi = (F2 + G0 + G2) / 4.0
    ub = tb / k
    terms = dict(
        mono=g * (Theta0 + psi + Pi / 4.0),
        isw=expmk * (phidot + psidot),
        dopp=np.gradient(g * ub, tau) / k,
        pol=0.75 / k**2 * np.gradient(np.gradient(g * Pi, tau), tau),
    )
    out = {}
    for name, s in terms.items():
        full = np.zeros_like(taus)
        full[mask] = s
        out[name] = full
    return out


def main():
    _init()
    ks = coarse_k_grid()
    with Pool(8, initializer=_init) as pool:
        rows = pool.map(_terms_one_k, ks)
    names = ["mono", "isw", "dopp", "pol"]
    S = {n: np.vstack([r[n] for r in rows]) for n in names}
    ells = np.array([31, 62, 126, 230, 345, 470, 700, 1400])
    p = Params()
    combos = {
        "full": S["mono"] + S["isw"] + S["dopp"] + S["pol"],
        "no_isw": S["mono"] + S["dopp"] + S["pol"],
        "mono_only": S["mono"],
        "no_dopp": S["mono"] + S["isw"] + S["pol"],
        "no_pol": S["mono"] + S["isw"] + S["dopp"],
    }
    ref = np.genfromtxt("/Users/muellerberndt/Projects/oph-meta/oph-physics-sim/runs/"
                        "oph_universe_64k_final_audited_20260711/"
                        "finite_repair_clock_cmb_tt_curves.csv",
                        delimiter=",", names=True)
    camb = np.interp(ells, ref["ell"], ref["camb_lcdm_powerlaw_D_ell"])
    print("ell:      " + "  ".join(f"{l:8d}" for l in ells))
    print("camb:     " + "  ".join(f"{v:8.1f}" for v in camb))
    for name, Sc in combos.items():
        D = cl_from_sources(_bg, ks, _taus, Sc, ells, n_s=0.9649, A_s=p.A_s)
        print(f"{name:9s} " + "  ".join(f"{v:8.1f}" for v in D))


if __name__ == "__main__":
    main()
