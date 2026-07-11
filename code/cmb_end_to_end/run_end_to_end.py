"""OPH end-to-end CMB pipeline driver (BOLTZ-V1 lane).

Chain: OPH tilt candidate -> primordial power law -> own Einstein-Boltzmann
hierarchy (perturbations.py, Theorem B3) -> own line-of-sight projection
(los.py, Proposition B5) -> D_ell -> comparison against (a) the CAMB curve
stored in the audited 64k run (cross-check only, never in the loop) and
(b) the 83 Planck 2018 binned TT points.

Declared imports (receipt-tagged): Planck 2018 background parameters,
atomic data in recombination (BOLTZ-X1), sigma_T normalization (BOLTZ-N1),
amplitude A_s (gated, #330). Massive-neutrino mass and lensing are NOT
implemented; both are recorded as known-omission terms in the receipt.
"""
from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from background import Params, Background, Recombination
from perturbations import PerturbationSolver
from los import tau_grid, coarse_k_grid, cl_from_sources

HERE = Path(__file__).resolve().parent
RUN = Path("/Users/muellerberndt/Projects/oph-meta/oph-physics-sim/runs/"
           "oph_universe_64k_final_audited_20260711")

P_SRC = 1.630972095694329
PHI = (1 + 5**0.5) / 2
E = float(np.e)

TILT_CANDIDATES = {
    "baseline_planck": 0.9649,                      # validation against CAMB only
    "oph_p_over_48": 1.0 - P_SRC / 48.0,            # repair-depth ansatz
    "oph_e_dP": 1.0 - E * (P_SRC - PHI),            # clock branch at kappa_rep = e
    "oph_empirical_clock": 0.978924,                # 64k lattice, declared unit step
}

ELLS = np.unique(np.round(np.geomspace(2, 2600, 72)).astype(int))

_ps = None
_taus = None


def _init():
    global _ps, _taus
    p = Params()
    bg = Background(p)
    rec = Recombination(bg)
    _ps = PerturbationSolver(bg, rec)
    _taus = tau_grid(bg, rec)
    return p, bg, rec


def _one_k(k):
    try:
        return _ps.sources(k, _taus)
    except Exception as exc:  # record failures, don't kill the pool
        sys.stderr.write(f"k={k}: {exc}\n")
        return np.zeros_like(_taus)


def main():
    p, bg, rec = _init()
    ks = coarse_k_grid()
    cache = HERE / "runtime" / "sources_cache.npz"
    cache.parent.mkdir(exist_ok=True)
    if cache.exists():
        dat = np.load(cache)
        if dat["ks"].size == ks.size and np.allclose(dat["ks"], ks):
            S = dat["S"]
        else:
            S = None
    else:
        S = None
    if S is None:
        with Pool(processes=8, initializer=_init) as pool:
            rows = pool.map(_one_k, ks)
        S = np.vstack(rows)
        np.savez_compressed(cache, ks=ks, S=S, taus=_taus)
        print("sources cached:", cache)

    # reference data from the audited run (cross-check + Planck bins)
    bins = np.genfromtxt(RUN / "finite_repair_clock_cmb_tt_bins.csv",
                         delimiter=",", names=True)
    curves = np.genfromtxt(RUN / "finite_repair_clock_cmb_tt_curves.csv",
                           delimiter=",", names=True)

    report = {
        "mode": "oph_end_to_end_cmb_v0",
        "solver": "own Einstein-Boltzmann hierarchy + LOS (no CAMB/CLASS in loop)",
        "imports": {
            "background": "Planck 2018 baseline H0/ombh2/omch2/tau (declared)",
            "BOLTZ-X1": "Saha+Peebles recombination, atomic data import",
            "BOLTZ-N1": "sigma_T normalization (absolute m_e import)",
            "A_s": "2.1e-9 (gated, #330)",
            "known_omissions": ["mnu=0.06 treated massless", "no lensing",
                                "no RECFAST corrections beyond Peebles"],
        },
        "candidates": {},
    }

    from scipy.interpolate import CubicSpline
    ell_bins = bins["ell"]
    curves_out = {}
    for name, ns in TILT_CANDIDATES.items():
        D = cl_from_sources(bg, ks, _taus, S, ELLS, n_s=ns, A_s=p.A_s)
        curves_out[name] = D
        spl = CubicSpline(ELLS, D)
        D_at_bins = spl(ell_bins)
        obs, sig = bins["observed_D_ell"], bins["minus_dD_ell"]
        chi2 = float(np.sum(((obs - D_at_bins) / sig) ** 2) / len(obs))
        A_fit = float(np.sum(obs * D_at_bins / sig**2)
                      / np.sum(D_at_bins**2 / sig**2))
        chi2_scaled = float(np.sum(((obs - A_fit * D_at_bins) / sig) ** 2) / len(obs))
        entry = {"n_s": ns, "chi2_per_pt_vs_planck_unscaled": chi2,
                 "amplitude_fit": A_fit, "chi2_per_pt_vs_planck_scaled": chi2_scaled}
        if name == "baseline_planck":
            camb = np.interp(ELLS, curves["ell"], curves["camb_lcdm_powerlaw_D_ell"])
            sel = ELLS >= 30
            rel = np.abs(D[sel] - camb[sel]) / camb[sel]
            entry["camb_crosscheck_median_reldev_ell_ge_30"] = float(np.median(rel))
            entry["camb_crosscheck_max_reldev_ell_ge_30"] = float(np.max(rel))
        report["candidates"][name] = entry
        print(f"{name}: ns={ns:.6f} chi2/pt={chi2:.3f} "
              f"(scaled {chi2_scaled:.3f}, A={A_fit:.4f})")

    np.savez_compressed(HERE / "runtime" / "dl_curves.npz",
                        ells=ELLS, **curves_out)
    out = HERE / "runtime" / "end_to_end_receipt.json"
    out.write_text(json.dumps(report, indent=2))
    print("receipt:", out)


if __name__ == "__main__":
    main()
