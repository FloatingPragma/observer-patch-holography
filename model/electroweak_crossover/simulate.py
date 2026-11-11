"""Run one point of the lattice theory and return the measurement series.

    python -m model.electroweak_crossover.simulate --n 24 --beta-G 8 --x 0.296 --y -0.09 --sweeps 20000 --out run.npz

The compiled kernels are used when numba is installed. Otherwise the numpy implementation runs.
"""
from __future__ import annotations

import argparse
import math

import numpy as np

from .lattice import Lattice, lattice_couplings


def run_point(n_side: int, beta_G: float, x: float, y: float, n_therm: int, n_meas: int, seed: int,
              n_overrelax: int = 4, start: str = "cold", radius: float = 1.0, backend: str = "auto") -> dict:
    """Series of (H, S2, S4) per compound sweep at the continuum point (x, y) and lattice spacing 4/(beta_G g3^2)."""
    beta_H, beta_R = lattice_couplings(x, y, beta_G)
    if backend == "auto":
        try:
            from . import kernels
            backend = "numba" if kernels.AVAILABLE else "numpy"
        except Exception:
            backend = "numpy"
    if backend == "numba":
        from . import kernels
        rng = np.random.default_rng(seed)
        kernels.seed_numba(seed)
        U, phi = kernels.cold_start(n_side, radius) if start == "cold" else kernels.hot_start(n_side, rng, radius)
        width = 0.5 / math.sqrt(n_side ** 3)
        kernels.run_series(U, phi, n_side, beta_G, beta_H, beta_R, n_overrelax, width, n_therm, 0, True, True)
        series, _ = kernels.run_series(U, phi, n_side, beta_G, beta_H, beta_R, n_overrelax, width, n_meas, 0, True, True)
        series = series[:, :3]
    else:
        lat = Lattice(n_side, beta_G, beta_H, beta_R, np.random.default_rng(seed), start=start, radius=radius)
        for _ in range(n_therm):
            lat.sweep(n_overrelax)
        series = lat.run(n_meas, n_overrelax)
    return dict(n_side=n_side, beta_G=beta_G, x=x, y=y, beta_H=beta_H, beta_R=beta_R, series=series, backend=backend)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--beta-G", type=float, required=True)
    p.add_argument("--x", type=float, required=True)
    p.add_argument("--y", type=float, required=True)
    p.add_argument("--therm", type=int, default=500)
    p.add_argument("--sweeps", type=int, default=1000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--backend", default="auto")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    r = run_point(a.n, a.beta_G, a.x, a.y, a.therm, a.sweeps, a.seed, backend=a.backend)
    np.savez_compressed(a.out, **{k: v for k, v in r.items()})
    r2 = r["series"][:, 1] / a.n ** 3
    print(f"backend {r['backend']}  beta_H {r['beta_H']:.7f}  beta_R {r['beta_R']:.7f}  <R^2> {r2.mean():.5f}")


if __name__ == "__main__":
    main()
