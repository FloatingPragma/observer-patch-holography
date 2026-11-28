"""Independent reconstruction and repeated-anchor check. No imports from replay.py."""

from math import exp, sqrt
from pathlib import Path
import json

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import hyp2f1
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
DATA = ROOT / "code/cosmology/paid_trail_bao"


def evaluate(z, kind, a, b):
    # Analytic antiderivative for dust plus constant vacuum density.
    s = 1 + z
    ratio = np.asarray(a)[..., None] / b

    def primitive(x):
        return x * hyp2f1(1 / 3, 1 / 2, 4 / 3, -ratio * x**3) / sqrt(b)

    dm = primitive(s) - primitive(1.0)
    dh = 1 / np.sqrt(b * (1 + ratio * s**3))
    return np.select(
        (kind == "DM_over_rs", kind == "DH_over_rs"), (dm, dh), np.cbrt(z * dm**2 * dh)
    )


def main():
    raw = np.loadtxt(
        DATA / "desi_gaussian_bao_ALL_GCcomb_mean.txt", dtype=str, skiprows=1
    )
    z, y = raw[:, :2].T.astype(float)
    kind = raw[:, 2]
    cov = np.loadtxt(DATA / "desi_gaussian_bao_ALL_GCcomb_cov.txt")
    inv = np.linalg.inv(cov)
    report = json.loads((HERE / "results.json").read_text())
    rc = report["models"]["rc_load"]
    # Invert the independent archived solution rather than using the pricing inputs.
    cert = json.loads((HERE / "rc_load_source_only_certificate.json").read_text())
    x = float(cert["fixed_point"]["X"]["lo"])
    b = exp(-x) * (147.09 * 3.0856775814913673e22 / 1.616255e-35) ** 2

    def objective(a, mask):
        r = y[mask] - evaluate(z[mask], kind[mask], a, b)
        return r @ np.linalg.solve(cov[np.ix_(mask, mask)], r)

    mask = np.ones(len(z), dtype=bool)
    fit = minimize_scalar(
        lambda a: objective(a, mask),
        bracket=(0.0003, 0.00035, 0.0004),
        method="brent",
        options={"xtol": 1e-11},
    )
    assert abs(fit.fun - rc["chi2_full"]) < 2e-6
    assert abs(fit.x - rc["A_full"]) < 1e-10
    (i,) = np.flatnonzero((z == 2.33) & (kind == "DH_over_rs"))
    anchor = (y[i] ** -2 - b) / (1 + z[i]) ** 3
    fixed = objective(anchor, mask)
    assert abs(fixed - rc["chi2_at_anchor"]) < 2e-6

    # A full 13-dimensional Gaussian draw with each mock anchor recalibrated.
    # Its full quadratic form is the conditional quadratic form because the
    # simulated predicted high-z radial reading is exactly that draw's anchor.
    rng = np.random.default_rng(280964)
    chol = np.linalg.cholesky(cov)
    mu = evaluate(z, kind, fit.x, b)
    extreme = 0
    count = 65536
    for _ in range(count // 1024):
        simulated = mu + rng.standard_normal((1024, len(z))) @ chol.T
        a = (simulated[:, i] ** -2 - b) / (1 + z[i]) ** 3
        residual = simulated - evaluate(z, kind, a, b)
        score = np.einsum("bi,ij,bj->b", residual, inv, residual)
        extreme += np.count_nonzero(score >= fixed)
    integral = rc["repeated_anchor_tail"][1]["p_repeating_anchor"]
    stderr = sqrt(integral * (1 - integral) / count)
    assert abs(extreme / count - integral) < 5 * stderr

    # Dust-only shape obeys exact linear Gaussian contrasts at paired coordinates.
    rows = []
    for j in np.flatnonzero(kind == "DM_over_rs"):
        (radial,) = np.flatnonzero((z == z[j]) & (kind == "DH_over_rs"))
        w = (
            np.eye(len(z))[j]
            - 2 * (1 + z[j]) * (sqrt(1 + z[j]) - 1) * np.eye(len(z))[radial]
        )
        rows.append(w)
    w = np.array(rows)
    d = w @ y
    q = d @ np.linalg.solve(w @ cov @ w.T, d)
    assert abs(q - 406.988328) < 1e-5

    # Entire-redshift omissions are separate from a single-coordinate anchor.
    for value in np.unique(z):
        use = z != value
        f = minimize_scalar(
            lambda a: objective(a, use),
            bracket=(0.00025, 0.00035, 0.00045),
            method="brent",
            options={"xtol": 1e-11},
        )
        check = rc["checks"][f"omit_{value}"]
        assert abs(f.x - check["A_training"]) < 2e-10
        # Full quadratic at the training estimate factors into train + conditional held.
        assert abs(objective(f.x, mask) - f.fun - check["fixed_A_chi2"]) < 1e-5
    print(
        f"Hypergeometric replay: {fixed:.6f}, {fit.fun:.6f}; "
        f"remeasured-anchor MC {extreme}/{count} = {extreme / count:.5f} (SE {stderr:.5f}); "
        f"flat dust-only shape {q:.5f}/6, p={chi2.sf(q, 6):.3g}; seven splits verified."
    )


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without -O to enable the scientific checks.")
    main()
