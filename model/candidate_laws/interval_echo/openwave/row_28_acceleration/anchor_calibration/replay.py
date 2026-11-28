"""Reproduce the fixed-amplitude OPH / DESI comparison from bundled primary records."""

import argparse
import hashlib
import json
from math import cos, exp, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2, ncx2, norm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
DATA = ROOT / "code/cosmology/paid_trail_bao"
BASE = "desi_gaussian_bao_ALL_GCcomb_"
MEAN_HASH = "9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585"
COV_HASH = "252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509"
MPC_M = 149597870700 * 648000 / pi * 1e6
C_KMS = 299792.458
R_DRAG = 147.09
PLANCK_M = 1.616255e-35


class Distances:
    def __init__(self, order=48):
        for part, expected in (("mean", MEAN_HASH), ("cov", COV_HASH)):
            path = DATA / (BASE + part + ".txt")
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError(f"Unexpected DESI {part} file")
        table = np.loadtxt(DATA / (BASE + "mean.txt"), dtype=str, skiprows=1)
        self.z, self.y = table[:, :2].astype(float).T
        self.kind = table[:, 2]
        self.cov = np.loadtxt(DATA / (BASE + "cov.txt"))
        self.s = 1 + self.z
        nodes, weights = np.polynomial.legendre.leggauss(order)
        self.quad_s = 1 + self.z[:, None] * (1 + nodes) / 2
        self.weights = self.z[:, None] * weights / 2

    def theory(self, a, b, c):
        """Return dimensionless DESI distance predictions for any scalar or vector A."""
        a = np.asarray(a)[..., None, None]
        u = self.quad_s
        dm = np.sum(self.weights / np.sqrt(a * u**3 + b + c * u**4), axis=-1)
        dh = 1 / np.sqrt(a[..., 0] * self.s**3 + b + c * self.s**4)
        dv = np.cbrt(self.z * dm**2 * dh)
        return np.where(
            self.kind == "DM_over_rs", dm, np.where(self.kind == "DH_over_rs", dh, dv)
        )


def capacities():
    rc = json.loads((HERE / "rc_load_source_only_certificate.json").read_text())[
        "inputs"
    ]
    p = float(rc["P"]["lo"])
    budget = 6 * pi / (p * float(rc["alpha_U"]["lo"]))
    x = (budget + sqrt(budget**2 - 24 * pi)) / 2
    other = json.loads(
        (
            ROOT / "code/capacity_readback/manifests/n_closure_branch_certificate.json"
        ).read_text()
    )["branches"]
    return dict(rc_load=pi * exp(x), **{d["branch_id"]: float(d["N"]) for d in other})


def fit(data, b, c, idx):
    chol = np.linalg.cholesky(data.cov[np.ix_(idx, idx)])

    def score(a):
        r = data.y[idx] - data.theory(a, b, c)[idx]
        white = solve_triangular(chol, r, lower=True)
        return float(white @ white)

    opt = minimize_scalar(
        lambda log_a: score(exp(log_a)),
        bounds=(-18, 0),
        method="bounded",
        options=dict(xatol=1e-12),
    )
    if not opt.success or not -17.9 < opt.x < -0.1:
        raise ValueError("Matter minimum is not resolved inside the search")
    return exp(opt.x), score


def conditional(data, a, b, c, train, hold, beta=0.01):
    """Condition correlated measurements; include uncertainty in the training calibration."""
    ct = data.cov[np.ix_(train, train)]
    cross = data.cov[np.ix_(hold, train)]
    transfer = np.linalg.solve(ct, cross.T).T
    ch = data.cov[np.ix_(hold, hold)] - transfer @ cross.T
    cholesky = np.linalg.cholesky(ch)
    transform = np.eye(len(data.z))[hold] - transfer @ np.eye(len(data.z))[train]
    white = solve_triangular(cholesky, transform, lower=True)
    f = data.theory(a, b, c)
    hres = transform @ (data.y - f)
    rs = white @ (data.y - f)
    chi_cond = rs @ rs

    # First-order forecast covariance with one GLS-estimated training coefficient.
    slope = (data.theory(a * 1.00001, b, c) - data.theory(a * 0.99999, b, c)) / (
        2e-5 * a
    )
    da = 1 / sqrt(slope[train] @ np.linalg.solve(ct, slope[train]))
    held_slope = transform @ slope
    predcov = ch + da**2 * np.outer(held_slope, held_slope)
    qpred = hres @ np.linalg.solve(predcov, hres)

    # Finite-size conservative cross-check. A positive weighted training sum is
    # Gaussian with known variance and a strictly decreasing expectation for all A>0.
    # Inverting its central interval covers the actual A with probability 1-beta.
    w = 1 / np.sqrt(ct.diagonal())
    observed = w @ data.y[train]
    margin = norm.isf(beta / 2) * sqrt(w @ ct @ w)

    def root(target):
        return brentq(
            lambda value: w @ data.theory(value, b, c)[train] - target,
            1e-14,
            1.0,
            xtol=1e-16,
        )

    low, high = root(observed + margin), root(observed - margin)

    # With the true A the conditional residual has exactly chi-square(h) norm.
    # Enclose its minimum over the training interval. Within each subdivision
    # every raw distance decreases, so the signed linear map can be bounded by
    # interval arithmetic even when a whitened coordinate is not monotone.
    panel = np.linspace(low, high, 257)
    fgrid = data.theory(panel, b, c)
    rgrid = (data.y - fgrid) @ white.T
    pos, neg = np.maximum(white, 0), np.minimum(white, 0)
    raw_lo, raw_hi = data.y - fgrid[:-1] - 1e-10, data.y - fgrid[1:] + 1e-10
    rlo = raw_lo @ pos.T + raw_hi @ neg.T
    rhi = raw_hi @ pos.T + raw_lo @ neg.T
    lower = np.sum(np.maximum(0.0, np.maximum(rlo, -rhi)) ** 2, axis=1).min()
    upper = np.sum(rgrid**2, axis=1).min()
    p_upper = min(1.0, beta + chi2.sf(lower, len(hold)))
    return dict(
        held_coordinates=[f"{data.z[i]} {data.kind[i]}" for i in hold],
        A_training=a,
        A_sigma_linear=da,
        A_training_interval_99=[low, high],
        conditional_residuals=hres.tolist(),
        predictive_sigmas_linear=np.sqrt(predcov.diagonal()).tolist(),
        fixed_A_chi2=float(chi_cond),
        predictive_chi2_linear=float(qpred),
        predictive_p_linear=chi2.sf(qpred, len(hold)),
        conservative_Q_min_bounds=[lower, upper],
        conservative_predictive_p_upper=p_upper,
    )


def repeated_anchor_tail(data, b, c, anchor, rest, a_true, threshold, order=128):
    """Tail of the plug-in residual when the high-redshift anchor is remeasured."""
    i = anchor[0]
    variance = data.cov[i, i]
    coeff = data.cov[rest, i] / variance
    cov = data.cov[np.ix_(rest, rest)] - variance * np.outer(coeff, coeff)
    chol = np.linalg.cholesky(cov)
    actual = data.theory(a_true, b, c)
    # Condition on a standard-normal anchor realization; its calibrated A and
    # predictive bias are then deterministic. Integrate out this one scalar.
    node, weight = np.polynomial.legendre.leggauss(order)
    zscore = 9 * node
    perturb = sqrt(variance) * zscore
    observed = actual[i] + perturb
    s = data.s[i]
    estimated = (observed**-2 - b - c * s**4) / s**3
    if np.min(observed) <= 0 or np.min(estimated) <= 0:
        raise ValueError("Unphysical calibrated anchor inside integration range")
    predicted = data.theory(estimated, b, c)
    bias = actual[rest] + perturb[:, None] * coeff - predicted[:, rest]
    delta = np.sum(solve_triangular(chol, bias.T, lower=True) ** 2, axis=0)
    probability = 9 * np.dot(
        weight, norm.pdf(zscore) * ncx2.sf(threshold, len(rest), delta)
    )
    return float(probability)


def summarize(data, n, ruler, length, omega_r):
    b = pi * (ruler * MPC_M / length) ** 2 / n if n else 0.0
    c = omega_r * (100 * ruler / C_KMS) ** 2
    idx = np.arange(len(data.z))
    full, score = fit(data, b, c, idx)
    anchor = np.flatnonzero((data.z == max(data.z)) & (data.kind == "DH_over_rs"))
    rest = np.setdiff1d(idx, anchor)
    z = data.s[anchor[0]]
    measured = data.y[anchor[0]]
    a_anchor = (measured**-2 - b - c * z**4) / z**3
    checks = {}
    for label, train, hold, a in [
        ("high_radial_forecasts_other_12", anchor, rest, a_anchor),
        ("other_12_forecast_high_radial", rest, anchor, fit(data, b, c, rest)[0]),
    ]:
        checks[label] = conditional(data, a, b, c, train, hold)
    for redshift in np.unique(data.z):
        hold = idx[data.z == redshift]
        train = idx[data.z != redshift]
        a_train, _ = fit(data, b, c, train)
        checks[f"omit_{redshift}"] = conditional(data, a_train, b, c, train, hold)
    at_full = data.theory(full, b, c)
    probe = checks["high_radial_forecasts_other_12"]
    threshold = probe["fixed_A_chi2"]
    interval = probe["A_training_interval_99"]
    trial_densities = [a_anchor, full, *np.linspace(*interval, 5)]
    integral = []
    for a in trial_densities:
        tail = repeated_anchor_tail(data, b, c, anchor, rest, a, threshold)
        refined = repeated_anchor_tail(data, b, c, anchor, rest, a, threshold, 192)
        if abs(tail - refined) > 1e-9:
            raise ValueError("Repeated-anchor integration did not converge")
        integral.append(dict(A_assumed_true=a, p_repeating_anchor=refined))
    return dict(
        B=b,
        C=c,
        A_full=full,
        A_anchor=a_anchor,
        chi2_full=score(full),
        chi2_at_anchor=score(a_anchor),
        nominal_profile_dof=12,
        nominal_profile_p=chi2.sf(score(full), 12),
        original_order=list(map(str, data.kind)),
        measurements=data.y.tolist(),
        forecast_at_anchor=data.theory(a_anchor, b, c).tolist(),
        forecast_at_profile=at_full.tolist(),
        omega_m_profile=full / (full + b + c),
        H0_profile=C_KMS * sqrt(full + b + c) / ruler,
        repeated_anchor_tail=integral,
        checks=checks,
    )


def horizon():
    receipt = json.loads((HERE.parent / "cut_signal.json").read_text())
    results = {}
    for label, strength in (("native", 8), ("exchange_only", 0)):
        rows = [
            x
            for x in receipt["records"]
            if x[2] is None and x[3] == strength and x[4] == 2
        ]
        pair = {(r[0], r[1]) for r in rows}
        if len(pair) != 50 or len(pair) != len(rows):
            raise ValueError("Missing enrolled crossing measurements")
        observed_endpoints = {int(key) for key in receipt["endpoints"]}
        if pair != {
            (edge, direction)
            for edge in observed_endpoints
            for direction in ("left", "right")
        }:
            raise ValueError("Crossing measurements do not match the enrolled boundary")
        duration = pi / 32
        # An active arm sees the even/odd written energies 4*kappa and 4*kappa+1.
        # Interference with a parked reference averages their two click probabilities.
        p_one = (
            2 - cos(4 * strength * duration) - cos((4 * strength + 1) * duration)
        ) / 4
        results[label] = dict(
            seams=len({r[0] for r in rows}),
            launches_per_bit=len(rows) * receipt["shots_per_bit"],
            zero_clicks=sum(r[5] for r in rows),
            one_clicks=sum(r[6] for r in rows),
            exact_zero_click_probability=0.0,
            ideal_one_click_probability=p_one,
        )
    return results


def source_status():
    direct = json.loads(
        (
            ROOT / "code/capacity_readback/runtime/direct_n_closure_verdict.json"
        ).read_text()
    )
    attachment = json.loads(
        (
            ROOT
            / "code/capacity_readback/runtime/horizon_record_attachment_verdict.json"
        ).read_text()
    )
    finite = direct["fixed_cutoff_result"]
    return dict(
        finite_dimension=finite["D"],
        measured_correctable_messages=finite["M0"],
        cosmic_value_selected=finite["cosmic_value_selected"],
        direct_source_verdict=direct["status"],
        horizon_verdict=attachment["status"],
    )


def run():
    data = Distances()
    caps = capacities()
    options = [
        ("rc_load", caps["rc_load"], R_DRAG, PLANCK_M, 0.0),
        ("rc_low_B", caps["rc_load"], R_DRAG - 0.26, PLANCK_M + 0.000018e-35, 0.0),
        ("rc_high_B", caps["rc_load"], R_DRAG + 0.26, PLANCK_M - 0.000018e-35, 0.0),
        ("rc_radiation", caps["rc_load"], R_DRAG, PLANCK_M, 4.18e-5),
        ("finite_presence", caps["finite_presence"], R_DRAG, PLANCK_M, 0.0),
        (
            "poisson_projective_limit",
            caps["poisson_projective_limit"],
            R_DRAG,
            PLANCK_M,
            0.0,
        ),
        ("dust_only", None, R_DRAG, PLANCK_M, 0.0),
    ]
    return dict(
        capacities_nats=caps,
        r_drag_Mpc=R_DRAG,
        horizon_crossing=horizon(),
        finite_source_attachment=source_status(),
        models={
            name: summarize(data, n, rd, lp, rad) for name, n, rd, lp, rad in options
        },
    )


def formatted(obj):
    if isinstance(obj, float):
        return float(format(obj, ".9g"))
    if isinstance(obj, dict):
        return {k: formatted(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [formatted(v) for v in obj]
    return obj


def main():
    opt = argparse.ArgumentParser(description=__doc__)
    opt.add_argument("--write", action="store_true")
    args = opt.parse_args()
    report = formatted(run())
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    else:

        def equal(a, b):
            if isinstance(a, dict):
                return a.keys() == b.keys() and all(
                    equal(v, b[k]) for k, v in a.items()
                )
            if isinstance(a, list):
                return len(a) == len(b) and all(equal(i, j) for i, j in zip(a, b))
            if isinstance(a, float):
                return np.isclose(a, b, rtol=1e-5, atol=2e-10)
            return a == b

        if not equal(report, json.loads(path.read_text())):
            raise ValueError("Recomputed result differs from receipt")
    rc = report["models"]["rc_load"]
    print(
        f"DESI chi² at high-redshift calibration {rc['chi2_at_anchor']:.3f}; "
        f"all-data density profile {rc['chi2_full']:.3f}. "
        f"Verified {len(report['models'])} models on seven whole-redshift omissions."
    )


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without -O to enable the scientific checks.")
    main()
