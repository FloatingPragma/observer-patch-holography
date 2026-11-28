"""Compare each physical drag ruler and flat background with DESI DR2 BAO."""

import hashlib
import json
from functools import lru_cache
from pathlib import Path

import camb
import numpy as np
from scipy.linalg import cholesky, cho_solve, solve_triangular
from scipy.optimize import brentq, minimize_scalar

from cosmology.cold_acoustic_polarization.boltzmann import BARYONS, RATIOS
from cosmology.cold_acoustic_polarization.trails import RECEIPT, prepare

HERE = Path(__file__).resolve().parent
PREFIX = "desi_gaussian_bao_ALL_GCcomb_"
C = 299792.458
REV = "b7b8a36e9bccb063081f811f323cada21ab5fbdd"
DATA = f"https://github.com/CobayaSampler/bao_data/tree/{REV}/desi_bao_dr2"


def whiten(chol, vec):
    return solve_triangular(chol, vec, lower=True, check_finite=False)


def quadratic(chol, vec):
    w = whiten(chol, vec)
    return float(w @ w)


class DistanceData:
    def __init__(self):
        table = np.loadtxt(HERE / f"{PREFIX}mean.txt", dtype=str, skiprows=1)
        self.z, self.obs = table[:, :2].T.astype(float)
        self.kind = table[:, 2]
        self.cov = np.loadtxt(HERE / f"{PREFIX}cov.txt")
        self.full = cholesky(self.cov, lower=True)
        self.train = self.z <= 1
        self.hold = ~self.train
        tt = self.cov[np.ix_(self.train, self.train)]
        ht = self.cov[np.ix_(self.hold, self.train)]
        hh = self.cov[np.ix_(self.hold, self.hold)]
        self.low = cholesky(tt, lower=True)
        self.B = cho_solve((self.low, True), ht.T).T
        self.conditional_cov = hh - self.B @ ht.T
        self.high = cholesky(self.conditional_cov, lower=True)

    def score(self, pred):
        res = self.obs - pred
        train = res[self.train]
        return {
            "full": quadratic(self.full, res),
            "train": quadratic(self.low, train),
            "held_given_train": quadratic(self.high, res[self.hold] - self.B @ train),
        }

    def free_gain(self, pred):
        shape = whiten(self.full, pred)
        return float((shape @ whiten(self.full, self.obs)) / (shape @ shape))


def physical_background(data, baryon, ratio, h0):
    p = camb.CAMBparams()
    helium = float(np.asarray(camb.bbn.get_predictor().Y_He(baryon, 0.0)).item())
    p.set_cosmology(
        H0=h0,
        ombh2=baryon,
        omch2=baryon * ratio,
        mnu=0.06,
        omk=0,
        tau=0.055,
        YHe=helium,
    )
    bg = camb.get_background(p)
    zs, indices = np.unique(data.z, return_inverse=True)
    dm = ((1 + zs) * bg.angular_diameter_distance(zs))[indices]
    dh = (C / bg.hubble_parameter(zs))[indices]
    distances = {
        "DM_over_rs": dm,
        "DH_over_rs": dh,
        "DV_over_rs": np.cbrt(data.z * dm**2 * dh),
    }
    vec = np.array([distances[k][i] for i, k in enumerate(data.kind)])
    params = bg.get_derived_params()
    rd = float(params["rdrag"])
    return {
        "H0": float(h0),
        "rd_Mpc": rd,
        "zdrag": float(params["zdrag"]),
        "Omega_m": sum(bg.get_Omega(x) for x in ("baryon", "cdm", "nu")),
        "Omega_de": bg.get_Omega("de"),
        "rdrag_times_h_Mpc": rd * h0 / 100,
    }, vec / rd


def fit(objective, bounds):
    result = minimize_scalar(
        objective, method="bounded", bounds=bounds, options={"xatol": 1e-5}
    )
    if not result.success:
        raise RuntimeError(result.message)
    return float(result.x)


def profile_errors(objective, h0, bounds):
    goal = objective(h0) + 1
    return [
        brentq(lambda h: objective(h) - goal, a, b, xtol=1e-5)
        for a, b in [(bounds[0], h0), (h0, bounds[1])]
    ]


def evaluate(data, baryon, board):
    @lru_cache(maxsize=None)
    def model(h0):
        return physical_background(data, baryon, board["ratio"], h0)

    def summary(h0, free=False):
        bg, vec = model(h0)
        gain = data.free_gain(vec) if free else 1.0
        return dict(
            **bg,
            gain=gain,
            effective_ruler_Mpc=bg["rd_Mpc"] / gain,
            predicted=(gain * vec).tolist(),
            scores=data.score(gain * vec),
        )

    ref = summary(70.0)
    h_min = max(35.0, 70 * np.sqrt(1 - ref["Omega_de"]) + 0.01)
    bounds = (h_min, 110.0)
    fixed = {str(h): summary(h) for h in (67.4, 70.0, 73.0)}
    profiles = {}
    for name, key, free in (
        ("all_H0", "full", False),
        ("low_z_H0", "train", False),
        ("all_H0_and_free_ruler", "full", True),
    ):

        def objective(h0):
            if key == "train":
                _, vector = model(h0)
                return quadratic(data.low, data.obs[data.train] - vector[data.train])
            return summary(h0, free)["scores"][key]

        h0 = fit(objective, bounds)
        best = summary(h0, free)
        interval = profile_errors(objective, h0, bounds)
        best["H0_delta_chi2_1_interval"] = interval
        if key == "train":
            best["held_scores_at_train_H0_endpoints"] = [
                summary(h)["scores"]["held_given_train"] for h in interval
            ]
        profiles[name] = best
    return {
        "fixed_H0_physical_ruler": fixed,
        "fixed_H0_70_free_ruler": summary(70.0, True),
        **profiles,
    }


def run():
    data = DistanceData()
    boards = [prepare(b) for b in json.loads(RECEIPT.read_text())["boards"]]
    stocks = {str(b["owners"]): b for b in boards}
    stocks["standard_control"] = {"ratio": RATIOS["standard_control"]}
    for stock in stocks.values():
        stock["runs"] = {str(b): evaluate(data, b, stock) for b in BARYONS}
        stock["BBN_calibration"] = {
            str(b): evaluate(data, b, stock) for b in (0.02218, 0.02163, 0.02273)
        }
    return {
        "CAMB": camb.__version__,
        "source_native_sha256": hashlib.sha256(RECEIPT.read_bytes()).hexdigest(),
        "source_DESI": DATA,
        "DESI_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(HERE.glob(PREFIX + "*.txt"))
        },
        "z": data.z.tolist(),
        "kind": data.kind.tolist(),
        "observed": data.obs.tolist(),
        "sigma": np.sqrt(np.diag(data.cov)).tolist(),
        "max_abs_released_cross_split_covariance": np.max(
            np.abs(data.cov[np.ix_(data.train, data.hold)])
        ),
        "stocks": stocks,
    }


if __name__ == "__main__":
    (HERE / "results.json").write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
