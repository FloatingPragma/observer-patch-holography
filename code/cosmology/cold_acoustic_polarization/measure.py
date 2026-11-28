"""Predict measured acoustic polarization after calibrating only on temperature."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_solve, cholesky, solve_triangular
from scipy.optimize import minimize

from .boltzmann import BARYONS, RATIOS, compute

HERE = Path(__file__).resolve().parent
CHANNELS = ("TT", "TE", "EE")
WINDOWS = [(120, 170), (270, 330), (400, 480), (730, 810), (1000, 1100)]


def whiten(chol, x):
    return solve_triangular(chol, x, lower=True, check_finite=False)


def norm(chol, x):
    s = whiten(chol, x)
    return float(s @ s)


def official_check(raw, cov):
    spec = raw["official_check_spectra"][[0, 2, 1]]  # clik input order TT,EE,TE
    parts = []
    for cl, n in zip(spec, (215, 199, 199), strict=True):
        vals = cl[30:2509] * raw["weight"]
        parts.append(np.add.reduceat(vals, raw["bin_start"])[:n])
    residual = (
        raw["observation"][:, 1]
        - np.concatenate(parts) / float(raw["official_check_calibration"]) ** 2
    )
    return -0.5 * norm(cholesky(cov, lower=True), residual)


class Likelihood:
    def __init__(self):
        with np.load(HERE / "plik_lite.npz") as raw:
            full = np.zeros((613, 613))
            full[np.tril_indices(613)] = raw["lower_triangle"]
            full = full + np.tril(full, -1).T
            self.check = official_check(raw, full)
            self.reference_check = float(raw["official_check_loglike"])
            if abs(self.check - self.reference_check) > 1e-5:
                raise ValueError(
                    f"Plik replay failed {self.check} vs {self.reference_check}"
                )
            lo, hi = 30 + raw["bin_start"], 30 + raw["bin_end"]
            index = np.where((lo >= 100) & (hi <= 1500))[0]
            self.low, self.high = lo[index], hi[index]
            self.ell = np.arange(self.low[0], self.high[-1] + 1)
            self.starts = self.low - self.low[0]
            self.weight = (
                raw["weight"][self.ell - 30] * 2 * np.pi / (self.ell * (self.ell + 1))
            )
            select = np.concatenate([index, index + 215, index + 414])
            self.y = raw["observation"][select, 1]
            self.eff = raw["observation"][index, 0]
            self.cov = full[np.ix_(select, select)]
        self.n = len(index)
        n = self.n
        self.T = cholesky(self.cov[:n, :n], lower=True)
        self.P = cholesky(self.cov[n:, n:], lower=True)
        self.B = cho_solve((self.T, True), self.cov[:n, n:]).T
        self.conditional_cov = self.cov[n:, n:] - self.B @ self.cov[:n, n:]
        self.cond = cholesky(self.conditional_cov, lower=True)
        self.log_l = np.log(self.ell / 500.0)
        self.temp = whiten(self.T, self.y[:n])

    def bands(self, spectrum, dilation, tilt):
        model = np.interp(self.ell / dilation, np.arange(len(spectrum)), spectrum)
        return np.add.reduceat(
            model * self.weight * np.exp(tilt * self.log_l), self.starts
        )

    def transformed(self, theory, pars):
        amplitude, dilation, tilt = pars
        return amplitude * np.concatenate(
            [self.bands(spec, dilation, tilt) for spec in theory.T]
        )

    def calibrate(self, tt, allow_tilt, limit):
        # Profile the linear amplitude exactly for each pair of nonlinear settings.
        def objective(x, amplitude=False):
            dilation = x[0]
            tilt = x[1] if allow_tilt else 0.0
            model = whiten(self.T, self.bands(tt, dilation, tilt))
            gain = max(0.0, (model @ self.temp) / (model @ model))
            diff = self.temp - gain * model
            return gain if amplitude else float(diff @ diff)

        bounds = [(0.55, 1.6)] + ([(-limit, limit)] if allow_tilt else [])
        # Acoustic zero crossings can create multiple dilation extrema. Refine the
        # best grid starting point plus every distinct one-dimensional local basin.
        ds = np.linspace(0.55, 1.6, 85)
        ts = (-limit, 0, limit) if allow_tilt else (0,)
        seeds = [
            min(([d, t] if allow_tilt else [d] for t in ts), key=objective) for d in ds
        ]
        val = np.array([objective(x) for x in seeds])
        candidates = [
            i
            for i in range(len(val))
            if (i == 0 or val[i] <= val[i - 1])
            and (i == len(val) - 1 or val[i] <= val[i + 1])
        ]
        fits = [
            minimize(
                objective,
                seeds[i],
                bounds=bounds,
                method="Nelder-Mead",
                options={"xatol": 1e-8, "fatol": 1e-6, "maxfev": 3000},
            )
            for i in candidates
        ]
        best = min(fits, key=lambda f: f.fun)
        if not best.success:
            raise RuntimeError(best.message)
        return [objective(best.x, True), best.x[0], best.x[1] if allow_tilt else 0.0]

    def check_polarization(self, model):
        n = self.n
        residual = self.y - model
        condres = residual[n:] - self.B @ residual[:n]
        result = {
            "TT": norm(self.T, residual[:n]),
            "joint_P": norm(self.P, residual[n:]),
            "conditional_P": norm(self.cond, condres),
        }
        for i, key in enumerate(("TE", "EE")):
            sl = slice(i * n, (i + 1) * n)
            block = slice((i + 1) * n, (i + 2) * n)
            result[key] = norm(
                cholesky(self.cov[block, block], lower=True), residual[block]
            )
            result["conditional_" + key] = norm(
                cholesky(self.conditional_cov[sl, sl], lower=True), condres[sl]
            )
        return result

    def summary_bands(self, vector, errors=False):
        result = {}
        fac = self.eff * (self.eff + 1) / (2 * np.pi)
        for i, ch in enumerate(CHANNELS[1:], start=1):
            block = slice(i * self.n, (i + 1) * self.n)
            d = vector[block] * fac
            var = np.diag(self.cov)[block] * fac**2
            result[ch] = {}
            for lo, hi in WINDOWS:
                selected = (self.low >= lo) & (self.high <= hi)
                w = 1 / var[selected]
                w = w / w.sum()
                if errors:
                    sub = self.cov[block, block][np.ix_(selected, selected)]
                    weights = w * fac[selected]
                    val = np.sqrt(weights @ sub @ weights)
                else:
                    val = w @ d[selected]
                result[ch][f"{lo}-{hi}"] = float(val)
        return result

    def evaluate(self, theory, allow_tilt, limit=0.2):
        pars = self.calibrate(theory[:, 0], allow_tilt, limit)
        fitted = self.transformed(theory, pars)
        result = self.check_polarization(fitted)
        n = self.n
        tolerances = [1e-4, 1e-5] + ([1e-4] if allow_tilt else [])
        jac = []
        for i, step in enumerate(tolerances):
            plus = np.array(pars)
            minus = np.array(pars)
            plus[i] += step
            minus[i] -= step
            jac.append(
                (self.transformed(theory, plus) - self.transformed(theory, minus))
                / (2 * step)
            )
        J = np.column_stack(jac)
        tJ = whiten(self.T, J[:n])
        nuisance_cov = np.linalg.inv(tJ.T @ tJ)
        Jp = J[n:] - self.B @ J[:n]
        pred_cov = self.conditional_cov + Jp @ nuisance_cov @ Jp.T
        res = self.y - fitted
        result["conditional_P_with_TT_parameter_errors"] = norm(
            cholesky(pred_cov, lower=True), res[n:] - self.B @ res[:n]
        )
        efficiency = {}
        for e in (0.95, 1.0, 1.05):
            modified = fitted.copy()
            modified[n : 2 * n] *= e
            modified[2 * n :] *= e**2
            efficiency[str(e)] = self.check_polarization(modified)["conditional_P"]
        return {
            "temperature_calibration": pars,
            "local_TT_parameter_sd": np.sqrt(np.diag(nuisance_cov)).tolist(),
            "nuisance_at_limit": bool(
                abs(pars[2]) > limit - 1e-5 or min(pars[1] - 0.55, 1.6 - pars[1]) < 1e-5
            ),
            "scores": result,
            "conditional_score_by_polarization_efficiency": efficiency,
            "predicted_windows_uK2": self.summary_bands(fitted),
        }


def run(cache=None):
    import camb

    lk = Likelihood()
    models = {}
    plot_predictions = {}
    fac = np.tile(lk.eff * (lk.eff + 1) / (2 * np.pi), 3)
    for b in BARYONS:
        models[str(b)] = {}
        for label, ratio in RATIOS.items():
            if (
                cache is not None
                and (path := cache / f"model_{b}_{label}.npz").exists()
            ):
                with np.load(path) as f:
                    theory = f["theory"]
                    zeq = float(f["zeq"])
            else:
                theory, zeq = compute(b, ratio)
            models[str(b)][label] = {
                "ratio": ratio,
                "z_equality": zeq,
                "amplitude_scale_tilt": lk.evaluate(theory, True),
                "amplitude_scale": lk.evaluate(theory, False),
                "broad_tilt_sensitivity": lk.evaluate(theory, True, 1.0),
            }
            if b == BARYONS[0] and label in (
                "federation_matching",
                "max_service",
                "standard_control",
            ):
                pars = models[str(b)][label]["amplitude_scale_tilt"][
                    "temperature_calibration"
                ]
                plot_predictions[label] = (
                    (lk.transformed(theory, pars) * fac).reshape(3, -1).tolist()
                )
    bands = {
        "ell": lk.eff.tolist(),
        "low": lk.low.tolist(),
        "high": lk.high.tolist(),
        "measured": (lk.y * fac).reshape(3, -1).tolist(),
        "sd": (np.sqrt(np.diag(lk.cov)) * fac).reshape(3, -1).tolist(),
        "predictions": plot_predictions,
    }
    return {
        "selected_bands_uK2": bands,
        "CAMB": camb.__version__,
        "n_per_channel": lk.n,
        "full_multipoles": [int(lk.low[0]), int(lk.high[-1])],
        "official_loglike_check": [lk.check, lk.reference_check],
        "measured_windows_uK2": lk.summary_bands(lk.y),
        "measured_window_sd_uK2": lk.summary_bands(lk.y, errors=True),
        "models": models,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-cache", type=Path)
    args = parser.parse_args()
    (HERE / "results.json").write_text(
        json.dumps(run(args.model_cache), indent=2, allow_nan=False) + "\n"
    )
