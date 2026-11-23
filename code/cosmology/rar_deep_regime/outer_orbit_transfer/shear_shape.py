"""Six-bin shear residuals for the actual orbit witnesses and a free spherical scale."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import lsq_linear
from scipy.stats import chi2

from .. import lensing_profile

HERE = Path(__file__).resolve().parent


def shear(acceleration, amplitude):
    return (
        4 / math.pi * acceleration + amplitude * np.sqrt(acceleration)
    ) / lensing_profile.FOUR_G


def six_bins(profile):
    low = np.asarray(profile["lower_bin_edge_m_s2"])
    selected = low >= 1e-13
    centre = np.asarray(profile["g_bar_m_s2"])[selected]
    low = low[selected]
    high = 2 * centre - low
    observed = np.asarray(profile["tangential_esd"])[selected]
    covariance = np.asarray(profile["covariance"])[np.ix_(selected, selected)]
    if len(centre) != 6 or np.any(low <= 0) or np.any(centre <= low):
        raise ValueError("Expected six positive isolated acceleration bins")
    return low, centre, high, observed, covariance


def fixed_amplitude(low, centre, high, observed, covariance, amplitude):
    # In this metric Euclidean residual length squared equals r.T @ inv(C) @ r.
    whitening = solve_triangular(
        np.linalg.cholesky(covariance), np.eye(len(centre)), lower=True
    )
    predicted = shear(centre, amplitude)
    residual = observed - predicted
    floor, ceiling = shear(low, amplitude), shear(high, amplitude)
    fit = lsq_linear(
        whitening,
        whitening @ observed,
        bounds=(floor, ceiling),
        method="bvls",
        tol=1e-13,
    )
    if not fit.success:
        raise RuntimeError(fit.message)
    return dict(
        a0_m_s2=float(amplitude**2),
        centres=dict(
            predicted_msun_psqpc=predicted.tolist(),
            observed_minus_predicted_msun_psqpc=residual.tolist(),
            marginal_residual_in_errors=(
                residual / np.sqrt(covariance.diagonal())
            ).tolist(),
            chi_square=float(np.linalg.norm(whitening @ residual) ** 2),
        ),
        arbitrary_within_bin_weights=dict(
            min_predicted_msun_psqpc=floor.tolist(),
            max_predicted_msun_psqpc=ceiling.tolist(),
            optimizing_mean_msun_psqpc=fit.x.tolist(),
            active_lower_or_upper=fit.active_mask.tolist(),
            min_squared_covariance_distance=float(
                np.linalg.norm(whitening @ (observed - fit.x)) ** 2
            ),
        ),
    )


def free_scale(centre, observed, covariance):
    lower_cholesky = np.linalg.cholesky(covariance)
    target = solve_triangular(lower_cholesky, observed - shear(centre, 0), lower=True)
    slope = solve_triangular(
        lower_cholesky, np.sqrt(centre) / lensing_profile.FOUR_G, lower=True
    )
    amplitude = max(0.0, float(slope @ target / (slope @ slope)))
    deviance = float(np.linalg.norm(target - amplitude * slope) ** 2)
    degrees = len(centre) - 1
    tail = float(chi2.sf(deviance, degrees)) if amplitude > 0 else None
    return dict(
        a0_m_s2=amplitude**2,
        amplitude_se=1 / math.sqrt(slope @ slope),
        chi_square=deviance,
        residual_degrees_of_freedom=degrees,
        gaussian_covariance_only_tail=tail,
    )


def run():
    linear_path = HERE / "linear.json"
    witnesses = json.loads(linear_path.read_bytes())["panels"]
    profiles = json.loads(lensing_profile.RECEIPT.read_bytes())["profiles"]
    surveys = {}
    for name in ("KiDS-bright", "GAMA"):
        low, centre, high, observed, covariance = six_bins(profiles[name])
        panels = {}
        for panel, choices in witnesses.items():
            panels[panel] = {}
            for choice, result in choices.items():
                if result["feasible"]:
                    amplitude = math.sqrt(result["a0_m_s2"])
                    panels[panel][choice] = fixed_amplitude(
                        low, centre, high, observed, covariance, amplitude
                    )
        surveys[name] = dict(
            central_bin_best_fit=free_scale(centre, observed, covariance),
            observed_msun_psqpc=observed.tolist(),
            witnesses=panels,
        )
    return dict(
        linear_witness_sha256=hashlib.sha256(linear_path.read_bytes()).hexdigest(),
        lens_profile_sha256=hashlib.sha256(
            lensing_profile.RECEIPT.read_bytes()
        ).hexdigest(),
        surveys=surveys,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as out:
        json.dump(run(), out, indent=2, allow_nan=False)
        out.write("\n")
