"""Test the fixed RC-LOAD horizon amplitude against DESI DR2 distance records."""

import argparse
import hashlib
import json
from math import exp, pi, sin, sqrt
from pathlib import Path
from statistics import NormalDist

import numpy as np
from scipy.constants import parsec
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2, norm

from cosmology.paid_trail_bao.measure import (
    DATA,
    HERE as DATA_DIR,
    PREFIX,
    DistanceData,
)

HERE = Path(__file__).resolve().parent
# Midpoints rounded from the archived input intervals. The declared macroscopic
# domain selects the larger quadratic root without using a distance measurement.
P = 1.6309720958588974
ALPHA_U = 0.041124247441816685
PRICE_BUDGET = 6 * pi / (P * ALPHA_U)
X = (PRICE_BUDGET + sqrt(PRICE_BUDGET**2 - 24 * pi)) / 2
L_PLANCK = 1.616255e-35  # CODATA 2022, metres
L_SIGMA = 0.000018e-35
R_DRAG = 147.09  # Planck TT,TE,EE+lowE+lensing base LCDM, Mpc
R_SIGMA = 0.26
MPC = 1e6 * parsec


def distance_vector(data, matter, vacuum, radiation=0.0):
    """A, B, C multiply (1+z)^3, 1, (1+z)^4 in (H*r_d/c)^2."""

    def radial(z):
        s = 1 + z
        return 1 / sqrt(matter * s**3 + vacuum + radiation * s**4)

    unique, ix = np.unique(data.z, return_inverse=True)
    dm = np.array([quad(radial, 0, z, epsabs=1e-11)[0] for z in unique])[ix]
    dh = np.array([radial(z) for z in data.z])
    dv = np.cbrt(data.z * dm**2 * dh)
    readings = {"DM_over_rs": dm, "DH_over_rs": dh, "DV_over_rs": dv}
    return np.array([readings[k][i] for i, k in enumerate(data.kind)])


def dust_contrasts(data):
    # For any positive A, integrating the dust Hubble law gives
    # D_M / D_H = 2 (1+z) (sqrt(1+z)-1). Work in the original Gaussian
    # distance coordinates, avoiding a Gaussian approximation to their ratio.
    rows = []
    contrasts = []
    for z in data.z[data.kind == "DH_over_rs"]:
        s = 1 + z
        factor = 2 * s * (sqrt(s) - 1)
        v = np.zeros(len(data.z))
        v[(data.z == z) & (data.kind == "DM_over_rs")] = 1
        v[(data.z == z) & (data.kind == "DH_over_rs")] = -factor
        rows.append(v)
        contrasts.append(dict(z=z, predicted_DM_over_DH=factor))
    design = np.array(rows)
    diff = design @ data.obs
    cov = design @ data.cov @ design.T
    statistic = float(diff @ np.linalg.solve(cov, diff))
    for item, value, error in zip(contrasts, diff, np.sqrt(np.diag(cov))):
        item.update(DM_minus_factor_DH=value, sigma=error)
    return dict(
        rows=contrasts,
        chi2=statistic,
        dof=len(diff),
        gaussian_tail=chi2.sf(statistic, len(diff)),
    )


def capacity_B(ruler, length):
    return exp(-X) * (ruler * MPC / length) ** 2


def nonacceleration_test(data):
    """Necessary inequalities for any always-decelerating or coasting flat FLRW history."""
    units = np.eye(len(data.z))
    s = 1 + data.z
    x = np.log1p(data.z)
    (volume,) = np.flatnonzero(data.kind == "DV_over_rs")
    radial = {z: i for i, z in enumerate(data.z) if data.kind[i] == "DH_over_rs"}
    transverse = {z: i for i, z in enumerate(data.z) if data.kind[i] == "DM_over_rs"}
    rows = []

    def check(kind, z1, z2, weights):
        value = weights @ data.obs
        error = sqrt(weights @ data.cov @ weights)
        rows.append(
            dict(
                kind=kind,
                z1=z1,
                z2=z2,
                contrast=value,
                sigma=error,
                standard_errors=value / error,
            )
        )

    # In x = ln(1+z), the slope of D_M/r_d is u = (1+z)D_H/r_d.
    # For q>=0, u is nonincreasing, so any earlier average bounds a later u.
    volume_scale = (data.z[volume] * x[volume] ** 2 / s[volume]) ** (1 / 3)
    for z2, j in radial.items():
        u2 = s[j] * units[j]
        if data.z[volume] <= z2:
            check(
                "volume_to_later_slope",
                data.z[volume],
                z2,
                units[volume] / volume_scale - u2,
            )
        for z1, i in transverse.items():
            if z1 <= z2:
                check("mean_to_later_slope", z1, z2, units[i] / x[i] - u2)
            if z1 < z2:
                u1 = s[radial[z1]] * units[radial[z1]]
                chord = (units[transverse[z2]] - units[i]) / (x[j] - x[i])
                check("slope_order", z1, z2, u1 - u2)
                check("upper_chord", z1, z2, u1 - chord)
                check("lower_chord", z1, z2, chord - u2)
    minimum = min(rows, key=lambda row: row["standard_errors"])
    # Gaussian marginals and the union bound suffice regardless of their correlations.
    return dict(
        rows=rows,
        strongest=minimum,
        family_size=len(rows),
        gaussian_family_p_upper=min(
            1.0, len(rows) * norm.cdf(minimum["standard_errors"])
        ),
    )


def joint_matter_box(data, vacuum, radiation, k, anchor_interval):
    """Invert monotone predictions; return a single shared A interval in all 13 bands."""
    lower, upper = anchor_interval
    sigma = np.sqrt(np.diag(data.cov))
    distance_max = distance_vector(data, lower, vacuum, radiation)
    distance_min = distance_vector(data, upper, vacuum, radiation)
    surviving = list(anchor_interval)
    for i, (measured, sd) in enumerate(zip(data.obs, sigma)):
        low, high = measured - k * sd, measured + k * sd
        if distance_min[i] > high or distance_max[i] < low:
            return None

        def difference(a, target):
            return distance_vector(data, a, vacuum, radiation)[i] - target

        if distance_max[i] > high:
            surviving[0] = max(
                surviving[0], brentq(difference, lower, upper, args=(high,), xtol=1e-16)
            )
        if distance_min[i] < low:
            surviving[1] = min(
                surviving[1], brentq(difference, lower, upper, args=(low,), xtol=1e-16)
            )
    return surviving if surviving[0] <= surviving[1] else None


def curved_nonacceleration_test(data):
    """One-sided, simultaneous comparison allowing arbitrary FLRW spatial curvature."""
    (volume,) = np.flatnonzero(data.kind == "DV_over_rs")
    transverse = np.flatnonzero(data.kind == "DM_over_rs")
    far = transverse[np.argmax(data.z[transverse])]
    z = data.z[volume]
    sigma = np.sqrt(np.diag(data.cov))
    k = norm.isf(0.05 / len(data.obs))
    v_upper = data.obs[volume] + k * sigma[volume]
    far_lower = data.obs[far] - k * sigma[far]
    rows = []
    for j in np.flatnonzero(data.kind == "DH_over_rs"):
        zb = data.z[j]
        if not z <= zb <= data.z[far]:
            continue
        u = (1 + zb) * (data.obs[j] - k * sigma[j])
        radial_min = np.log1p(z) * u
        # On a direct ray prior to the first conjugate point, a decreasing
        # D_M at the near object already exceeds the far one. On the rising
        # branch R sin(chi/R) increases with chi then increases with R at fixed chi.
        if not 0 < radial_min < pi * far_lower / 2:
            raise ValueError("Rising-branch lower-bound domain not certified")
        transverse_min = far_lower * sin(radial_min / far_lower)
        v_min = (z * transverse_min**2 * u / (1 + z)) ** (1 / 3)
        rows.append(
            dict(
                z_radial=zb,
                u_lower=u,
                near_comoving_radial_lower=radial_min,
                near_transverse_lower=transverse_min,
                nonaccelerating_DV_lower=v_min,
                gap=v_min - v_upper,
            )
        )
    return dict(
        family_one_sided_k=k,
        confidence_at_least=0.95,
        z_volume=z,
        measured_DV_upper=v_upper,
        z_far=data.z[far],
        measured_far_DM_lower=far_lower,
        rows=rows,
    )


def prediction(data, vacuum, radiation=0.0):
    anchor = int(
        np.flatnonzero((data.z == max(data.z)) & (data.kind == "DH_over_rs")).item()
    )
    sigma = np.sqrt(np.diag(data.cov))
    s = 1 + data.z[anchor]

    def matter_from_anchor(dh):
        return (dh**-2 - vacuum - radiation * s**4) / s**3

    matter = matter_from_anchor(data.obs[anchor])
    forecast = distance_vector(data, matter, vacuum, radiation)
    # Simultaneous coverage at least 95% for 13 published Gaussian
    # marginals by the union bound, including the calibration measurement.
    k = NormalDist().inv_cdf(1 - 0.05 / (2 * len(data.obs)))
    lower = distance_vector(
        data,
        matter_from_anchor(data.obs[anchor] - k * sigma[anchor]),
        vacuum,
        radiation,
    )
    upper = distance_vector(
        data,
        matter_from_anchor(data.obs[anchor] + k * sigma[anchor]),
        vacuum,
        radiation,
    )
    bounds_met = (upper >= data.obs - k * sigma) & (lower <= data.obs + k * sigma)
    anchor_interval = (
        matter_from_anchor(data.obs[anchor] + k * sigma[anchor]),
        matter_from_anchor(data.obs[anchor] - k * sigma[anchor]),
    )
    joint_interval = joint_matter_box(data, vacuum, radiation, k, anchor_interval)

    # Profile ordinary + additional dust density as a nuisance on all points
    # in a separate check. Optimize in log A, far beyond the observational range.
    def objective(log_a):
        res = data.obs - distance_vector(data, exp(log_a), vacuum, radiation)
        return float(res @ np.linalg.solve(data.cov, res))

    fit = minimize_scalar(
        objective, bounds=(-18, 0), method="bounded", options={"xatol": 1e-12}
    )
    if not fit.success or abs(fit.x + 18) < 1e-5 or abs(fit.x) < 1e-5:
        raise ValueError("Matter-profile minimum unresolved")
    records = [
        dict(
            z=z,
            kind=kind,
            observed=obs,
            sigma=sd,
            forecast=val,
            anchor_interval_forecast=[lo, hi],
            intervals_overlap=bool(ok),
        )
        for z, kind, obs, sd, val, lo, hi, ok in zip(
            data.z, data.kind, data.obs, sigma, forecast, lower, upper, bounds_met
        )
    ]
    return dict(
        B=vacuum,
        C=radiation,
        anchor_row=anchor,
        A_from_anchor=matter,
        omega_m_from_anchor=matter / (matter + vacuum + radiation),
        q0_from_anchor=(matter / 2 + radiation - vacuum)
        / (matter + vacuum + radiation),
        rows=records,
        simultaneous_box_k=k,
        disjoint_held_rows=int(sum(~bounds_met)),
        shared_A_in_all_13_simultaneous_bands=joint_interval,
        gaussian_13D_95_threshold=chi2.ppf(0.95, len(forecast)),
        all_data_A=exp(fit.x),
        all_data_profile_chi2=fit.fun,
        nominal_dof=len(forecast) - 1,
        nominal_tail=chi2.sf(fit.fun, len(forecast) - 1),
    )


def run():
    data = DistanceData()
    vacuum = capacity_B(R_DRAG, L_PLANCK)
    fixed = prediction(data, vacuum)
    # Sensitivities, not joint statistical intervals: extremes of the external
    # ruler and length at one quoted error, and an illustrative radiation term.
    low = capacity_B(R_DRAG - R_SIGMA, L_PLANCK + L_SIGMA)
    high = capacity_B(R_DRAG + R_SIGMA, L_PLANCK - L_SIGMA)
    radiation = 4.18e-5 * (100 * R_DRAG / 299792.458) ** 2
    sensitivity = {}
    for label, b, c in (
        ("low_B", low, 0),
        ("high_B", high, 0),
        ("physical_omega_r_4.18e-5", vacuum, radiation),
    ):
        check = prediction(data, b, c)
        sensitivity[label] = {
            k: check[k]
            for k in (
                "B",
                "C",
                "all_data_profile_chi2",
                "q0_from_anchor",
                "disjoint_held_rows",
            )
        }
    return dict(
        data=DATA,
        hashes={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [DATA_DIR / (PREFIX + s + ".txt") for s in ("mean", "cov")]
        },
        archive="374b2390801842d38de9b9c326e79a8ab29a6be5",
        X=X,
        N=pi * exp(X),
        lambda_lstar_squared=3 * exp(-X),
        lstar_assumed_m=L_PLANCK,
        ruler_assumed_Mpc=R_DRAG,
        unconditional_flat_dust_shape=dust_contrasts(data),
        flat_nonaccelerating_history=nonacceleration_test(data),
        arbitrary_curvature_nonaccelerating_history=curved_nonacceleration_test(data),
        fixed_capacity=fixed,
        zero_vacuum=prediction(data, 0),
        sensitivity=sensitivity,
    )


def verify(actual, archived):
    if isinstance(actual, dict):
        return actual.keys() == archived.keys() and all(
            verify(value, archived[key]) for key, value in actual.items()
        )
    if isinstance(actual, list):
        return len(actual) == len(archived) and all(map(verify, actual, archived))
    if isinstance(actual, float):
        return bool(np.isclose(actual, archived, rtol=2e-7, atol=0.0))
    return actual == archived


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    result = run()
    target = HERE / "capacity_distances.json"

    def rounded(x):
        if isinstance(x, float):
            return float(format(x, ".10g"))
        if isinstance(x, dict):
            return {k: rounded(v) for k, v in x.items()}
        if isinstance(x, list):
            return [rounded(v) for v in x]
        return x

    serialized = rounded(result)
    if opt.write:
        target.write_text(json.dumps(serialized, indent=2, allow_nan=False) + "\n")
    elif not verify(serialized, json.loads(target.read_text())):
        raise ValueError("Capacity distance report differs; inspect inputs and output")
    print("DESI DR2 dust contrasts and fixed capacity replayed")


if __name__ == "__main__":
    main()
