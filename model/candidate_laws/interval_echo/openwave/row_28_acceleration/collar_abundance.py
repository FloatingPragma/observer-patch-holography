"""Spend the spherical source's independently measured mass ratio in the BAO test."""

import argparse
import json
from fractions import Fraction
from itertools import combinations
from statistics import NormalDist
from math import exp
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo_partition_repair.published_strings.archive import (
    posted_words,
)
from .capacity_distances import (
    DistanceData,
    L_PLANCK,
    L_SIGMA,
    MPC,
    R_DRAG,
    R_SIGMA,
    X,
    capacity_B,
    distance_vector,
    verify,
)

HERE = Path(__file__).resolve().parent
C_LIGHT = 299792.458  # km/s
OMEGA_B = 0.02218  # DESI DR2 II Eq. 14, BBN external calibration
BBN_ERROR = 0.00055


def omega_vac(length):
    # H_100 in km/s/Mpc; multiply the inverse horizon length by c/H_100.
    return exp(-X) * (MPC * (C_LIGHT / 100) / length) ** 2


def required(dh, z, ruler, length, baryon, neutrino, radiation):
    hubble2 = (C_LIGHT / (100 * ruler * dh)) ** 2
    total_dust = (hubble2 - omega_vac(length) - radiation * (1 + z) ** 4) / (1 + z) ** 3
    return (total_dust - neutrino) / baryon - 1


def curved_background_floor(
    data, ruler_upper, baryon_upper, neutrino_upper, radiation_upper
):
    """Eliminate arbitrary constant spatial curvature and vacuum in all radial triples."""
    radial = data.kind == "DH_over_rs"
    redshifts = data.z[radial]
    observed = data.obs[radial]
    sigma = np.sqrt(np.diag(data.cov)[radial])
    # Protect the whole six-coordinate family before selecting the sharpest triple.
    k = NormalDist().inv_cdf(1 - 0.04 / (2 * len(redshifts)))
    lo = (observed + k * sigma) ** -2
    hi = (observed - k * sigma) ** -2
    scale_upper = (100 * ruler_upper / C_LIGHT) ** 2
    checks = []
    for choice in combinations(range(len(redshifts)), 3):
        idx = list(choice)
        t = 1 + redshifts[idx]
        # Sum w*(H*r_d/c)^2 eliminates both the constant and curvature.
        s = [1 + Fraction(str(z)) for z in redshifts[idx]]
        cross = [s[(i + 1) % 3] ** 2 - s[(i + 2) % 3] ** 2 for i in range(3)]
        scale = sum(x**3 * c for x, c in zip(s, cross))
        exact = [c / scale for c in cross]
        assert [
            sum(c * x**degree for x, c in zip(s, exact)) for degree in (0, 2, 3)
        ] == [0, 0, 1]
        w = np.array([float(c) for c in exact])
        dimensionless_min = w @ np.where(w > 0, lo[idx], hi[idx])
        # A nonpositive numerator alone supplies no bound beyond nonnegative mass.
        min_dust = (
            max(
                0.0,
                dimensionless_min / scale_upper - max(0.0, w @ t**4) * radiation_upper,
            )
            if dimensionless_min > 0
            else 0.0
        )
        checks.append(
            dict(
                z=redshifts[idx].tolist(),
                weights=list(map(str, exact)),
                physical_dust_lower=min_dust,
                extra_per_baryon_lower=(min_dust - neutrino_upper) / baryon_upper - 1,
            )
        )
    return dict(
        total_triples=len(list(combinations(redshifts, 3))),
        k=k,
        simultaneous_radial_coverage_at_least=0.96,
        triples=checks,
        strongest=max(checks, key=lambda item: item["extra_per_baryon_lower"]),
    )


def run():
    words = posted_words()
    traveller_count = sum(left for left, _ in words)
    coincidence_count = sum(left * right for left, right in words)
    # Read mass weights from the declared one-mass-unit-per-accepted-post
    # source construction. A coincidence is a subset of traveller posts.
    measured = Fraction(coincidence_count, traveller_count)
    data = DistanceData()
    idx = int(
        np.flatnonzero((data.z == max(data.z)) & (data.kind == "DH_over_rs")).item()
    )
    z, dh, sd = data.z[idx], data.obs[idx], np.sqrt(data.cov[idx, idx])
    nominal_nu, nominal_r = 0.00064, 0.000036  # declared late-time background
    scale = (100 * R_DRAG / C_LIGHT) ** 2
    cases = {}
    for name, ratio in (
        ("observed_pair_frequency", float(measured)),
        ("maximal_coincidence", 1.0),
    ):
        pred = distance_vector(
            data,
            (OMEGA_B * (1 + ratio) + nominal_nu) * scale,
            capacity_B(R_DRAG, L_PLANCK),
            nominal_r * scale,
        )
        err = data.obs - pred
        cases[name] = dict(
            ratio=ratio,
            prediction=pred.tolist(),
            all_data_chi2=float(err @ np.linalg.solve(data.cov, err)),
            predicted_radial_at_high_z=pred[idx],
        )

    # Monotonic directions follow from H^2 > 0 and a positive required
    # excess throughout this box. Enlarge both non-baryonic known components.
    # Bounds are assumed explicitly and are not a joint posterior probability.
    envelope = dict(
        dh_upper=dh + 3 * sd,
        ruler_upper=R_DRAG + 3 * R_SIGMA,
        length_lower=L_PLANCK - 3 * L_SIGMA,
        baryon_upper=OMEGA_B + 3 * BBN_ERROR,
        neutrino_upper=0.003,
        radiation_upper=0.00005,
    )
    floor = required(
        envelope["dh_upper"],
        z,
        envelope["ruler_upper"],
        envelope["length_lower"],
        envelope["baryon_upper"],
        envelope["neutrino_upper"],
        envelope["radiation_upper"],
    )
    # A high-H endpoint minimizes the forecast D_H/r_d at the coincidence ceiling.
    hubble2_max = (
        float(2 * envelope["baryon_upper"] + envelope["neutrino_upper"]) * (1 + z) ** 3
        + envelope["radiation_upper"] * (1 + z) ** 4
        + omega_vac(envelope["length_lower"])
    )
    min_dh_at_ceiling = C_LIGHT / (100 * envelope["ruler_upper"] * np.sqrt(hubble2_max))
    return dict(
        traveller=traveller_count,
        coincidence=coincidence_count,
        addresses=len(words),
        measured_ratio=str(measured),
        universal_unit_weight_ceiling=1,
        physical_baryon_density=OMEGA_B,
        nominal_neutrino=nominal_nu,
        nominal_radiation=nominal_r,
        z=z,
        DH_over_rd=dh,
        sigma_DH=sd,
        cases=cases,
        ratio_required_at_central_inputs=required(
            dh, z, R_DRAG, L_PLANCK, OMEGA_B, nominal_nu, nominal_r
        ),
        external_box=envelope,
        minimum_required_ratio=floor,
        minimum_DH_at_ceiling=min_dh_at_ceiling,
        curved_constant_vacuum=curved_background_floor(
            data,
            envelope["ruler_upper"],
            envelope["baryon_upper"],
            envelope["neutrino_upper"],
            envelope["radiation_upper"],
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    result = run()
    target = HERE / "collar_abundance.json"
    if opt.write:
        target.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    elif not verify(result, json.loads(target.read_text())):
        raise ValueError(
            "Published-source extrapolation differs from archived comparison"
        )
    print("Published-source abundance comparison replayed")


if __name__ == "__main__":
    main()
