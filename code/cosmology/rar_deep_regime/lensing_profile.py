"""Test the projected deep profile against published KiDS shear with a frozen SPARC scale."""

import argparse
import hashlib
import io
import json
import tarfile
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear

HERE = Path(__file__).resolve().parent
CALIBRATION = HERE / "runtime/matched_observable_receipt.json"
RECEIPT = HERE / "receipts/lensing_profile.json"
SPARC_DEEP = HERE / "receipts/sparc_deep_regime_diagnostic.json"
ARCHIVE = "https://kids.strw.leidenuniv.nl/sci_data/brouwer2021_rar.tar"
ARCHIVE_SHA256 = "73adc2d4e848fa0a6a43187f7b5447e749b92eaec4a795b9902b0beaf8a78733"
PROFILES = {
    "KiDS-bright": "Fig-4-5-C1_RAR-KiDS-isolated",
    "GAMA": "Fig-4-C1_RAR-GAMA-isolated",
    "KiDS-bright, hot-gas mass sensitivity": "Fig-4_RAR-KiDS-isolated_hotgas",
}
# Brouwer et al.'s h70=1 conversion of one Msun/pc² of excess surface density to 4G DeltaSigma.
FOUR_G = 4 * 4.52e-30 * 3.086e16


def load_profile(archive, prefix):
    data = {}
    pins = {}
    for kind in ("Nobins", "covmatrix"):
        name = f"{prefix}_{kind}.txt"
        with archive.extractfile(name) as member:
            payload = member.read()
        pins[name] = hashlib.sha256(payload).hexdigest()
        data[kind] = np.loadtxt(io.BytesIO(payload))
    bins = data["Nobins"]
    raw = data["covmatrix"]
    centres, i = np.unique(raw[:, 2], return_inverse=True)
    other, j = np.unique(raw[:, 3], return_inverse=True)
    n = len(bins)
    np.testing.assert_array_equal(centres, other)
    np.testing.assert_allclose(centres, bins[:, 0], rtol=5e-5, atol=0)
    np.testing.assert_array_equal(np.sort(n * i + j), np.arange(n * n))
    cov = np.empty((n, n))
    cov[i, j] = raw[:, 4] / raw[:, 6]
    np.testing.assert_allclose(cov, cov.T, rtol=5e-10, atol=1e-11)
    np.testing.assert_allclose(
        np.sqrt(cov.diagonal()), bins[:, 3] / bins[:, 4], rtol=2e-4
    )
    # Fig. 4 uses midpoints of 15 logarithmically spaced acceleration intervals.
    edges = np.geomspace(1e-15, 5e-12, 16)
    np.testing.assert_allclose(
        centres, (edges[:-1] + edges[1:]) / 2, rtol=2e-10, atol=0
    )
    keep = (centres >= 1e-13) & (centres < 0.1 * 1.2e-10)
    return dict(
        members=pins,
        supplied_bins=n,
        g_bar_m_s2=centres[keep].tolist(),
        lower_bin_edge_m_s2=edges[:-1][keep].tolist(),
        tangential_esd=(bins[keep, 1] / bins[keep, 4]).tolist(),
        cross_esd=(bins[keep, 2] / bins[keep, 4]).tolist(),
        covariance=cov[np.ix_(keep, keep)].tolist(),
    )


def mass_and_bin_controls(b, lower, density, cross, cov, a0, sparc_interval):
    """Fixed-membership sensitivity on bins wholly above the isolation floor."""
    upper = 2 * b - lower
    whitener = np.linalg.solve(np.linalg.cholesky(cov), np.eye(len(b)))

    def projected(mean_b, mean_sqrt_b, shift, acceleration_scale):
        m = 10.0**shift
        return (
            4 / np.pi * m * mean_b + np.sqrt(m * acceleration_scale) * mean_sqrt_b
        ) / FOUR_G

    def chi(pred):
        white = whitener @ (density - pred)
        return float(white @ white)

    # Uniform pairs per sky area: R dR is proportional to db / b² for each lens.
    normalization = 1 / lower - 1 / upper
    area_b = np.log(upper / lower) / normalization
    area_sqrt_b = 2 * (1 / np.sqrt(lower) - 1 / np.sqrt(upper)) / normalization
    centers = (b, np.sqrt(b))
    area = (area_b, area_sqrt_b)
    shifts = np.linspace(-0.2, 0.2, 81)
    samples = {"centers": [], "uniform_sky_area": [], "independent_bin_means": []}
    edge_mixtures = []
    for shift in shifts:
        for name, moments in (("centers", centers), ("uniform_sky_area", area)):
            samples[name].append(chi(projected(*moments, shift, a0)))
        lo = projected(lower, np.sqrt(lower), shift, a0)
        hi = projected(upper, np.sqrt(upper), shift, a0)
        bounded = lsq_linear(
            whitener, whitener @ density, bounds=(lo, hi), method="bvls"
        )
        if not bounded.success:
            raise RuntimeError(bounded.message)
        samples["independent_bin_means"].append(chi(bounded.x))
        edge_mixtures.append((bounded.x - lo) / (hi - lo))

    results = {}
    for name, likelihood in samples.items():
        best = int(np.argmin(likelihood))
        results[name] = {
            "nominal_mass_chisq": likelihood[40],
            "lowest_sampled_chisq": likelihood[best],
            "total_mass_shift_dex_at_lowest": round(float(shifts[best]), 3),
        }
        if name == "independent_bin_means":
            results[name]["upper_edge_mixture_at_lowest"] = np.clip(
                edge_mixtures[best], 0, 1
            ).tolist()

    # For a common distribution of b/bin_center the normalized anomaly is
    # constant, and positive compact baryons give a term increasing in sqrt(b).
    # Compare equal-count halves, using only bin order to define the contrast.
    if len(b) != 6:
        raise ValueError("The isolation-floor contrast requires six bins.")
    contrast = np.repeat([1.0, -1.0], 3) / (3 * np.sqrt(b / 1e-12))
    sigma = float(np.sqrt(contrast @ cov @ contrast))
    results["low_minus_high_normalized_esd"] = {
        "tangential": float(contrast @ density),
        "assigned_statistical_se": sigma,
        "standardized_excess_over_zero_upper_bound": float(contrast @ density) / sigma,
        "cross_over_same_statistical_se": float(contrast @ cross) / sigma,
    }
    results["center_chisq_nominal_mass_at_sparc_interval_endpoints"] = [
        chi(projected(*centers, 0.0, scale)) for scale in sparc_interval
    ]
    return results


def score(profile, a0, sparc_interval):
    b = np.array(profile["g_bar_m_s2"])
    density = np.array(profile["tangential_esd"])
    cov = np.array(profile["covariance"])
    chol = np.linalg.cholesky(cov)
    # DeltaSigma_b = M_b/(pi R²). Abel projection of rho_A proportional to 1/r²
    # gives DeltaSigma_A = sqrt(a0 g_bar)/(4G). Both terms assume spherical weak-field lensing.
    compact = 4 / np.pi * b / FOUR_G
    envelope = np.sqrt(a0 * b) / FOUR_G

    def quadratic(residual):
        white = np.linalg.solve(chol, residual)
        return float(white @ white)

    isolated = np.array(profile["lower_bin_edge_m_s2"]) >= 1e-13
    residual = (density - compact - envelope)[isolated]
    return dict(
        bins=len(b),
        wholly_above_isolation_floor=dict(
            bins=int(isolated.sum()),
            chisq_profile=float(
                residual @ np.linalg.solve(cov[np.ix_(isolated, isolated)], residual)
            ),
            controls=mass_and_bin_controls(
                b[isolated],
                np.array(profile["lower_bin_edge_m_s2"])[isolated],
                density[isolated],
                np.array(profile["cross_esd"])[isolated],
                cov[np.ix_(isolated, isolated)],
                a0,
                sparc_interval,
            ),
        ),
        projected_baryons=compact.tolist(),
        projected_anomaly=envelope.tolist(),
        chisq_profile=quadratic(density - compact - envelope),
        chisq_baryons_only=quadratic(density - compact),
        chisq_cross_zero=quadratic(profile["cross_esd"]),
    )


def build(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != ARCHIVE_SHA256:
        raise ValueError("KiDS public archive digest changed")
    calibration = json.loads(CALIBRATION.read_text())
    a0 = calibration["primary_fraction_0p1"]["channel_a_rar"]["a0_m_s2"]
    rotation = json.loads(SPARC_DEEP.read_text())["deep_radial_acceleration"]
    np.testing.assert_allclose(
        rotation["a0_fixed_exponent_m_s2"], a0, rtol=1e-14, atol=0
    )
    interval = rotation["a0_galaxy_cluster_bootstrap_95pct_m_s2"]
    with tarfile.open(path) as archive:
        profiles = {
            name: load_profile(archive, prefix) for name, prefix in PROFILES.items()
        }
    return dict(
        public_data=ARCHIVE,
        archive_sha256=ARCHIVE_SHA256,
        calibration_source=str(CALIBRATION.relative_to(HERE)),
        fitted_sparc_a0_m_s2=a0,
        sparc_cluster_bootstrap_95pct_m_s2=interval,
        mass_and_bin_protocol={
            "systematic": "Brouwer et al. Fig. 4 quotes ±0.2 dex on stellar masses. Scan a common total compact-baryon rescaling -0.200 to +0.200 dex by 0.005 dex, holding observed memberships and covariances fixed. Gas fractions and galaxy-by-galaxy rebinning are unavailable here.",
            "bins": "Use the six wholly above the isolation floor. Uniform sky area assumes uniform background pair density and weights. Independent bin means allow any nonnegative within-bin weights, including different edge mixtures in each bin.",
            "scores": "Lowest scores on the mass grid and bin-edge mixtures use the observed shear. These are feasibility/sensitivity diagnostics, with no fitted degrees of freedom assigned.",
            "shape": "Subtract the mean of the upper three from the lower three ESDs after dividing by sqrt(gbar_center / 1e-12 m/s²). A common fractional pair distribution in every bin and positive compact and isothermal amplitudes predict a nonpositive difference.",
        },
        four_G_m_s2_per_Msun_pc2=FOUR_G,
        minimum_inclusive_gbar_m_s2=1e-13,
        maximum_exclusive_gbar_m_s2=1.2e-11,
        profiles={
            key: dict(**profile, **score(profile, a0, interval))
            for key, profile in profiles.items()
        },
    )


def check(saved, actual):
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            raise ValueError("Mismatched lensing receipt keys")
        for key, value in actual.items():
            check(saved[key], value)
    elif isinstance(actual, (float, list)):
        np.testing.assert_allclose(saved, actual, atol=0, rtol=2e-10)
    elif actual != saved:
        raise ValueError("Mismatched lensing data provenance")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "archive", type=Path, help="Brouwer et al. public Figure 4 tarball"
    )
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = build(args.archive)
    if args.write:
        RECEIPT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    else:
        check(json.loads(RECEIPT.read_text()), result)
    for key, data in result["profiles"].items():
        print(
            key,
            {
                name: data[name]
                for name in (
                    "bins",
                    "chisq_profile",
                    "chisq_baryons_only",
                    "chisq_cross_zero",
                )
            },
        )
