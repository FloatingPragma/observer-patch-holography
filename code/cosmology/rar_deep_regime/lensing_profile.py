"""Test the projected deep profile against published KiDS shear with a frozen SPARC scale."""

import argparse
import hashlib
import io
import json
import tarfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CALIBRATION = HERE / "runtime/matched_observable_receipt.json"
RECEIPT = HERE / "receipts/lensing_profile.json"
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


def score(profile, a0):
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
    with tarfile.open(path) as archive:
        profiles = {
            name: load_profile(archive, prefix) for name, prefix in PROFILES.items()
        }
    return dict(
        public_data=ARCHIVE,
        archive_sha256=ARCHIVE_SHA256,
        calibration_source=str(CALIBRATION.relative_to(HERE)),
        fitted_sparc_a0_m_s2=a0,
        four_G_m_s2_per_Msun_pc2=FOUR_G,
        minimum_inclusive_gbar_m_s2=1e-13,
        maximum_exclusive_gbar_m_s2=1.2e-11,
        profiles={
            key: dict(**profile, **score(profile, a0))
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
