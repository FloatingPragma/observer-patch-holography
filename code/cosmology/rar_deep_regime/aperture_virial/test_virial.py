"""Independent aperture, lens-pressure and catalogue controls."""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, quad, simpson
from scipy.special import exp1

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cosmology.rar_deep_regime.aperture_virial import measure


def test_top_hat_for_anisotropic_stellar_orbits():
    # nu=exp(-r), v_c²=1 and Osipkov-Merritt anisotropy with r_a=0.6.
    for aperture in (0.35, 1.3, 4.0):

        def pressures(r):
            density = np.exp(-r)
            radial = ((r + 1) * density + 0.6**2 * exp1(r)) / (r * r + 0.6**2)
            return radial, radial * 0.6**2 / (r * r + 0.6**2)

        def projected(r):
            lower = np.sqrt(max(0.0, 1 - (aperture / r) ** 2))
            a = (1 - lower**3) / 3
            b = 1 - lower - a
            pr, pt = pressures(r)
            return (
                r
                * r
                * np.array(
                    [
                        pr * a + pt * b,
                        np.exp(-r) * b / 2,
                        (r + 1) * np.exp(-r) * a / r**2,
                    ]
                )
            )

        num, low, high = [
            quad(lambda r: projected(r)[i], 0, aperture)[0]
            + quad(lambda r: projected(r)[i], aperture, np.inf)[0]
            for i in range(3)
        ]
        assert low < num < high
        surface = (
            aperture**2
            / 2
            * quad(lambda z: pressures(np.hypot(aperture, z))[0], 0, np.inf)[0]
        )
        np.testing.assert_allclose(num - low, surface, rtol=2e-9)


def test_pressure_column_and_inversion_for_signed_finite_pressure():
    # At c=G=1 take p=epsilon*(a-r²)*exp(-r²) and M_RC-M_W=2π r³p.
    a, eps = 2.0, 0.002

    def pressure(r):
        return eps * (a - r * r) * np.exp(-r * r)

    def dmass(r):
        z = r * r
        return 2 * np.pi * eps * z * np.exp(-z) * (3 * (a - z) - 2 * z * (1 + a - z))

    def column(R):
        # Integrate independent spherical mass-shell fractions inside a cylinder.
        core = quad(dmass, 0, R)[0]
        halo = quad(
            lambda r: dmass(r) * (R / r) ** 2 / (1 + np.sqrt(1 - (R / r) ** 2)),
            R,
            np.inf,
            epsabs=1e-13,
        )[0]
        return (core + halo) / (np.pi * R * R)

    def delta(R):
        return R**2 * np.sqrt(np.pi) * eps * np.exp(-(R**2)) * (a + 0.5 - R**2)

    for R in (0.45, 1.0, 1.7):
        expected = 2 * quad(lambda z: pressure(np.hypot(R, z)), 0, np.inf)[0]
        np.testing.assert_allclose(column(R), expected, atol=2e-12)
        inverted = (
            2 / (np.pi * R) * quad(lambda t: delta(R / np.cos(t)), 0, np.pi / 2)[0]
        )
        np.testing.assert_allclose(inverted, pressure(R), atol=2e-12)
    annulus = 2 * quad(lambda R: delta(R) / R, 0.45, 1.7)[0]
    np.testing.assert_allclose(annulus, column(0.45) - column(1.7), atol=2e-12)


def test_isotropic_aperture_reproduces_the_separate_slacs_reader():
    tracer = measure.Tracer()
    isotropic = cumulative_simpson(tracer.nu[::-1], x=-tracer.logr[::-1], initial=0)[
        ::-1
    ]
    archive = json.loads(
        (measure.HERE.parent / "receipts/stellar_lensing.json").read_text()
    )
    rows = json.loads((measure.HERE / "sample.json").read_text())
    for row, entry in zip(rows, archive["objects"], strict=True):
        ap = tracer.aperture(row["Re"])
        weight = tracer.radius**3 * sum(ap)
        k = simpson(weight * isotropic, x=tracer.logr) / simpson(
            weight * tracer.nu, x=tracer.logr
        )
        low, high = tracer.moment(ap, 1.0)
        assert low < k < high
        efficiency = 1 - measure.slacs.comoving(row["zFG"]) / measure.slacs.comoving(
            row["zBG"]
        )
        prediction = (
            2
            * np.pi
            * measure.slacs.ARCSEC_PER_RADIAN
            * efficiency
            * (row["sigma"] / measure.slacs.C_KM_S) ** 2
            / k
        )
        np.testing.assert_allclose(prediction, entry[1], atol=1.2e-5, rtol=0)


def test_refined_survey_comparison_and_quadratic_mass_inverse():
    saved = json.loads((measure.HERE / "results.json").read_text())
    inputs = json.loads((measure.HERE / "sample.json").read_text())
    high = measure.calc(inputs, size=1801, order=160)
    assert high["n_spectra"] == 53 and high["n_matching_stellar_masses"] == 18
    for group in high["extended_square_root"]["bSIE"]:
        for image in ("bSIE", "bLTM"):
            wide = high["extended_square_root"][image][group]
            base = saved["extended_square_root"][image][group]
            for edge in ("maximum_lower", "minimum_upper"):
                assert base[edge][0] == wide[edge][0]
                np.testing.assert_allclose(
                    base[edge][1], wide[edge][1], atol=6e-5, rtol=0
                )
    for mass in (0.025, 1.7, 26):
        np.testing.assert_allclose(
            measure.positive_mass(3 * mass + 4 * np.sqrt(mass), 3, 4), mass
        )
