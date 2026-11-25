import json
import sys
from pathlib import Path
from math import pi, sqrt

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cosmology.rar_deep_regime.double_ring_pressure import measure as m


def test_shell_extrema_and_direct_positive_exterior():
    # Direct solid-angle integration, independently of the integration by parts
    # in the certificate. For radius < R2 the strongest extra mass is a core.
    r1, r2 = 1.43, 2.07

    def admitted(R, r):
        return 1 - sqrt(max(0.0, 1 - (R / r) ** 2))

    def weight(r):
        return (admitted(r1, r) / r1**2 - admitted(r2, r) / r2**2) / pi

    core = (r1**-2 - r2**-2) / pi
    for r in np.geomspace(0.01, 100.0, 101):
        assert weight(r) <= core + 1e-14
        if r >= r2:
            assert weight(r) < 0

    # Smooth finite M(r)=r^3/(1+r^2)^(3/2), increasing on all radii.
    def M(r):
        return r**3 / (1 + r * r) ** 1.5

    def density(r):
        return 3 * r**2 / (1 + r * r) ** 2.5

    s = sqrt(1 - (r1 / r2) ** 2)
    upper = (
        s
        / pi
        * (
            quad(lambda x: M(r1 / sqrt(1 - (s * x) ** 2)), 0, 1)[0] / r1**2
            - M(r2) / (r2**2 * (1 + s))
        )
    )
    inner = quad(lambda r: density(r) * weight(r), 0, r2, points=[r1])[0]
    outer = quad(lambda r: density(r) * weight(r), r2, np.inf)[0]
    assert abs(inner - upper) < 1e-10 and outer < 0


def test_two_plane_sign_and_primary_sheet_invariance():
    # Source-plane deflection enters the second lens equation once and affects
    # neither the observed velocity nor the first lensing critical density.
    scale, crit, eff = m.geometry()
    unrestricted = m.nominal(m.LOG_MASSES["Salpeter"], m.acceleration()[0])
    S1, S2 = unrestricted["unlimited_orbital_means"]
    sig = unrestricted["zero_pressure_secondary_equivalent_km_s"]
    b2 = 4 * pi * (sig / m.C_KM_S) ** 2 * eff * m.ARCSEC_PER_RADIAN
    sheet = unrestricted["zero_pressure_primary_sheet_kappa1"] * crit[0]
    assert abs((S1 + sheet) / crit[0] - 1) < 1e-13
    assert abs((S2 + sheet) / crit[1] + b2 / m.THETA[1] - 1) < 1e-13
    assert m.THETA[1] * (1 - (S2 + sheet) / crit[0]) > 0
    primary = unrestricted["stellar_fibre"]
    spectrum = next(
        row
        for row in json.loads(
            (m.HERE.parent / "aperture_virial/sample.json").read_text()
        )
        if row["name"] == "J0946+1006"
    )
    low, high = primary["anisotropy_bounds_km_s"]
    assert low < primary["isotropic_km_s"] < high
    assert abs(primary["isotropic_km_s"] - spectrum["sigma"]) < spectrum["e_sigma"]
    for foreground in (-200, 0, 100):
        c1 = crit[0] - foreground
        c2 = crit[1] * (1 - b2 / m.THETA[1]) - foreground
        assert abs((S1 - S2) - (c1 - c2)) < 1e-8
    assert scale > 0


def test_positive_contrast_bound_and_source_inputs():
    document = json.loads((m.HERE / "results.json").read_text())
    assert json.loads(json.dumps(m.certified())) == document["all_positive_exteriors"]
    for panel in document["all_positive_exteriors"].values():
        assert panel["pressure_contrast_upper"] < 0
        assert panel["additional_interior_mass_lower_msun"] > 0
    assert document["nominal"]["Chabrier"]["zero_pressure_secondary_bend_arcsec"] < 0
    np.testing.assert_allclose(
        document["radii_kpc"], np.array(m.THETA) * m.geometry()[0]
    )
