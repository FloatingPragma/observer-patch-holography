import json
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm

from cosmology.rar_deep_regime.compensation import measure as m


def test_uniform_shell_conserves_energy_and_saturates_radius_bound():
    R, M, rho = 0.291, 3.4e12, 2.1e10
    T = m.donor_radius(R, M, rho)
    assert 4 * np.pi / 3 * rho * (T**3 - R**3) == pytest.approx(M)
    # Construct the post-transfer physical density over the entire donor:
    # the core gains M, the shell is exactly empty, the exterior is unchanged.
    core = 3 * M / (4 * np.pi * R**3)
    assert (rho + core) * R**3 == pytest.approx(rho * T**3)
    assert m.compensated_shear_mass(1.1 * T, R, T, M) == 0
    assert m.compensated_shear_mass(T, R, T, M) == 0
    assert m.compensated_shear_mass(R, R, T, M) > M


def test_actual_deprojection_of_finite_signed_source_recovers_enclosed_excess():
    R, T, M = 0.5, 3.0, 2.0
    for r in (0.25, 0.5, 1.5, 2.99, 3.0, 4.0):

        def integrand(theta):
            projected_r = r / np.sin(theta)
            return m.compensated_shear_mass(projected_r, R, T, M) / (
                np.pi * projected_r**2
            )

        knots = [np.arcsin(r / a) for a in (R, T) if r < a]
        mass = 4 * r**2 * quad(integrand, 0, np.pi / 2, points=knots, epsabs=1e-12)[0]
        expected = M * (min(r, R) / R) ** 3 - M * (
            max(0, min(r, T) ** 3 - R**3) / (T**3 - R**3)
        )
        assert mass == pytest.approx(expected, abs=2e-10)


def test_aperture_mass_and_error_propagation_against_si():
    inputs = json.loads((m.HERE / "inputs.json").read_text())
    out = m.run(inputs)
    radius_kpc = out["apertures"]["outer"]["radius_kpc"]
    assert radius_kpc == 291.07
    rows = np.array(inputs["mistele_all"]["outer"]["v_stat_syst_kms"])
    Rsi = radius_kpc * 3.085677581491367e19
    v, s, d = rows.T
    # Reproduce using physical constants and the mass-squared-speed Jacobian.
    bound = (
        Rsi * 1e6 * (v * v - 2 * v * (2 * s + d)) / (6.67430e-11 * m.MSUN_KG)
        - 2 * 10**m.LOG_BARYON_UPPER
    )
    stored = out["apertures"]["outer"]["cases"]["two_systematic_double_baryons"]
    np.testing.assert_allclose(
        bound, stored["min_net_excess_and_exterior_debit_msun"], rtol=6e-5
    )
    assert 2.2414027276 == pytest.approx(norm.ppf(1 - 0.05 / 4), abs=1e-10)
    assert np.all(
        np.array(stored["min_exterior_debit_over_allowed_galaxy_baryons"]) > 5
    )


def test_stacked_bound_allows_a_minority_of_large_donors():
    # Signed profiles can compensate fully inside the aperture for some lenses.
    R, rho = 0.3, 2.0
    radii = np.array([0.05, 0.10, 1.20])
    weights = np.array([0.4, 0.3, 0.3])
    deficits = 4 * np.pi / 3 * rho * np.maximum(radii**3 - R**3, 0)
    lower = m.donor_radius(R, weights @ deficits, rho)
    assert lower == pytest.approx((weights @ np.maximum(radii, R) ** 3) ** (1 / 3))
    assert min(radii) < lower < max(radii)
    assert (weights @ radii**3) ** (1 / 3) < lower


def test_shear_bound_and_additive_mass_sheet():
    x = json.loads((m.HERE / "inputs.json").read_text())
    shear = m.direct_shear(x)
    cap = shear["caps"]["catalog_mean_source_efficiency"]["compact_ceiling_h_msun_pc2"]
    sigma = 100.0
    mass = shear["max_ordinary_budget_msun"]
    r = shear["largest_pair_radius_lower_limit_kpc"] * 1000
    assert (mass / (np.pi * r * r) + sigma - sigma) / 0.7 == pytest.approx(cap)
    assert cap < 60.0 - 2 * 17
    assert (
        shear["caps"]["arbitrary_source_efficiencies_up_to_one"][
            "compact_ceiling_h_msun_pc2"
        ]
        > 60.0 - 2 * 17
    )
