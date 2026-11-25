"""Source custody, absolute gas inventory and independent compact-halo cross-checks."""

import hashlib
from math import sqrt

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

from cosmology.rar_deep_regime.warped_outer_rings import acquire, measure as m


def test_published_columns_and_enrollment():
    path = m.HERE / "rings.json"
    assert (
        hashlib.sha256(path.read_bytes()).hexdigest()
        == "7b1d4c54a3aa1ccb2de4497aeabd467120183d350f04eeb76ed70cf9d1436bed"
    )
    r = m.read(path)
    assert r["source"] == acquire.SOURCE and len(r["rows"]) == 21
    profile = r["rows"]
    assert [t["speed_kms"] for t in profile[-4:]] == [143.9, 139.6, 146.7, 168.8]
    for k in (0, 1, 2):
        x = np.array([p["radius_arcsec"] for p in profile])
        y = np.array([max(0, p["NHI_1e19_cm2"] - k * p["error_NHI"]) for p in profile])
        for target in x[-4:]:
            upper = 0.9 * target
            z = quad(
                lambda u: u * np.interp(u, x, y),
                0,
                upper,
                points=x[(x > 0) & (x < upper)],
                epsabs=1e-7,
            )[0]
            g = z * 2 * np.pi * (16.1e6 * m.ARCSEC) ** 2 * m.SURFACE * 1.4
            np.testing.assert_allclose(
                m.gas_inside(profile, target, 16.1, k), g, rtol=5e-15
            )
    result = m.run()
    for rr in result["sample"]:
        b = rr["enclosure"]["1"]
        assert 0.903 < b["fraction_retaining_95pct"] < 0.904
        assert 0.90 < b["minimum_retention_for_90pct"] < 0.91
        assert rr["enclosure"]["2"]["fraction_retaining_95pct"] < 0.89
        assert rr["gbar_m_s2"] < 1.2e-11
    assert result["prior_enclosure_lower"] < 0.83


def test_independent_edge_phase_distance_and_ordering():
    saved = m.read(m.HERE / "results.json")
    rebuilt = m.run()

    # Integer effects are absent; the finite population model uses elementary
    # floating arithmetic, so replay its numerical fields with a rounding allowance.
    def same(a, b):
        if isinstance(a, dict):
            assert a.keys() == b.keys()
            for k in a:
                same(a[k], b[k])
        elif isinstance(a, list):
            assert len(a) == len(b)
            for x, y in zip(a, b, strict=True):
                same(x, y)
        elif isinstance(a, float):
            np.testing.assert_allclose(a, b, rtol=2e-13, atol=1e-14)
        else:
            assert a == b

    same(saved, rebuilt)
    r = saved["sample"]
    gm = m.G_KPC * saved["mass_msun"]
    a0 = saved["reference_a0_m_s2"]
    distances = np.array([x["radius_kpc"] for x in r])
    v = np.array([x["speed_kms"] for x in r])
    e = np.array([x["error_speed"] for x in r])
    plateau = sqrt(gm * a0 * m.KPC / 1e6)

    def worst(L):
        pred = np.sqrt(gm / distances + plateau * np.minimum(1, L / distances))
        return np.max(abs(pred - v) / e)

    least = saved["sharp_edge_nominal"]["2"]["lower_kpc"]
    assert 29.25 < least < 29.27 and abs(worst(least) - 2) < 2e-13
    assert worst(least - 0.001) > 2 and worst(least + 0.001) < 2
    assert not saved["sharp_edge_nominal"]["1"]["feasible"] and worst(1e4) > 1.5
    pred = np.array([x["forecast_kms"] for x in r])
    fit = minimize_scalar(
        lambda z: np.max(abs(z * pred - v) / e),
        bounds=(0.5, 1.5),
        method="bounded",
        options={"xatol": 1e-12},
    )
    assert (
        abs(fit.fun - saved["common_distance_control"]["minimax_error_multiplier"])
        < 3e-7
    )
    clock = saved["finite_spacing_calibration"]
    sp = clock["spacing_kpc"]
    u = distances / sp
    scores = []
    for phi in np.arange(2000) / 2000:
        pr = np.sqrt(
            gm / distances
            + sqrt(gm * clock["a0_m_s2"] * m.KPC / 1e6) * np.floor(u + phi) / u
        )
        scores.append(max(abs(pr - v) / e))
    np.testing.assert_allclose(
        saved["shared_cut_phase"]["minimax_error_multiplier_range"],
        [min(scores), max(scores)],
    )
    d = saved["cross_galaxy_at_outer_ring"]
    assert d["nominal"]["required_minimum_excess_ratio"] > 7
    assert d["one_error"]["against"]["tilt_and_lens_mass"]["minimum_gap"] > 0.70
    assert d["two_errors"]["against"]["tilt_and_lens_mass"]["minimum_gap"] < 0
    assert all(
        c["radial_order_preserved"]
        for row in d.values()
        for c in row["against"].values()
    )
