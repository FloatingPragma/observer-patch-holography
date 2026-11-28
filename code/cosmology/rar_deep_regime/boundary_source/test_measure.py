"""Independent coarea, distance, aperture and out-of-sample checks."""

import json

import numpy as np
from scipy.integrate import quad

from .distance import coarse_tv_envelope, run as distance_run
from .measure import run, scores
from .source import AREA_UNITS, HERE, THRESHOLD_MAX, rebin, shell_factor, variation


def test_shared_threshold_integral_matches_direct_perimeters():
    edges = np.array([0, 1, 2.5, 4.0])
    sigma = np.array([12.0, 5.0, 9.0])
    levels = sorted(set([0.0, *sigma, THRESHOLD_MAX]))
    exact = 0.0
    for lo, hi in zip(levels, levels[1:]):
        inside = np.r_[sigma > (lo + hi) / 2, False]
        perimeter = (
            sum(
                r
                for r, change in zip(edges[1:], inside[:-1] != inside[1:], strict=True)
                if change
            )
            * 2
            * np.pi
        )
        exact += (hi - lo) * perimeter
    np.testing.assert_allclose(sum(variation(edges, sigma)), exact * AREA_UNITS)
    top_hat = np.array([6.0, 6.0, 6.0])
    np.testing.assert_allclose(
        sum(variation(edges, top_hat)), 6 * 8 * np.pi * AREA_UNITS
    )


def test_exponential_at_fine_resolution():
    rd = 1.25
    e = np.linspace(0, 20 * rd, 5001)
    sigma = 100 * np.exp(-(e[1:] + e[:-1]) / (2 * rd))
    mass = np.pi * AREA_UNITS * np.diff(e**2) @ sigma
    np.testing.assert_allclose(sum(variation(e, sigma)), mass / rd, rtol=5e-6)


def test_transport_integral_on_both_sides_of_source():
    shells, weight = np.array([0.25, 1.2, 4.1]), np.array([12.0, 1.1, 9.0])
    for r in (0.12, 2.3, 17):

        def integrand(s):
            return (
                sum(w * s / max(s, b) for b, w in zip(shells, weight, strict=True))
                / sum(weight)
            ) ** 2

        integral = quad(integrand, 0, r, points=shells[shells < r])[0] / r
        np.testing.assert_allclose(
            shell_factor(r, shells, weight), integral, atol=2e-13
        )


def test_rebin_mass_and_all_distance_envelope():
    rng = np.random.default_rng(3198)
    for _ in range(20):
        edges = np.r_[0, np.cumsum(rng.uniform(0.02, 1, 12))]
        sigma = rng.uniform(0, 220, 12)
        lo, hi = coarse_tv_envelope(edges[1:], sigma, 0.8, 0.6, 1.6)
        for scale in np.linspace(0.6, 1.6, 61):
            new, values = rebin(edges * scale, sigma, 0.8)
            np.testing.assert_allclose(
                np.diff(new**2) @ values,
                np.diff((edges * scale) ** 2) @ sigma,
                rtol=1e-12,
            )
            assert lo - 1e-5 <= sum(variation(new, values)) <= hi + 1e-5


def test_held_labels_cannot_change_gain_or_prediction():
    receipt = run()
    targets = {k: v.copy() for k, v in receipt["targets"].items()}
    for n in receipt["held"]:
        targets[n]["excess_v2"] *= 99
        targets[n]["observed_kms"] += 500
    other = scores(receipt["profiles"], targets, receipt["train"], receipt["held"])
    for law, result in receipt["held_predictions"].items():
        assert other[law]["gain"] == result["gain"]
        for n, p in result["held"].items():
            assert (
                p["predicted_total_kms"] == other[law]["held"][n]["predicted_total_kms"]
            )


def test_receipts_and_measurement_units():
    receipt = run()
    assert json.loads(json.dumps(receipt)) == json.loads(
        (HERE / "results.json").read_text()
    )
    assert len(receipt["train"]) == len(receipt["held"]) == 5
    assert set(receipt["train"]).isdisjoint(receipt["held"])
    assert all(
        v["max_sigma_msun_pc2"] <= THRESHOLD_MAX for v in receipt["profiles"].values()
    )
    assert set(receipt["pair_high_low"]) == {"NGC3198", "NGC3521"}
    dist = distance_run()
    assert dist == json.loads((HERE / "distance.json").read_text())
    for row in dist.values():
        assert row["required_gain_ratio_envelope"][1] < 1
        for n, v in row["galaxies"].items():
            assert v["S_envelope"][0] <= v["S"] <= v["S_envelope"][1]
            assert v["Y_box"][0] < v["Y"] < v["Y_box"][1]
            assert v["edge_kpc"] < v["orbit_kpc"]


def test_square_lattice_boundary_and_finite_shell_response():
    from .maps import mapped_shell_response, square_tv

    # Centred 2x2 square, four sides, with every boundary bond at the same radius.
    source = np.zeros((8, 8))
    source[3:5, 3:5] = 7
    width = 0.8
    expected = 8 * width * 7 * AREA_UNITS
    np.testing.assert_allclose(square_tv(source, width), expected)
    r = np.hypot(width, width / 2)
    for radius in (0.3, 2.2):
        enclosed, weighted = mapped_shell_response(source, width, 0, radius)
        np.testing.assert_allclose(enclosed, 0 if radius < r else expected)
        expected_factor = (
            radius**2 / (3 * r**2) if radius < r else 1 - 2 * r / (3 * radius)
        )
        np.testing.assert_allclose(weighted, expected * expected_factor)


def test_softened_disk_force_by_direct_orbital_quadrature():
    from .measure import G
    from .stellar_maps import ordinary_disk_v2

    radius = 8.0
    masses = np.array([3.0e8, 2.0e9, 1e9, 9e7])
    rs = np.array([0.0, 4.0, 8.0, 10.0])
    for height in (0.2, 0.4, 0.8):

        def direct(angle):
            return (
                G
                * radius
                * np.sum(
                    masses
                    * (radius - rs * np.cos(angle))
                    / (radius**2 + rs**2 + height**2 - 2 * radius * rs * np.cos(angle))
                    ** 1.5
                )
            )

        mean = quad(direct, 0, np.pi, epsabs=1e-8)[0] / np.pi
        np.testing.assert_allclose(
            ordinary_disk_v2(radius, masses, rs, height), mean, rtol=2e-12
        )


def test_all_grid_receipts_and_source_flux():
    from .maps import run as atom_report
    from .stellar_maps import run as phase_report

    gas = atom_report()
    phase = phase_report()
    assert gas == json.loads((HERE / "maps.json").read_text())
    assert phase == json.loads((HERE / "stellar_co.json").read_text())
    source = json.loads((HERE / "map_inputs.json").read_text())
    for n, entry in gas["profiles"].items():
        meta = source[n]
        np.testing.assert_allclose(
            entry["atomic_He_msun"],
            meta["positive_flux_jy_kms"] * 2.356e5 * meta["distance_mpc"] ** 2 * 1.36,
        )
        assert (
            entry["atomic_mass_fraction_beyond_orbit"]
            < entry["atomic_mass_fraction_beyond_optical_aperture"]
        )
    r = phase["ratios"]["0.8kpc_0offset"]
    np.testing.assert_allclose(r["mapped_all"], 4.08655, rtol=2e-5)
    assert r["mapped_both_gases"] < r["mapped_all"] < r["mapped_stars"]


def test_gas_transfer_replays_without_held_gains():
    from .gas_transfer import CALIBRATOR, OTHER_PAIR_MEMBER, run as radial_run
    from .gas_image_transfer import run as image_run

    sources = json.loads((HERE / "gas_image_inputs.json").read_text())
    target = run()["targets"]
    image = image_run(sources)
    radial = radial_run()
    assert image == json.loads((HERE / "gas_image_transfer.json").read_text())
    assert radial == json.loads((HERE / "gas_transfer.json").read_text())
    for law, item in image["predictions"].items():
        ids = set(item["held"])
        assert len(ids) == 8 and ids.isdisjoint([CALIBRATOR, OTHER_PAIR_MEMBER])
        for n, entry in item["held"].items():
            pred_v2 = (
                target[CALIBRATOR]["excess_v2"]
                * sources[n]["source_msun_kpc"][law]
                / sources[CALIBRATOR]["source_msun_kpc"][law]
            )
            np.testing.assert_allclose(entry["predicted_excess_v2"], pred_v2)
    nominal = image["predictions"]["0.8_both_gases_compact"]["held"]
    assert (
        sum(
            not v["common_gain_possible_with_quoted_speed_errors"]["2"]
            for v in nominal.values()
        )
        == 4
    )
