"""Check the physical readouts, finite sky certificates and continuum normalizations."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from cosmology.rar_deep_regime.brightness_pairs import midpoint_halo as sky
from openwave.interval_echo.row_19.selected_collar import midpoint_bulk as bulk

HERE = Path(__file__).resolve().parent


def receipt(name):
    z = json.loads((HERE / f"{name}.json").read_bytes())
    assert all(
        sha256((HERE / f).read_bytes()).hexdigest() == d for f, d in z["inputs"].items()
    )
    return z


def test_spherical_source_and_finite_reads():
    z = receipt("midpoint_bulk")
    for r in (0.25, 0.75, 15 / 16, 1.25, 1.5):
        # Independent integration of both positive angular caps at fixed
        # shell radius, including shells wholly inside larger apertures.
        def shell(t):
            if not t:
                return 0
            k = min(1, max(-1, (1 + t * t - r * r) / (2 * t)))
            return t * (
                (1 - k) * (1 - (1 + t * t) / (r * r)) + t * (1 - k * k) / (r * r)
            )

        weighted = quad(shell, 0, 1, points=[abs(1 - r)], epsabs=1e-13)[0]
        endpoint, _, expected = bulk.continuum(r, False)
        assert abs(weighted - expected) < 1e-11
        if r <= 1:
            assert abs(endpoint - 15 * r * weighted / 32) < 1e-12
        else:
            assert endpoint > 15 * r * weighted / 32
    for screen in z["levels"]:
        for row in screen["observations"]:
            mesh = screen["angular_mesh_upper"]
            r = row["radius"]
            for j, name in enumerate(("endpoint", "midpoint", "weighted_midpoint")):
                panel = row[name]
                p = panel["finite_admission"]
                starts = row["distinct_starts"]
                assert abs(p - row["independent_address_hits"][j] / starts) < 0.009
                assert all(
                    sum(t) == z["starts_per_word_arm"] and min(t) >= 0
                    for t in panel["counts"]
                )
                forecast = p * np.mean(z["conditional_charge"])
                assert panel["charge"][0] <= forecast <= panel["charge"][1]
            if row["center"] is None:
                lo, hi = row["separation_read"]["enclosed_fraction"]
                assert lo <= row["midpoint"]["finite_admission"] <= hi
            elif row["address_law"] == "solid_angle":
                p = row["endpoint"]["finite_admission"]
                selfpair = screen["angular_self_pair"]
                assert (
                    max(0, r - mesh) ** 4 / 16 - selfpair
                    <= p
                    <= min(2, r + mesh) ** 4 / 16 + selfpair
                )
                assert (
                    abs(
                        row["weighted_midpoint"]["finite_admission"]
                        - row["round_limit_admissions"][2]
                    )
                    <= 2 * mesh / r + selfpair
                )
    rows = z["levels"][-1]["observations"]
    for center in (0, 128, 255):
        need = [
            p["additional_central_mass_fraction"]
            for p in rows
            if p["address_law"] == "solid_angle"
            and p["center"] == center
            and p["radius"] > 1
        ]
        assert need[0][1] < need[1][0]


def test_compact_ordinary_and_exterior_read():
    z = receipt("settled_midpoint_force")
    o = z["ordinary_mode"]
    assert o["compact_laplacian_spectrum"] == [[4, 1], [6, 2]]
    probes = o["observer_reads"]
    for t, rows in probes["times"].items():
        assert rows[0][0] == z["shots"] and sum(rows[0]) == z["shots"]
    assert probes["times"]["64"][1][0] < 0.7 * z["shots"]
    assert all(e["energy"][0] <= 42 <= e["energy"][1] for e in probes["energy"])
    assert min(a[0] for a in z["potential_edges"][0]["pair_increments"]) > 0
    assert (
        0.286
        < z["potential_difference_from_records"][0]
        < z["potential_difference_from_records"][1]
        < 0.289
    )
    for q, offset in (("1", 8.5), ("2", 16.25)):
        m = z["detector_controls"][q]
        assert m["settled"]["held"][0] <= 0 <= m["settled"]["held"][1]
        assert m["neighbor"]["held"][0] <= offset <= m["neighbor"]["held"][1]
        assert (
            z["overlap"][q]["held_conditional"][0]
            <= 0
            <= z["overlap"][q]["held_conditional"][1]
        )


def verify(rows, cap, cert):
    if cert["possible"]:
        x = [F(p) / b for p, b in zip(cert["parameters"], cap, strict=True)]
        assert all(0 <= z <= 1 for z in x)
        assert all(sum(v * w for v, w in zip(a, x, strict=True)) <= b for a, b in rows)
    else:
        left = [F(0)] * len(cap)
        right = 0
        for index, word in cert["multipliers"]:
            y = F(word)
            assert y > 0
            a, b = rows[index]
            left = [u + y * v for u, v in zip(left, a, strict=True)]
            right += y * b
        assert sum(min(z, 0) for z in left) - right == F(cert["exact_gap"]) > 0


def test_midpoint_sky_witnesses(monkeypatch):
    data = json.loads((sky.HERE / "midpoint_halo.json").read_bytes())
    raw = (sky.sample.measure.DATA / "table2.dat").read_text().splitlines()
    meta = {
        v[:11].strip(): v
        for v in (sky.sample.measure.DATA / "table1.dat").read_text().splitlines()
    }
    monkeypatch.setattr(sky.sample.iv, "dps", 36)
    assert set(data["summary"]["arbitrary_centered_mixture_excluded"]) == {
        "NGC5907",
        "UGC00128",
        "UGC09133",
    }
    for name, row in data["galaxies"].items():
        entries = [{"line": i} for i in row["lines"]]
        lo, hi = sky.sample.slope_box(meta[name])

        def check(entries, cut, cert, extra=F(0)):
            rings = sky.rings_in(raw, entries, extra)
            lower, upper = (
                (sky.sample.upper, sky.sample.lower)
                if cert["possible"]
                else (sky.sample.lower, sky.sample.upper)
            )
            rows, cap = sky.problem(rings, max(0, F(lower(lo))), F(upper(hi)), cut)
            verify(rows, cap, cert)

        check(entries, None, row["nonnegative_mixture"])
        for i, cert in enumerate(row["by_number_of_inside_radii"]):
            check(entries, i, cert)
        if not row["nonnegative_mixture"]["possible"]:
            w = row["mixture_witness"]
            assert len(w["lines"]) == 3
            check([{"line": i} for i in w["lines"]], None, w["proof"])
            extra = row["additional_speed_halfwidth_kms"]
            check(entries, None, extra["low_certificate"], F(extra["excluded_at"]))
            check(entries, None, extra["high_certificate"], F(extra["feasible_at"]))


def test_lensing_velocity_certificates():
    result = json.loads((sky.HERE / "midpoint_halo.json").read_bytes())
    for path in (
        sky.HERE / "outer_requirement.json",
        sky.sample.measure.DATA / "table1.dat",
        sky.sample.measure.DATA / "table2.dat",
        sky.PROTOCOL,
        sky.LENS_INPUT,
        sky.COLD_INPUT,
        sky.WORD_INPUT,
    ):
        assert sha256(path.read_bytes()).hexdigest() == result["inputs"][path.name]
    cosmology = result["cosmological_dilution"]
    require_lo, require_hi = map(
        F, cosmology["required_power_of_total_pair_price_in_comoving_volume"]
    )
    assert require_lo == F("0.992860") and require_hi == F("1.002688")
    factor = cosmology["required_pair_price_factor_per_scale_doubling"]
    np.testing.assert_allclose(factor, [2 ** float(require_lo), 2 ** float(require_hi)])
    assert 1 < factor[0] < 2 < factor[1]
    src = cosmology["measured_pair_price_nats"]
    mass = cosmology["mass_times_radius_per_pair_in_natural_units"]
    np.testing.assert_allclose(mass, 15 * np.array(src) / (32 * np.pi))
    data = json.loads(sky.LENS_INPUT.read_bytes())
    assert (
        data["sha256"]
        == "2d2f446b70fbc74d8fa6fe66b51f69ac7615dcec1214a4f75bb62f900c947c71"
    )
    for col, mass in enumerate(("10.5", "10.8", "11.1", "11.5")):
        for factor in (1, 2):
            panel = result["deprojected_lens_stacks"][mass][str(factor)]
            cap = F(panel["compact_GM_cap"])
            # Verify the upper-mass rounding using rational tenths in log mass.
            rounded_mass = cap / F("0.000004300917270")
            assert rounded_mass.denominator == 1
            exact_tenth_power = 2**10 * 10 ** int(F(mass) * 10)
            assert (rounded_mass - 1) ** 10 < exact_tenth_power <= rounded_mass**10
            rings = []
            for i, row in enumerate(data["rows"]):
                r = F(str(row["radius_kpc"]))
                v, scatter, syst = map(
                    lambda a: F(str(a)), row["velocity_stat_sys"][col]
                )
                dv = factor * (scatter + syst)
                rings.append((r, max(0, v - dv) ** 2, (v + dv) ** 2, [1 / r], i))
            for cut, cert in [*enumerate(panel["cuts"]), (None, panel["mixture"])]:
                rows, cap_vector = sky.problem(rings, F(1), F(1), cut, cap)
                verify(rows, cap_vector, cert)
            assert panel["mixture"]["possible"]
            assert panel["any_edge"] == any(z["possible"] for z in panel["cuts"])
    assert [
        sum(
            result["deprojected_lens_stacks"][k][str(f)]["any_edge"]
            for k in result["deprojected_lens_stacks"]
        )
        for f in (1, 2)
    ] == [2, 3]


def test_direct_shear_projection_and_pairs():
    from cosmology.rar_deep_regime.brightness_pairs import midpoint_shear as lens

    # Integrate the column and projected cylinder independently in 3D,
    # then check the two sharp pair bounds throughout the allowed edge range.
    for radius in (0.05, 0.35, 0.93, 0.999, 1.0, 1.1, 2.0):
        edge = 1.0
        zmax = np.sqrt(max(0, edge**2 - radius**2))
        column = 2 * quad(lambda z: 1 / np.hypot(radius, z), 0, zmax)[0]

        def shell(r):
            fraction = 1 if r <= radius else 1 - np.sqrt(1 - (radius / r) ** 2)
            return 4 * np.pi * r * fraction

        cylinder = quad(shell, 0, edge, points=[min(radius, edge)])[0]
        shear = cylinder / (np.pi * radius**2) - column
        assert abs(shear - 2 * lens.esd(radius, edge)) < 1e-9
    for near in (0.2, 0.6, 0.93):
        low, high = lens.bounds(near**2)
        for edge in np.r_[np.geomspace(0.01, 100, 300), 1]:
            ratio = lens.esd(1, edge) / lens.esd(near, edge)
            assert low - 1e-12 <= ratio <= high + 1e-12
        assert abs(lens.esd(1, 1) / lens.esd(near, 1) - high) < 1e-13
    stored = json.loads((lens.HERE / "midpoint_shear.json").read_bytes())
    lens.lensing_profile.check(stored, lens.run())
    for survey, expected in (("KiDS-bright", -1.90072), ("GAMA", 0.16118)):
        panel = stored["surveys"][survey]["tangential_esd"]
        assert len(panel["pairs"]) == 15
        assert abs(panel["maxima"]["centre"]["z"] - expected) < 0.00001


def test_native_rest_clock_records():
    raw = receipt("midpoint_clocks")
    assert raw["exterior_sharp_pair_configurations"] == 486
    assert all(
        m["multiplicities_at_4_and_6"] == [1, 2]
        for m in raw["fixed_modes_per_cell"].values()
    )
    assert (
        sum(
            4 * len(data["counts"])
            for arms in raw["material_clock_reads"].values()
            for data in arms.values()
        )
        + 4 * len(raw["times"])
        == raw["number_of_settings"]
    )
    for owner, arms in raw["material_clock_reads"].items():
        mode = raw["fixed_modes_per_cell"][owner]
        a, b = np.asarray(mode["levels"], float)
        port = mode["port"] % 12
        u, v = a[port] / np.linalg.norm(a), b[port] / np.linalg.norm(b)
        for panel in arms.values():
            for t, phase, tapes in zip(
                raw["times"], panel["phases"], panel["counts"], strict=True
            ):
                assert phase["frequency"][0] <= 1 <= phase["frequency"][1]
                for j, (hit, miss) in enumerate(tapes):
                    chance = (u * u + v * v) / 2 + u * v * np.cos(t - j * np.pi / 2)
                    assert hit + miss == raw["shots_per_setting"]
                    assert (
                        abs(hit / (hit + miss) - chance)
                        <= raw["simultaneous_probability_error"]
                    )
    assert all(
        lo <= 0 <= hi and hi - lo < 0.00045
        for lo, hi in raw["change_of_clock_difference"].values()
    )
    assert raw["near_local_word"]["largest_absolute_occupation_change"] > 0.48


def test_three_comparable_regions():
    data = receipt("three_regions")
    h, b, c = (
        np.sqrt(float(F(data["radius_squared"][k])))
        for k in ("half", "boundary", "central")
    )
    assert abs(2 * h * h - b * b - c * c + 0.5) < 1e-14
    rng = np.random.default_rng(2026)
    for size in (0.1, 1, 13):
        locations = size * rng.normal(size=(1000, 3))

        def norm2(points):
            return np.sum(points**2, axis=1)

        def weight(distance, radius):
            return np.maximum(radius**2 - distance, 0)

        axis = np.array([1, 0, 0])
        kernel = (
            2 * weight(norm2(locations - axis / 2), h)
            - weight(norm2(locations - axis), b)
            - weight(norm2(locations), c)
        )
        assert np.max(kernel) <= 5e-15
    for row in data["levels"]:
        assert row["endpoint_residual"][0] > 0 > row["weighted_midpoint_residual"][1]
        assert row["central"]["endpoint_hits"] == 0
        for name in ("boundary", "central", "half"):
            region = row[name]
            for instrument in ("endpoint", "weighted_midpoint"):
                expected = region[instrument]["finite_admission"]
                assert (
                    abs(
                        region[instrument + "_hits"] / region["distinct_address_starts"]
                        - expected
                    )
                    < 0.01
                )
        small = row["half"]["endpoint"]["charge"][0]
        upper = (
            b * row["boundary"]["endpoint"]["charge"][1]
            + c * row["central"]["endpoint"]["charge"][1]
        )
        np.testing.assert_allclose(
            2 * h * small - upper, row["endpoint_residual"][0], rtol=1e-12
        )


def test_equal_volume_fillings_and_area_forecast():
    import sympy as sp

    r = sp.symbols("r", real=True)
    S = sp.sqrt(15) / 4
    U = r**5 * (1 - r**2) ** 4 * (r**2 - S**2)
    h = sp.diff(U, r)
    assert U.subs(r, 0) == U.subs(r, S) == U.subs(r, 1) == 0
    assert h.subs(r, 1) == 0
    assert sp.simplify(h.subs(r, S) - 2 * S**6 * (1 - S**2) ** 4) == 0
    ratio = sp.Poly(sp.cancel(h / r**2), r)
    assert all(k[0] % 2 == 0 for k in ratio.monoms())
    z = receipt("three_regions")["levels"][-1]
    for name, x in (
        ("half", 0.75),
        ("central", np.sqrt(15) / 4),
        ("boundary", np.sqrt(11) / 4),
    ):
        expected = -15 * x / 32 * np.array(z[name]["weighted_midpoint"]["charge"][::-1])
        np.testing.assert_allclose(
            z[name]["stationary_gain_times_area_forecast"], expected
        )
    assert z["missing_regional_term_lower_bound"] > 1.2029e-5
    assert z["negative_mass_lower_bound"] > 1.4430e-5
    raw = receipt("three_regions")
    for row in raw["signed_stress_sweep"]:
        halfwidth = float(F(row["halfwidth"]))
        assert max(raw["radial_normalization_errors"]) < halfwidth
        assert (
            row["below"]["finite_admission"] == 0
            and row["above"]["finite_admission"] == 1
        )
        low = row["above"]["charge"][0] - row["below"]["charge"][1]
        assert low > 0.0033049
        np.testing.assert_allclose(
            row["minimum_total_variation"], low / (4 * np.pi * halfwidth)
        )
    assert raw["signed_stress_sweep"][-1]["minimum_in_midpoint_masses"] > 2184.45
