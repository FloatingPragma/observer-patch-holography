"""Cross-check the actual source alphabet, shell implication and held clocks."""

from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np
from scipy.sparse.csgraph import shortest_path

from openwave.interval_echo.row_19.selected_collar import paired_ranging as survey


def receipt(name):
    data = json.loads((survey.HERE / (name + ".json")).read_bytes())
    for path, digest in data["inputs"].items():
        assert sha256((survey.HERE / path).read_bytes()).hexdigest() == digest
    return data


def contains(band, value):
    assert band[0] <= value <= band[1]


def test_hemisphere_and_actual_paired_ranging():
    z = receipt("paired_ranging")
    assert len(z["classes"]) == 171
    assert sum(map(len, z["classes"])) == 210
    assert len(z["ordinary"]["checked_blank_load_pairs"]) == 9
    rng = np.random.default_rng(1939)
    for _ in range(50):
        x, y = [rng.uniform() * a / np.linalg.norm(a) for a in rng.normal(size=(2, 3))]
        normal = np.cross(x, y)
        normal /= np.linalg.norm(normal)
        u, v = [normal * np.sqrt(1 - a @ a) for a in (x, y)]
        lifted = x @ y + np.sqrt((1 - x @ x) * (1 - y @ y))
        assert (
            np.max(np.abs(np.linalg.norm([x + u, x - u, y + v, y - v], axis=1) - 1))
            < 1e-12
        )
        assert abs(min((x + u) @ (y + v), (x - u) @ (y - v)) - lifted) < 1e-12
    for check in z["triangle_tests"].values():
        assert check["fixed"]["violated_triangles"] == 0
        h = check["best"]["witness"]
        assert sum(h["legs"]) < h["direct"] <= sum(h["legs"]) + h["switch"]
    slopes = z["shells"]["fixed"]["area_volumes"]
    for tag, sign in (("native", 1), ("reversed_second", -1)):
        band = slopes[tag]["count_enclosure"]["area_slope"]
        contains(band, slopes[tag]["compensated"]["fixed_volume_area_slope"])
        assert min(sign * v for v in band) > 0
    for level in z["angular_surveys"]:
        assert level["single_speed_minimax_relative_error"] > 0.35
    angular = receipt("hemisphere_source")["measurements"][-3:]
    for r in angular:
        assert abs(r["isotropic_enclosed"] - 2 * r["isotropic_endpoint"]) < 1e-14
    for k in ("mass_over_price_flat", "mass_over_price_cosine"):
        assert (
            angular[0]["rows"][k][0]
            > angular[1]["rows"][k][1]
            > angular[2]["rows"][k][1]
        )


def test_energy_and_local_clock_recover_reference_entropy():
    z = receipt("pair_clock_information")
    board = survey.bulk.federation(0)
    w, turned = survey.standards(board)
    e = np.asarray(
        survey.routes.volume.radar_source_and_hand(board)[-1]["source_word_energies"]
    )
    summed = (e[:, None] + e).ravel()
    reader = z["local_reader"]
    first = reader["calibrations"]["all"]
    assert len(first["ports"]) == 1
    i, j = first["ports"][0]
    # Independently remove the physical matching, then count the internal swaps.
    for k, matching in enumerate(survey.routes.codebook()[0][1:]):
        adj = survey.routes.cell().adjacency().toarray().copy()
        for a, b in matching:
            adj[a, b] = adj[b, a] = 0
        dist = int(shortest_path(adj, unweighted=True)[i, j])
        assert [summed[k], dist] == first["sharp_records_in_calibration_order"][k]
    assert (
        len(set(summed))
        < 9
        == len(set(map(tuple, first["sharp_records_in_calibration_order"])))
    )
    primary = reader["calibrations"][reader["primary"]]
    assert primary["ports"] == [[3, 6]] and primary["minimal_decks"] == 1
    for k, matching in enumerate(survey.routes.codebook()[0][1:]):
        result = int((3, 6) not in matching)
        assert primary["sharp_records_in_calibration_order"][k] == [summed[k], result]
        np.testing.assert_allclose(
            primary["one_and_two_load_arrival_probabilities"][k], result, atol=3e-13
        )
    assert len(set(map(tuple, primary["sharp_records_in_calibration_order"]))) == 9
    for tag, laws in (("native", w[:3]), ("reversed_second", turned)):
        ref, fresh, shared = np.array(laws) / 2592
        # The two signed arms have exactly the same two-station energy law.
        _, labels = np.unique(
            [[a, b] for a in summed for b in summed], axis=0, return_inverse=True
        )
        assert np.all(np.bincount(labels, weights=(laws[2] - laws[1]).ravel()) == 0)
        q = -np.sum((shared - fresh) * np.log(ref))
        h = np.sum(ref * np.log(ref) - fresh * np.log(fresh))
        obs = reader["acquisitions"][tag]
        assert all(
            sum(v) == reader["shots_per_arm"]
            for v in obs["reference_fresh_shared_counts"]
        )
        contains(obs["modular_slope"], q)
        contains(obs["finite_entropy_increase"], h)
        rows = z["placements"][tag]
        assert sum(r["clock"]["sum_absolute_count_change"] == 0 for r in rows) == 162
        assert all(r["clock_and_energy"]["sum_absolute_count_change"] > 0 for r in rows)
        held = z["held"][tag]
        amount = sum(max(0, v) for v in held["change_of_each_joint_class"]) / 2592
        contains(held["joint_event_change_interval"], amount)


def test_fixed_volume_and_independent_entropy_residual():
    z = receipt("paired_balls")
    assert (z["markers"], z["boundary_links"]) == (171, 570)
    rows = z["metric_results"]["fixed"]
    changes = []
    for row in rows:
        areas = []
        for obs in row["readings"].values():
            trial = obs["finite_dose"]
            lo, hi = obs["population"]
            a = F(trial["upper_deadline_fraction"])

            def finite(k):
                return (
                    (1 - a) * sum(map(F, lo[k]))
                    + a * sum(map(F, hi[k]))
                    - (F(lo[k][0]) + F(hi[k][0])) / 2
                )

            assert all(finite(k) == F(v) for k, v in trial["exact_changes"].items())
            assert finite("V") == 0 and trial["measured_entropy_change"][0] > 0
            areas.append(finite("A"))
        assert areas[0] * areas[1] < 0
        changes.append(areas[0])
    assert sum(a > 0 for a in changes) == 6
    cert = z["calibration_obstructions"]["fixed"]["native"]["finite"]["equality"]
    weights = {i: F(q) for i, q in cert["center_weights"]}
    residual = F(0)
    columns = dict.fromkeys(("A", "screen", "pair"), F(0))
    for row in rows:
        w = weights.get(row["center"], F(0))
        trial = row["readings"]["native"]["finite_dose"]
        for k in columns:
            columns[k] += w * F(trial["exact_changes"][k])
        residual += w * F(trial["measured_entropy_change"][0 if w >= 0 else 1])
    assert sum(map(abs, weights.values())) == 1 and all(
        v == 0 for v in columns.values()
    )
    assert residual == F(cert["certified_minimum_max_error_nats"]) > F("0.00091683")
    free = z["metric_results"]["best"][0]["readings"]["native"]["finite_dose"]
    assert all(F(v) == 0 for v in free["exact_changes"].values())
    assert free["measured_entropy_change"][0] > 0.00107727


def test_enclosed_mass_area_slope_and_cosmological_bound():
    z = receipt("curved_shells")
    cs = np.mean(z["cosine_shell_locations"], axis=1)
    weights = np.array(z["shell_multiplicities"]) / 190
    assert weights.sum() == 1
    for tag, rows in z["regimes"].items():
        w = weights.copy()
        if tag == "empty_center":
            w[0] = 0
        for i, row in enumerate(rows):
            mass = sum((cs * w)[: i + 1])
            for key in ("direct_normalized_mass", "normalized_mass_read"):
                contains(row[key], mass)
            r = np.sqrt(np.mean(row["squared_areal_radius"]))
            c = np.sqrt(1 - r * r)
            area = -8 * np.pi * (mass - c * w[: i + 1].sum()) / r
            contains(row["area_slope_per_G_mass_unit"], area)
            contains(row["inward_acceleration_per_G_mass_unit"], mass / r**2)
            if tag == "complete":
                contains(row["local_detector_normalized_mass"], mass)
                for cut in row["cuts"]:
                    assert all(
                        len(a) == 82 and sum(a) == 2**51
                        for a in cut["clock_and_energy"]["counts"]
                    )
            else:
                for test in row["nonlinear_spatial_constraint"]:
                    assert all(
                        abs(t["volume_difference"]) < 1e-11
                        for t in test["finite_volume_checks"]
                    )
                assert abs(test["centered_area_slope"] - area) < 2e-7
    cosmic = z["cosmic_doubling"]
    assert (
        cosmic["largest_conserved_price_density_ratio"][1]
        < cosmic["observed_density_ratio"][0]
    )
    assert cosmic["minimum_price_growth_measured_profile"][0] > 1.38611
    assert cosmic["minimum_price_growth_isotropic"][0] > 1.32675


def test_spherical_chord_population_and_sparc_certificates():
    from scipy.integrate import quad
    from openwave.interval_echo.row_19.selected_collar import curved_orbits as orbits
    from cosmology.rar_deep_regime.brightness_pairs import midpoint_halo as sky

    z = receipt("curved_orbits")
    old = json.loads((sky.HERE / "midpoint_halo.json").read_bytes())
    assert set(old["summary"]["arbitrary_centered_mixture_excluded"]) == set(
        z["groups"]["excluded"]
    )
    assert len(z["groups"]["feasible"]) == 133 and not z["groups"]["unresolved"]
    for x in np.linspace(0.01, 0.99, 11):
        # t=1-s^2 removes the endpoint singularity in the density of flat screens.
        mixture = quad(
            lambda s: 1.5 * (1 - s * s) * min(x / (1 - s * s), 1),
            0,
            1,
            points=[np.sqrt(1 - x)],
            epsabs=1e-12,
        )[0]
        assert abs(mixture - (1 - (1 - x) ** 1.5)) < 1e-11
    assert [r["stations"] for r in z["chord_survey"]] == [20, 80, 320, 1280, 5120]
    for trial in z["chord_survey"]:
        for row in trial["radii"]:
            contains(row["refinement_envelope"], row["population_active_fraction"])
            observed = row["address_gate_hits"] / trial["distinct_address_starts"]
            assert abs(observed - row["population_active_fraction"]) < 0.005
    raw = (sky.sample.measure.DATA / "table2.dat").read_text().splitlines()
    selected = json.loads((sky.HERE / "outer_requirement.json").read_bytes())[
        "galaxies"
    ]
    meta = {
        v[:11].strip(): v
        for v in (sky.sample.measure.DATA / "table1.dat").read_text().splitlines()
    }
    assert (
        sum(v["one_screen"]["outcome"] == "feasible" for v in z["galaxies"].values())
        == 120
    )
    assert sum(len(v["catalogue_lines"]) for v in z["galaxies"].values()) == 1288
    for name, trial in z["galaxies"].items():
        rings = sky.rings_in(raw, selected[name]["rings"])
        bounds = sky.sample.slope_box(meta[name])
        for choice in (trial, trial["one_screen"]):
            if choice["outcome"] == "feasible":
                h, d, b = map(F, choice["baryonic_parameters"])
                assert (
                    F(sky.sample.lower(bounds[0]))
                    <= h
                    <= F(sky.sample.upper(bounds[1]))
                )
                assert h / 10 <= d <= h and h / 10 <= b <= 7 * h / 5
                screens = [
                    (F(t), F(w))
                    for t, w in choice[
                        "edge_coordinate_and_squared_speed_at_first_ring"
                    ]
                ]
                assert all(0 < t <= 1 and w >= 0 for t, w in screens)
                for r, lo, hi, parts, _ in rings:
                    k = r / rings[0][0]
                    bary = sum(a * v for a, v in zip(parts, [h, d, b], strict=True))
                    bands = [(w, orbits.shape(k, t, t)) for t, w in screens]
                    assert lo <= bary + sum(w * lr[0] for w, lr in bands)
                    assert bary + sum(w * lr[1] for w, lr in bands) <= hi
        single = trial["one_screen"]
        if single["outcome"] == "excluded":
            specs, cap, dark = orbits.geometry(rings, bounds, False)
            partition = sorted(single["cover"], key=lambda row: F(row["ends"][0]))
            last = F(0)
            for box in partition:
                a, b = map(F, box["ends"])
                assert a == last and b > a
                last = b
                rhs = F(0)
                left = [F(0)] * 4
                for j, w in box["multipliers"]:
                    v, r, k, sign = specs[j]
                    coef = (
                        (dark * sign * orbits.shape(k, a, b)[0 if sign > 0 else 1])
                        if sign
                        else 0
                    )
                    col = [*(c * x for c, x in zip(cap, v, strict=True)), coef]
                    scale = max(F(1), abs(r), *map(abs, col))
                    weight = F(w) / scale
                    assert weight >= 0
                    rhs += weight * r
                    left = [x + weight * y for x, y in zip(left, col, strict=True)]
                assert sum(min(0, v) for v in left) - rhs == F(box["gap"]) > 0
            assert last == 1
