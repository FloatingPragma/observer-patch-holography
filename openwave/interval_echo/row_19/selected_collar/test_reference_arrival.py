"""Check charge transmission, archived histories, ordinary flux and sky exclusions."""

from fractions import Fraction as F
from hashlib import sha256
from math import log, sqrt

import numpy as np
import sympy as sp

from openwave.interval_echo.row_19.selected_collar import archived_overlap as archive
from openwave.interval_echo.row_19.selected_collar import modular_stock_sky as stock

local = archive.local


def overlap(a, b):
    return max(a[0], b[0]) <= min(a[1], b[1])


def test_receipts_bind_measured_inputs_and_protocols():
    for name in (
        "reference_arrival",
        "dust_source_read",
        "archived_overlap",
        "archive_partitions",
        "stationary_work",
        "phase_and_flux",
        "modular_stock_sky",
        "modular_arrival_flight",
    ):
        for path, expected in local.inputs(name)["inputs"].items():
            assert sha256((local.HERE / path).read_bytes()).hexdigest() == expected


def test_source_information_required_by_the_pressureless_tower():
    row = local.inputs("dust_source_read")
    n = row["tickets"]
    block = np.array(row["sharp_counts"], dtype=float)
    weights = -np.log(block / block.sum())
    interaction = weights - weights[:, [0]] - weights[[0], :] + weights[0, 0]
    determinant = np.linalg.det(interaction[1:, 1:])
    assert 0 < row["observed_minor"][0] <= determinant <= row["observed_minor"][1]
    for expected, ticks, gain in zip(
        interaction[1:], row["thresholds"], row["gains"], strict=True
    ):
        scale = np.mean(gain)
        np.testing.assert_allclose(
            expected - min(expected),
            np.array(ticks) * scale / n,
            rtol=0,
            atol=scale / n + 1e-15,
        )
    arms = row["arms"]
    assert (
        arms["same_serial"]["arrival_pair_counts_per_2592"]
        == arms["fresh"]["arrival_pair_counts_per_2592"]
    )
    assert arms["same_serial"]["full_charge"][0] > 0.00156
    for arm in arms.values():
        assert overlap(arm["full_charge"], arm["charge_from_two_counts"])
    assert (
        row["dust_inferred_change_in_counted_mass"][0]
        <= F(8, 33)
        <= row["dust_inferred_change_in_counted_mass"][1]
    )
    rods = row["dust_held_rod_comparison"]
    assert [v["inverse_radius"] for v in rods] == [512, 2048, 8192]
    assert all(abs(v["predicted_rest_component"] - 19 / 22) < 1e-8 for v in rods)
    errors = [v["maximum_component_difference"] for v in rods]
    assert errors[2] < errors[1] / 15 < errors[0] / 225 and errors[2] < 8.2e-8


def test_reference_repeat_measures_every_sharp_modular_cost():
    z = local.inputs("reference_arrival")
    ref = local.survey.standards(local.survey.bulk.federation(0))[0][0].ravel() / 2592
    h = sum(F(1, i) for i in range(1, z["panel_limit"] + 1))
    assert z["harmonic_normalization"][0] <= h <= z["harmonic_normalization"][1]
    assert 0 < z["max_charge_error_nats"] < 2e-14
    for p, a, seen, records in zip(
        ref,
        z["pass_rate_by_word"],
        z["literal_passes_by_word"],
        z["record_minima_by_word"],
        strict=True,
    ):
        assert abs(-log(p) - float(h) * a) < z["max_charge_error_nats"] + 2e-15
        assert abs(seen / z["literal_starts_per_word"] - a) < sqrt(
            log(2 * 81 / 0.001) / (2 * z["literal_starts_per_word"])
        )
        assert abs(records / z["literal_starts_per_word"] + log(p)) < 0.15
    assert z["literal_reference_comparisons"] >= 81 * z["literal_starts_per_word"]


def test_positive_charges_with_zero_or_opposite_serial_transfer():
    arms = local.inputs("reference_arrival")["arms"]
    base = np.array(arms["fresh"]["source_program"])
    for k in ("shared", "same_serial", "opposed_serial"):
        a = arms[k]
        d = np.array(a["source_program"]) - base
        assert d.sum() == 0 and not d.sum(axis=0).any() and not d.sum(axis=1).any()
        assert a["measured_modular_contrast"][0] > 0
        assert overlap(a["sharp_law_modular_contrast"], a["measured_modular_contrast"])
        assert overlap(
            a["measured_arrival_change"], a["arrival_from_word_measurements"]
        )
    assert F(arms["same_serial"]["serial_acceptance"]) == F(25, 36)
    assert F(arms["opposed_serial"]["serial_acceptance"]) - F(25, 36) == -F(1, 144)


def test_two_source_records_reach_dispatch_using_blank_swaps():
    z = local.inputs("reference_bus")
    assert z["complete_symbol_pairs"] == 81 and z["paths"][0] == [0]
    assert z["paths"][1][-1] == z["collector"] and len(z["paths"][1]) == 8
    for file, digest in z["inputs"].items():
        assert sha256((local.HERE / file).read_bytes()).hexdigest() == digest
    board = local.survey.bulk.federation(0)
    seams = set(board.seams)
    assert all(
        tuple(sorted((a, b))) in seams for a, b in zip(z["paths"][1], z["paths"][1][1:])
    )
    arms = local.inputs("reference_arrival")["arms"]
    for key, row in z["forward_paths"].items():
        assert sum(row["flag_counts_all_starts"]) == z["starts_per_arm"]
        assert (
            abs(
                row["accepted"] / z["starts_per_arm"]
                - arms[key]["collision_acceptance"]
            )
            < 0.015
        )


def test_twenty_four_held_instrument_predictions():
    z = local.inputs("modular_arrival_flight")
    assert z["gates"] == [62, 229] and z["time"] == 0.125
    n = 0
    for key, row in z["outcomes"].items():
        assert sum(row["available_counts"]) == z["starts_per_setting"]
        assert all(sum(k) == z["starts_per_setting"] for k in row["held_counts"])
        if key.endswith("fresh"):
            continue
        for j, (a, b) in enumerate(
            zip(
                row["source_only_prediction"],
                row["observed_position_and_flight_change"],
                strict=True,
            )
        ):
            assert overlap(a, b)
            sign = (1 if key.startswith("forward") else -1) * (-1 if j == 3 else 1)
            assert min(sign * x for x in a + b) > 0
            n += 1
    assert n == 24


def test_quiet_null_event_and_entropy_obstruction():
    z = local.inputs("archived_overlap")
    assert z["quiet_global_rank"] == 4
    row = z["surveys"][0]
    null = z["null_reads"][row["null_measurements"][0]]
    assert F(null["ordinary_hit_given_quiet"]) == F(7, 24)
    assert sum(F(x) ** 2 for x in null["exact_null_vector"]) == F(7, 24)
    assert F(row["ordinary_log_coefficient_lower"]) == F(null["ordinary_joint_hit"]) / 2
    for r in z["null_reads"].values():
        for dose, k in r["hit_miss_other_counts"].items():
            assert sum(k) == z["starts"]
            assert (
                abs(k[0] / z["starts"] - float(F(dose) * F(r["ordinary_joint_hit"])))
                < 2e-10
            )
        assert r["hit_miss_other_counts"]["0"][0] == 0
    assert sum("null_measurements" in r for r in z["surveys"]) == 10
    assert all(
        r["finite"]["overlap_gain"][1] < 0 < r["finite"]["ordinary_gain"][0]
        for r in z["surveys"]
        if "finite" in r
    )


def test_every_locally_readable_partition_in_both_radius_choices():
    z = local.inputs("archive_partitions")
    counts = []
    for r in z["surveys"]:
        counts.append(len(r["gains"]))
        assert r["largest_upper"] == max(v[1] for v in r["gains"].values()) < -7.28
    assert sorted(counts) == [30] * 4 + [225] * 8 and sum(counts) == 1920
    for r in z["history_conditional_tomography"].values():
        assert sum(r["positions"]) == z["starts_per_setting"]
        assert len({(p["path"][0], p["path"][-1]) for p in r["pairs"]}) == 24 * 23 // 2
        assert all(
            sum(v) == z["starts_per_setting"]
            for p in r["pairs"]
            for v in p["quadratures"]
        )


def test_phase_calibration_needs_a_contracting_area():
    z = local.inputs("phase_and_flux")
    assert 0.193 < z["joined_entropy_derivative_lower"] < 0.194
    assert z["required_area_derivative_upper_at_material_calibration"] < -0.023
    assert z["phase_marker_volume_and_area_change"] == [0, 0]
    work = z["phase_energy_change"]
    a, b = z["phase_release_flux_change"]
    assert overlap(work, [-2 * b, -2 * a]) and a < b < 0 < work[0]
    assert z["phase_finite_flight_flux_change"][1] < 0
    assert z["flux_common_coefficient_determinant"][1] < 0
    # Sharp native bridge read swaps the two phase outcomes and preserves absence.
    for e in z["cuts"][0]:
        k = z["arms"]["fresh"]["seam_counts"][str(tuple(e))]
        q = z["arms"]["phase"]["seam_counts"][str(tuple(e))]
        assert abs(k[0] - q[1]) / z["starts"] < 3e-10


def test_mass_area_and_entropy_after_free_work_on_both_charge_readers():
    z = local.inputs("stationary_work")
    n = z["starts"]
    assert F(z["dose"]) == F(1, 4) and z["finite_wait_trace_error_upper"] < 3e-11
    for arm in z["source"].values():
        seams = arm["work_seam_counts"]
        assert len(seams["fresh"]) == 5
        work = (
            sum(seams["phase"][key][1] - k[1] for key, k in seams["fresh"].items()) / n
        )
        assert arm["kinetic_work"][0] < work < arm["kinetic_work"][1]
        assert all(sum(k) == n for records in seams.values() for k in records.values())
    registered = {
        (r["center"], r["direction"], tuple(r["hop_radii"])): r
        for r in local.inputs("restored_balls")["surveys"]
    }
    seen = set()
    rate = sum(z["clock_ratio"]) / 2
    for r in z["surveys"]:
        key = r["center"], r["direction"], tuple(r["radii"])
        seen.add((r["instrument"], key))
        zones = registered[key]
        assert len(r["fresh_counts"]) == 8 and sum(r["fresh_counts"]) == n
        raw_area = [
            zones["cuts"][j] * rate ** (2 * slow)
            for j in (0, 1)
            for slow in (0, 1)
            for _ in (0, 1)
        ]
        for name, v in r["changes"].items():
            counts = v["held_counts"]
            assert (
                len(counts) == 8
                and sum(counts) == n
                and 0 < v["public_probability"] < 1
            )
            delta = [
                (x - y) / n for x, y in zip(counts, r["fresh_counts"], strict=True)
            ]
            measured = sum(d * a for d, a in zip(delta, raw_area, strict=True))
            assert (
                v["held_area_difference"][0] <= measured <= v["held_area_difference"][1]
            )
            mass = sum(delta[1::2])
            assert v["held_mass_difference"][0] <= mass <= v["held_mass_difference"][1]
            assert v["held_volume_difference"][0] < 0 < v["held_volume_difference"][1]
            assert v["entropy_tangent"][0] > 0
            expected_area_sign = (
                -1
                if name == "work" or (name == "ordinary" and r["center"] == 62)
                else 1
            )
            for field in ("held_equal_volume_area_change", "area_tangent"):
                assert min(expected_area_sign * a for a in v[field]) > 0
            if name == "work":
                assert v["held_equal_volume_mass_change"][0] > 0
            if name == "ordinary":
                more_inside = r["center"] == 62 and r["radii"] != [2, 3]
                assert (
                    min(
                        (1 if more_inside else -1) * a
                        for a in v["held_equal_volume_mass_change"]
                    )
                    > 0
                )
        for field in ("work_ordinary_determinant", "overlap_ordinary_determinant"):
            assert r[field][0] * r[field][1] > 0
        missing = r["additional_entropy_at_localized_ordinary_calibration"]
        if r["center"] == 62:
            assert missing["overlap"][1] < 0
        else:
            assert missing is None
    assert seen == {
        (instrument, key)
        for instrument in ("serial", "reference")
        for key in registered
    }
    assert len(seen) == len(z["surveys"]) == 24


def test_equal_enclosed_mass_and_positive_boundary_response():
    z = local.inputs("equal_mass_survey")
    assert z["far_site"] == 229 and len(z["records"]) == 4
    for file, digest in z["inputs"].items():
        assert sha256((local.HERE / file).read_bytes()).hexdigest() == digest
    for row in z["records"].values():
        assert (
            row["mass_neutral_area_tangent"][0] > 0
            and row["mass_neutral_entropy_tangent"][0] > 0
        )
        assert row["held_volume_change"][0] < 0 < row["held_volume_change"][1]
        assert row["held_mass_change"][0] < 0 < row["held_mass_change"][1]
        a, b = row["predicted_matched_quarter_area_change"], row["held_area_change"]
        assert min(a + b) > 0 and overlap(a, b)
        assert (
            len(row["held_branch_counts"]) == 16
            and sum(row["held_branch_counts"])
            == sum(row["fresh_counts"])
            == z["starts"]
        )
        calibrated = row["matched_finite_localization_fraction"]
        actual = sum(row["held_branch_counts"][8:]) / z["starts"]
        assert calibrated[0] < actual < calibrated[1]


def test_stock_bound_even_with_independent_ring_and_inventory_nuisances():
    z = local.inputs("modular_stock_sky")
    c = z["source_caps"]
    assert z["rings"] == 1288 and len(z["galaxies"]) == 136
    for name in ("min_vertex", "max_vertex"):
        m = sp.Matrix(3, 3, list(map(F, c[name])))
        p = sp.Matrix(list(map(F, c["class_masses"])))
        assert m * sp.ones(3, 1) == m.T * sp.ones(3, 1) == p
    assert (
        c["fixed_marginal_span_upper"]
        < 0.008878
        < c["any_population_span_upper"]
        < 0.246
    )
    assert list(z["summary"]["expanded_HI"].values()) == [85, 61, 17]
    assert list(z["summary"]["catalogue"].values()) == [105, 83, 37]
    # Rebuild the strongest expanded-HI ring without running the catalogue classifier.
    g = z["galaxies"]["NGC3741"]
    r = g["bounds"]["expanded_HI"]
    meta, rot = [
        (stock.sky.sample.measure.DATA / n).read_text().splitlines()
        for n in ("table1.dat", "table2.dat")
    ]
    t = meta[g["table1_line"] - 1]
    radius, slow, _, (gas, disk, bulge), _ = stock.sky.rings_in(
        rot, [{"line": r["witness_line"]}]
    )[0]
    assert t[:11].strip() == "NGC3741"
    mass = 10**9 * (
        F(7, 5) * (F(t[40:47]) + 2 * F(t[48:55])) + F(266, 100) * F(t[86:93])
    )
    _, hmax = stock.sky.sample.slope_box(t)
    maxim = max(
        gas * a + disk * b + bulge * d
        for a in (F(1, 2), 2)
        for b, d in stock.sky.sample.RECT
    )
    req = (
        radius
        * max(0, slow / F(stock.sky.sample.upper(hmax)) - maxim)
        / (stock.G_HIGH * mass)
    )
    assert abs(float(req - F(r["required_fraction_lower"]))) < 1e-12 and req > 4
