import json

from . import late_read, measure, placements


def test_same_local_stress_and_geometry_have_distinct_force():
    r = measure.run()
    assert r["highest_common_position_diagonal_jet"] == 8
    assert r["ports"] == [0, 229]
    for q, rec in r["charges"].items():
        assert rec["maximum_direct_optical_table_and_jet_difference"] < 1e-12
        assert min(rec["reference_minimum_image_singular_values"]) > 1e-5
        for name in ("ordinary", "pairing"):
            a, b = json.loads(rec["after_original_flight"][name])
            assert (b < 0) if q == "-1" else (a > 0)
            observed = rec["count_slopes"][name]
            assert abs(observed["estimate"]) > observed["halfwidth"]


def test_both_native_mark_addresses_must_translate():
    r = placements.run()
    assert r["translated_addresses"] == (228, 229, 218, 219)
    for result in r["charges"].values():
        for row in result["squared_hilbert_schmidt_residual_by_parity"]:
            assert row["both_moved"] == 0 < min(row["left_moved"], row["right_moved"])
        for row in result["independent_head_residual_moments"]:
            assert all(x == [0] * 7 + [32] for x in row.values())


def test_calibrated_late_local_tomography_recovers_all_marked_force():
    rr = late_read.run()
    assert (
        rr["pair_count"] == 8 and rr["calibration_product_settings_per_channel"] == 64
    )
    for q, r in rr["charges"].items():
        assert len(r["parities"]["0"]) == 8
        assert r["parities"]["0"]["0_229"]["integer_determinant"] == 2**20
        assert (
            r["reconstructed_full_slopes"]["pairing"] < -5e-5
            if q == "-1"
            else r["reconstructed_full_slopes"]["pairing"] > 6e-5
        )


def test_ordinary_calibrated_sharp_focusing_and_held_source():
    from . import sharp_focus

    data = sharp_focus.run()
    parse = sharp_focus.parse_interval
    for q, panels in data["charges"].items():
        for direction in ("forward", "reverse", "odd", "even"):
            held = parse(
                panels["sharp_focusing"][direction]["observed_over_force_prediction"]
            )
            assert 1 not in held
        clock_test = data["clock_calibrations"]["odd"]["held_constraints"][q]
        assert 0 not in parse(clock_test["ordinary_calibrated_compatibility_residual"])
    for fit in data["two_meter_ordinary_calibrations"].values():
        assert 0 not in parse(fit["ordinary_minor"])
        assert all(
            0 not in parse(v["observed_minus_predicted_ratio"]) for v in fit["held"]
        )
    for direction in ("forward", "reverse", "odd", "even"):
        clock = data["clock_calibrations"][direction]["held_constraints"]
        a, b = (
            parse(clock[q]["required_clock_acceleration_with_free_species_gain"])
            for q in ("-1", "1")
        )
        assert a.b < b.a or b.b < a.a


def test_direct_signed_flights_and_source_differences_match_held_certificate():
    import numpy as np

    from . import sharp_focus

    board = measure.translation.Atlas(1).board
    route = measure.bridge.native.iso.pairing.formation.region(board)
    rings = measure.bridge.beam.screens(board, route)
    states = measure.bridge.tangents()
    optical = sharp_focus.screen.measure
    basis = np.c_[np.ones(5), measure.bridge.beam.BEARINGS]
    certified = sharp_focus.archived("directional_stress/connection/sharp_focus.json")
    measured = sharp_focus.archived("directional_stress/connection/late_read.json")
    step = 2**-22
    for q in (-1, 1):
        reads = {}
        for label, direction in (("forward", rings), ("reverse", rings[::-1])):
            signal = optical.signal_jets(board, route, direction, q, 4)
            slopes = {}
            for name in ("ordinary", "pairing"):
                moments = []
                for sign in (-1, 1):
                    source = states["reference"] + sign * step * states[name]
                    assert all(np.linalg.eigvalsh(block).min() > 0 for block in source)
                    probs = np.array(
                        [measure.bridge.beam.probabilities(s, source) for s in signal]
                    )
                    jac = optical.matrix_jets(basis.T @ probs @ basis)
                    moments.append(-np.trace(np.linalg.solve(jac[0].T, jac[2].T)))
                slopes[name] = (moments[1] - moments[0]) / (2 * step)
            reads[label] = slopes
        for label, sign in (("odd", -1), ("even", 1)):
            reads[label] = {
                k: v + sign * reads["reverse"][k] for k, v in reads["forward"].items()
            }
        m = measured["charges"][str(q)]["reconstructed_full_slopes"]
        for name, row in reads.items():
            held = (row["pairing"] / row["ordinary"]) / (m["pairing"] / m["ordinary"])
            bounds = certified["charges"][str(q)]["sharp_focusing"][name][
                "observed_over_force_prediction"
            ]
            assert held in sharp_focus.parse_interval(bounds)


def test_full_material_force_separates_the_exact_equal_radar_pair():
    from . import equal_radar

    result = equal_radar.run()
    assert result["range_rows_identical"] == 12
    for q, lower in ((-1, 1.2543e-8), (1, 4.8911e-7)):
        cell = result["charges"][str(q)]
        slope = json.loads(cell["force_derivative"])
        assert (
            slope[0] > lower
            and slope[0] <= float(cell["nominal_late_local_derivative"]) <= slope[1]
        )
        change = json.loads(cell["force_change_on_certified_chord"])
        assert change[0] > lower * 1e-7
        for branch in cell["pair_contributions"].values():
            assert len(branch) == 8


def test_entropy_changes_at_exactly_equal_radar_and_force():
    from . import force_entropy

    result = force_entropy.run()
    assert result["rank"] == 14
    assert result["exact_kernel_dimension"] == 4
    assert json.loads(result["signed_entropy_slope_nats"])[0] > 1.19e-5
    assert json.loads(result["entropy_increase_nats_lower"])[0] > 1e-12
    witness = result["independent_joint_read"]
    assert witness["product"] == "XX"
    assert (
        witness["observed_plus_difference"] + witness["simultaneous_99pct_halfwidth"]
        < 0
    )
