"""Independent checks of the angular kernel, source perturbation and affine-clock test."""

import json

import numpy as np
import pytest

from . import measure
from ..action_response.four_ticks import parse_interval


def receipt(name):
    return json.loads((measure.HERE / name).read_text())


def test_common_tidal_matrix_with_nonaffine_clock():
    tidal = np.array([[1.2, 0.35], [0.35, 0.6]])
    rate, acceleration = 1.7, -0.8
    data = []
    for image, affine_velocity in [
        (np.array([[2.0, -0.3], [0.5, 1.0]]), np.array([[0.7, 0.1], [-0.4, 0.2]])),
        (np.array([[1.0, 0.2], [0.6, 1.2]]), np.array([[0.1, -0.5], [0.3, 0.4]])),
    ]:
        jets = [
            image,
            rate * affine_velocity,
            -(rate**2) * tidal @ image + acceleration * affine_velocity,
        ]
        _, b, c = measure.optical_matrices(jets)
        np.testing.assert_allclose(-(c - acceleration / rate * b) / rate**2, tidal)
        data.append((b, c))
    np.testing.assert_allclose(
        (data[0][1] - data[1][1]) / (data[0][0] - data[1][0]), acceleration / rate
    )


@pytest.mark.parametrize("charge", [-1, 1])
def test_opposite_launches_and_same_source_deformation(charge):
    # Independent reverse propagation with reversed launch columns. Compare
    # analytic first variations against centered perturbations of the density.
    f = measure.beam.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = measure.beam.screens(board, route)
    state, directions, *_ = measure.beam.aperture.source_and_changes()
    rows = receipt("source_response.json")["species"][str(charge)]
    basis = np.c_[np.ones(5), measure.beam.BEARINGS]
    native = receipt("certificate.json")["reference"]["4.0"][str(charge)]
    for name, apertures in [("forward", rings), ("reverse", rings[::-1])]:
        effects = measure.signal_jets(board, route, apertures, charge, 4)

        def moments(density):
            return (
                basis.T
                @ np.array([measure.beam.probabilities(s, density) for s in effects])
                @ basis
            )

        def read(density):
            _, b, c = measure.optical_matrices(measure.matrix_jets(moments(density)))
            return np.array([b.trace(), -c.trace()])

        jets = measure.matrix_jets(moments(state))
        if name == "forward":
            det, b, c = measure.optical_matrices(jets)
            for value, key in [
                (det, "signed_area"),
                (b.trace(), "expansion"),
                (-c.trace(), "trace_focusing_if_affine"),
            ]:
                assert value in parse_interval(native[key])
        for key, d in [
            ("ordinary", directions["ordinary"]),
            ("unpair", directions["unpairing"]),
        ]:
            dm = moments(d)
            dj = measure.response_jets(moments(state), dm)
            db, dc = measure.deformation_response(jets, dj)
            delta = np.array([db.trace(), -dc.trace()])
            h = 1.0e-6
            for endpoint in (state + h * d, state - h * d):
                assert np.min(np.linalg.eigvalsh(endpoint)) > 0
            centered = (read(state + h * d) - read(state - h * d)) / (2 * h)
            np.testing.assert_allclose(delta, centered, rtol=2e-6, atol=1e-5)
            for value, field in zip(delta, ["expansion", "clock_focusing"]):
                assert value in parse_interval(rows[name][key][field])


def test_replay_and_exclusions():
    strict = receipt("certificate.json")
    sampled = receipt("counts.json")
    exact = receipt("results.json")
    for snapshot in exact:
        t = str(snapshot["dwell"])
        counts = sampled["records"][t]
        assert abs(counts["dwell_from_calibrated_clock"] - float(t)) < 1e-4
        for q in ["-1", "1"]:
            g = snapshot["charges"][q]
            ref = strict["reference"][t][q]
            assert g["signed_area"] in parse_interval(ref["signed_area"])
            assert g["trace_focusing_if_affine"] in parse_interval(
                ref["trace_focusing_if_affine"]
            )
            assert (parse_interval(counts[q]["signed_area_99pct"]).a > 0) == (
                g["signed_area"] > 0
            )
            assert all(
                sum(row) == sampled["shots_per_pixel"]
                for row in counts[q]["counts_by_sender"]
            )
            b = np.array(g["optical_deformation"])
            shear, twist = np.array(g["shear"]), g["twist"]
            np.testing.assert_allclose(
                np.trace(b @ b),
                g["expansion"] ** 2 / 2 + np.trace(shear @ shear) - 2 * twist**2,
            )
    for row in sampled["population_comparisons"]:
        q = str(row["charge"])
        true = exact[0]["charges"][q]["trace_focusing_if_affine"]
        errors = [abs(x["clock_focusing"] - true) for x in row["stencils"]]
        assert errors[0] > errors[1] > errors[2]
        assert true < 0 < row["zero_contrast_dwell_four"]["clock_focusing"]
    first, second = map(
        parse_interval, strict["necessary_common_nonaffinity_11_and_12"]
    )
    assert first.b < 0 < second.a
    assert parse_interval(strict["shared_trace_times_squared_affine_rate"]).b < 0
    reversal = receipt("source_response.json")
    assert (
        parse_interval(reversal["ordinary_calibration"]["signed_current_factor"]).a > 0
    )
    free = reversal["source_dependent_calibration"]
    assert 0 not in parse_interval(free["clock_and_material_conversion_minor"])
    assert parse_interval(free["reference_signed_current_factor"]).b < 0
    for directions in free["source_tangents"].values():
        du = directions["required_du_for_fixed_conversion"]
        assert parse_interval(du["-1"]).b < parse_interval(du["1"]).a
        assert parse_interval(directions["fixed_volume_survey_area_derivative"]).b < 0
        assert (
            parse_interval(directions["required_conversion_derivative_if_fitted"]).a > 0
        )
    for row in reversal["species"].values():
        assert 1 not in parse_interval(
            row["forward_minus_reverse"]["held_focusing_over_current_prediction"]
        )
    for row in reversal["species"].values():
        assert (
            parse_interval(
                row["forward_plus_reverse"][
                    "held_focusing_over_whole_seam_work_prediction"
                ]
            ).b
            < 1
        )
    bounds = reversal["odd_even_comparison"][
        "common_collar_rescaling_minimax_relative_error_lower"
    ]
    for key, cut in [("-1", 0.027547), ("1", 0.132880), ("both_charges", 0.142882)]:
        assert parse_interval(bounds[key]).a > cut
