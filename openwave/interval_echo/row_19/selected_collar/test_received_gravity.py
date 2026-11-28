"""Independent checks on the serial source, exterior dynamics, geometry and sky bounds."""

from fractions import Fraction as F
from hashlib import sha256
from math import log, sqrt

import numpy as np
from scipy.sparse import coo_matrix

from openwave.interval_echo.row_19.selected_collar import (
    receiver_stationarity as station,
)
from openwave.interval_echo.row_19.selected_collar import routed_contacts as routes

local = routes.local


def check_inputs(row):
    for path, digest in row["inputs"].items():
        assert sha256((local.HERE / path).read_bytes()).hexdigest() == digest


def test_serial_source_is_native_two_contact_coincidence():
    z = local.inputs("shared_contact")
    check_inputs(z)
    truth = np.array(z["serial_route"]["sharp_arrivals"])
    word = np.array([1, 1, 0, 0, 1, 1, 1, 0, 1])
    np.testing.assert_array_equal(truth, np.outer(word, word))
    fresh, shared = [z["relay"]["serial"][a] for a in ("fresh", "shared")]
    assert F(shared["probability"]) - F(fresh["probability"]) == F(1, 144)
    assert shared["read"]["mean"][0] > fresh["read"]["mean"][1]
    assert F(z["relay"]["0"]["shared"]["probability"]) < F(
        z["relay"]["0"]["fresh"]["probability"]
    )


def test_odd_visits_including_stopping_visit_and_full_joint_counts():
    z = local.inputs("routed_contacts")
    check_inputs(z)
    for row in z["profiles"]:
        g = np.array(row["green_from_source"])
        d = np.array(row["diagonal_green"])
        u = np.array(row["odd_visit_probabilities"])
        np.testing.assert_allclose((2 * d - 1) * u, g, rtol=1e-14)
        near, far = row["held_edge"]
        clock = row["serial_near_minus_far_frequency_per_unit_acceptance"]
        assert abs(clock + 24 * (u[near] - u[far])) < 1e-14
        assert (
            abs(clock / 120 - row["serial_fraction_of_blank_beat_per_unit_acceptance"])
            < 1e-14
        )
        if row["stations"] != 20:
            assert row["raw_exterior_equation_max_abs"] > 1e-5
            continue
        a, b = row["held_edge"]
        prob = np.vectorize(lambda x: float(F(x)))(row["joint_parity_and_stop"])
        np.testing.assert_allclose([prob[1].sum(), prob[:, 1].sum()], u[[a, b]])
        np.testing.assert_allclose(
            prob.sum(axis=(0, 1))[:2], g[[a, b]] / 2 ** row["stop_bits"]
        )
        np.testing.assert_allclose(
            np.array(row["sampled_paths"]["two_receivers"]) / z["starts"],
            prob,
            atol=sqrt(log(2 * 12 / 0.01) / (2 * z["starts"])),
            rtol=0,
        )


def test_inward_free_launch_and_ordinary_scale():
    z = local.inputs("routed_contact_flight")["serial_agreement"]
    check_inputs(z)
    arms = z["arms"]
    delta = list(map(F, z["overlap_change"]))
    assert delta[0] > 0 and delta[1] - 2 * delta[0] > 0 and delta[2] < 0
    for i in (0, 1):
        assert (
            F(arms["shared"]["in_phase_derivatives"][i])
            > F(arms["fresh"]["in_phase_derivatives"][i])
            > F(arms["no_reports"]["in_phase_derivatives"][i])
            == 0
        )
    for t in (0, 1):
        d1, d2, dc = z["ordinary_doses"]["ordinary_one"]["finite_dose_intervals"][t]
        assert 0 < d1[0] < d1[1] < d2[0] < d2[1] < dc[0] < dc[1]
        # The in-phase single launch itself records extra inward motion.
        for q in (0, 1):
            assert (
                arms["shared"]["reads"][t]["position_counts"][q][0]["mean"][0]
                > arms["fresh"]["reads"][t]["position_counts"][q][0]["mean"][1]
            )
    for tag in ("ordinary_one", "ordinary_two"):
        v = list(map(F, z["ordinary_doses"][tag]["reference_change"]))
        assert v[1] <= 2 * v[0]


def test_source_only_finite_predictions_are_inside_independent_held_intervals():
    row = local.inputs("source_predictions")
    check_inputs(row)
    assert row["series_order"] == 26
    release = row["release_from_measured_source_and_route"]
    for panel in row["times"]:
        for k, (pred, observed) in enumerate(
            zip(panel["source_only"], panel["held"], strict=True)
        ):
            assert observed[0] < pred[0] <= pred[1] < observed[1]
            assert (
                (min(release[k]) > 0 and min(observed) > 0)
                if k < 2
                else (max(release[k]) < 0 and max(observed) < 0)
            )


def test_arbitrary_passive_mark_layouts_on_all_native_bridges():
    # Direct integer double commutators, retaining side counts on the whole
    # graph. The shape of these controls has no radial symmetry requirement.
    rng = np.random.default_rng(91930)
    for level in (0, 1, 2):
        board = local.survey.bulk.federation(level)
        exterior = set(map(tuple, board.external))
        for probability in (0.0, 0.31, 0.69):
            blank = np.array(
                [e for e in board.seams if e in exterior or rng.random() >= probability]
            )
            i, j = blank.T
            active = np.bincount(np.r_[i, j], minlength=board.n)
            marks = np.array(board.degree) - active
            for q in (1, 2):
                diag = active + 16 * q * q * (6 - marks)
                h = coo_matrix(
                    (
                        np.r_[diag, -np.ones(2 * len(i), dtype=int)],
                        (
                            np.r_[np.arange(board.n), i, j],
                            np.r_[np.arange(board.n), j, i],
                        ),
                    ),
                    shape=(board.n, board.n),
                ).tocsr()
                for a, b in board.external:

                    def column(k):
                        lo, hi = h.indptr[k : k + 2]
                        return dict(zip(h.indices[lo:hi], h.data[lo:hi], strict=True))

                    left, right = column(a), column(b)
                    assert left.keys() & right.keys() == {a, b}
                    plus = 0
                    for v in left.keys() | right.keys():
                        ha, hb = left.get(v, 0), right.get(v, 0)
                        plus += (ha + hb) ** 2 * (v // 12 == a // 12) - (ha + hb) * ha
                    # Four times native curvature. The two incident gateway
                    # ports each have original coordination six.
                    assert plus == (1 + 16 * q * q) * (marks[a] - marks[b])


def test_material_fixed_volume_tangent_and_finite_difference():
    z = local.inputs("receiver_stationarity")
    check_inputs(z)
    for row in z["surveys"]:
        assert len(row["eligible_tetrahedra"]) == 26
        assert sum(row["face_pair_counts"].values()) == 46
        assert row["face_pair_counts"]["positive_area_positive_entropy"] == 18
        for target, clock in zip((1, 0.8), row["sharp_clocks"], strict=True):
            lo, hi = clock["relative_frequency"]
            assert lo < target < hi
            assert all(
                sum(read["counts"]) == z["starts_per_acquisition"]
                for read in clock["reads"]
            )
        c = row["first_positive_pair"]
        for key in (
            "area_tangent",
            "entropy_tangent_nats",
            "delta_area",
            "delta_entropy_nats",
        ):
            lo, hi = c["measured"][key]
            assert 0 < lo < c["ideal"][key][0] <= c["ideal"][key][1] < hi
        a, b = c["stations"]
        assert len(set(a) & set(b)) == 3
        assert station.gram_data(a, np.array(row["earliest_swap_counts"][0]))[0] == 256


def test_outer_orbit_integral_bounds_and_conditional_full_dose():
    z = local.inputs("contact_lapse_sky")
    check_inputs(z)
    rows = [r for r in z["galaxies"].values() if "weak_log_lapse_span" in r]
    assert len(rows) == 128 and len(z["galaxies"]) == 136
    for row in rows:
        r = list(map(F, row["radii_kpc"]))
        lo = list(map(F, row["cumulative_lower"]))
        hi = list(map(F, row["suffix_upper"]))
        assert all(x <= y for x, y in zip(lo, hi, strict=True))
        assert sorted(lo) == lo and sorted(hi) == hi
        span = list(map(F, row["weak_log_lapse_span"]))
        for values, actual in ((lo[:-1], span[0]), (hi[1:], span[1])):
            integral = sum(
                x * (1 / a - 1 / b)
                for x, a, b in zip(values, r[:-1], r[1:], strict=True)
            )
            assert actual == integral / F("299792.458") ** 2
        assert float(span[1]) < z["native_log_lapse_full_dose"][0]
        assert 0 < row["maximum_embedded_dose_upper"] < 1


def test_coherent_shells_and_disjoint_ordinary_response():
    z = local.inputs("shell_receivers")
    check_inputs(z)
    assert z["native_two_load_comparisons"] == 400
    shells = z["shells"]
    assert len({(r["stations"], r["stop_bits"], r["shell"]) for r in shells}) == 78
    outward = []
    for row in shells:
        assert row["record_gradient"] > row["record_error"]
        assert all(p["amplitude_error"] < 2e-13 for p in row["launches"])
        a, b = row["record_curvatures"]
        assert b > 2 * a > 0
        ordinary = row["ordinary_gradient"]
        assert abs(ordinary) > row["ordinary_error"]
        for tag in ("ordinary_one_curvatures", "ordinary_two_curvatures"):
            c, d = row[tag]
            if ordinary > 0:
                assert 0 < d <= 2 * c
            else:
                assert max(c, d) < 0
        if ordinary < 0:
            outward.append((row["stations"], row["stop_bits"], row["shell"]))
    assert outward == [(20, 5, 4)]
