"""Verify the measured ordinary source, the calibrated read and the sky certificates."""

from fractions import Fraction as F
from hashlib import sha256
from itertools import product

import numpy as np

from openwave.interval_echo.row_19.selected_collar import local_trapping as trap
from openwave.interval_echo.row_19.selected_collar import trapped_orbits as skytest

local = trap.local
steady = trap.steady


def test_protocol_and_catalogue_inputs():
    names = (
        "contact_census spectator_census material_polarization stationary_polarization "
        "full_board_polarization restored_balls local_trapping polarization_predictions "
        "polarization_sky trapped_orbits coherent_balls"
    )
    for name in names.split():
        for path, digest in local.inputs(name)["inputs"].items():
            assert sha256((local.HERE / path).read_bytes()).hexdigest() == digest


def test_vacancy_identifies_two_sources_degenerate_at_release():
    saved = local.inputs("contact_census")
    truth = trap.exp.native_pulse()
    assert truth == saved["truth"]
    for q, j in product(range(3), repeat=2):
        assert truth[0][q][j] == j
        assert truth[1][q][j] == q
    record, moved = (saved["contrasts"][k] for k in ("record", "redistribution"))
    a1, a2, beat = map(F, record["exact_release"])
    d = F(saved["extra_mark_gradient"])
    assert (a1, a2, beat) == (17 * d / 4, 65 * d / 4, -24 * d)
    assert moved["exact_release"] == record["exact_release"]
    assert moved["measured_vacancy"][1] < 0 < record["measured_vacancy"][0]
    off = local.inputs("spectator_census")
    assert off["crossing_matrix_elements_checked"] == 3768
    assert off["offdiagonal_elements"] == 3120


def test_early_total_source_is_compensated_at_every_size():
    rows = local.inputs("material_polarization")["full_boards"]
    assert {r["cells"] for r in rows} == {20, 80, 320}
    assert sum(len(r["shells"]) for r in rows) == 78
    for row in rows:
        values, err = row["density_by_cell"], row["per_cell_error"]
        assert abs(sum(values)) <= len(values) * err
        assert max(values[1:]) + err < 0 < values[0] - err
        for shell in row["shells"]:
            assert shell["enclosed"] > shell["enclosed_error"]
            assert shell["next_shell_mass"] < -shell["next_shell_error"]


def test_exact_energy_average_keeps_degenerate_coherence():
    energy = steady.EnergyAverage(np.diag(np.array([0, 0, 4], dtype=object)))
    source = np.ones((3, 3), dtype=object)
    kept, den = energy.pinch((source, 3))
    target = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 1]])
    np.testing.assert_array_equal(3 * kept, den * target)
    np.testing.assert_allclose(
        energy.finite_wait(source.astype(complex) / 3),
        target / 3,
        atol=float(energy.timing_error) + 1e-14,
    )


def test_stationary_mirror_reversal_and_monitoring():
    row = local.inputs("stationary_polarization")["mirror_matched"]
    assert 0 < F(row["full_monitoring_contraction_upper"]) < F(97, 100)
    for key, data in row["outcomes"].items():
        fwd, rev = (data["arms"][label] for label in ("forward", "reversed"))
        b = F(fwd["exact_available"])
        assert b == -F(rev["exact_available"]) == F(data["exact_reversal_odd"]) > 0
        assert fwd["available"][0] <= float(b) <= fwd["available"][1]
        assert rev["available"][0] <= -float(b) <= rev["available"][1]
        assert F(row["combined_trace_error_upper"]) < b / 10000
        if key == "unread":
            modes = {k: F(v) for k, v in fwd["available_by_restored_modes"].items()}
            assert sum(modes.values()) == b
            assert modes["gateway_visible"] == 0 < modes["gateway_dark"]


def test_full_board_two_probes_clock_and_energy():
    z = local.inputs("full_board_polarization")
    predictions = local.inputs("polarization_predictions")
    assert z["gates"] == [62, 229]
    assert 0 < z["source_contrast_score_error_upper"] < 1.5e-8
    assert F(z["maximum_extra_energy_over_ordinary_rest"]) == (
        F(z["arrival_excess"]) * z["kinetic_bandwidth_upper"] / z["rest_energy"]
    )
    for label, arm in z["outcomes"].items():
        b = arm["steady_available"]
        sign = 1 if label == "forward" else -1
        assert sign * b > z["source_contrast_score_error_upper"]
        lo, hi = arm["available_from_counts"]
        assert lo <= b <= hi and lo * hi > 0
        assert abs(sum(arm["steady_mass_by_cell"])) < 1e-13
        np.testing.assert_allclose(
            arm["predicted_release"], np.array([8.5, 16.25, -15.5]) * b
        )
        meter = arm["kinetic_meter"]
        lo, hi = meter["contrast_counts"]
        error = z["kinetic_bandwidth_upper"] * z["source_contrast_score_error_upper"]
        assert (
            0
            < lo
            < meter["predicted_energy"] - error
            < meter["predicted_energy"] + error
            < hi
        )
        for t, (held, predicted) in enumerate(
            zip(
                arm["held_flight_enclosed"],
                predictions["position_predictions"][label],
                strict=True,
            )
        ):
            for k, ((a, bhi), (c, d)) in enumerate(zip(held, predicted, strict=True)):
                assert max(a, c) <= min(bhi, d)
                direction = sign * (-1 if k == 2 else 1)
                assert min(direction * c, direction * d) > 0
                if t == 1 or k == 2:
                    assert min(direction * a, direction * bhi) > 0
    for settings in predictions["calibration_counts"]:
        for positions in settings:
            assert all(sum(v) == predictions["calibration_starts"] for v in positions)


def test_every_fixed_volume_comparison_excludes_a_common_gain():
    z = local.inputs("restored_balls")
    assert z["proper_clock_ratio"][0] < 209 / 240 < z["proper_clock_ratio"][1]
    assert (
        len(
            {(r["center"], tuple(r["hop_radii"]), r["direction"]) for r in z["surveys"]}
        )
        == 12
    )
    for row in z["surveys"]:
        a, b = row["common_gain_determinant"]
        assert a <= b and a * b > 0
        assert set(row["ball_sites"][0]) < set(row["ball_sites"][1])


def test_regional_interference_measures_the_reference_logarithm():
    z = local.inputs("coherent_balls")
    assert len(z["regions"]) == 16 and len(z["comparisons"]) == 12
    number = z["starts_per_setting"]
    for r in z["regions"].values():
        n = len(r["sites"])
        assert len(r["fresh_phase_counts"]) == n * (n - 1) // 2
        assert all(sum(k) == number for k in r["mode_and_outside_counts"].values())
        for pair in r["fresh_phase_counts"]:
            assert len(pair["quadratures"]) == 2
            assert all(sum(k) == number for k in pair["quadratures"])
        assert r["reference_eigenvalue_lower"] > 0
        assert (
            r["log_norm_error_upper"]
            >= r["coherence_norm_upper"] / r["reference_eigenvalue_lower"]
        )
    for row in z["comparisons"]:
        assert row["record_gain"][1] < 0
        a, b = row["common_gain_determinant"]
        assert a * b > 0
    target = z["comparisons"][0]["required_additional_entropy_contrast"]
    assert -145.69e-6 < target[0] < target[1] < -142.68e-6
    # A localized one-site sector and its vacuum have ordinary binary entropy.
    single = next(v for v in z["regions"].values() if len(v["sites"]) == 1)
    p = np.array(single["mode_and_outside_counts"]["fresh"], float) / number
    e = -p @ np.log(p.astype(float))
    assert single["entropy"][0] < e < single["entropy"][1]


def test_compact_mass_from_independent_stopped_paths():
    board = local.survey.bulk.federation(0)
    steps, hist, _ = trap.exp.flight.routes.layout(board)
    codes = trap.exp.flight.routes.receiver_codes(board, steps)
    internal = sorted(
        set(local.survey.routes.cell().seams)
        - set(local.survey.routes.codebook()[0][0])
    )
    chosen = 5
    region, _, proj, scores = trap.arrangements(board, steps, codes, internal)[chosen]
    # Enumerate path mass by current cell and four visit parities. This uses
    # only forward probabilities, independently of the saved resolvent.
    visited = np.zeros((board.cells, 16))
    visited[0, 0] = 1
    mask = {v: 1 << (3 - i) for i, v in enumerate(region)}
    indices = np.array([np.arange(16) ^ mask.get(c, 0) for c in range(board.cells)])
    parity = np.zeros(16)
    for _ in range(1024):
        after = np.take_along_axis(visited, indices, axis=1)
        parity += after.sum(axis=0) / 32
        visited = (31 / 384) * hist.T @ after
    assert abs(parity.sum() - 1) < 2e-14
    assert np.trace(proj) == 3
    exact = local.inputs("local_trapping")["surveys"][0]["exact_cell_t6"][chosen]
    np.testing.assert_allclose(
        parity @ np.array(scores, float) / (9 * board.n * 144),
        float(F(exact)),
        atol=1e-19,
        rtol=1e-12,
    )


def test_trapping_is_positive_at_finite_time_and_conserves_all_reads():
    z = local.inputs("local_trapping")
    for r in z["surveys"]:
        assert len(r["compact_dimensions"]) == r["cells"]
        assert set(r["compact_dimensions"]) == ({3} if r["cells"] == 20 else {3, 4, 5})
        assert F(r["all_cell_excess_probability_lower"]) > 0
        assert all(
            x > e for x, e in zip(r["per_cell_t6"], r["per_cell_error"], strict=True)
        )
    base = z["surveys"][0]
    assert len(base["unit_time_path_witnesses"]) == 20
    assert (
        min(a["unconditional_excess_lower"] for a in base["unit_time_path_witnesses"])
        > 0
    )
    sampled = base["accepted_trajectory_sample"]
    assert (
        sum(sampled["meter_counts"])
        == sampled["paths"] * sampled["shots_per_path"]
        == 2**33
    )
    assert all(x > 0 for x in sampled["meter_counts"])


def test_sparc_exclusions_and_fitted_witness():
    smooth = local.inputs("polarization_sky")
    trapped = local.inputs("trapped_orbits")
    assert smooth["summary"]["ring_count"] == 1288
    assert (
        len(smooth["summary"]["excluded"]) == 81
        and smooth["summary"]["compatible"] == 55
    )
    assert trapped["summary"] == dict(
        leading_compatible=93, finite_envelope_excludes=25, galaxies=136
    )
    profile = list(map(F, trapped["cumulative_leading"]))
    assert 0 < profile[0] and sorted(set(profile)) == profile and profile[-1] == 1
    sky = skytest.sky
    raw = (sky.sample.measure.DATA / "table2.dat").read_text().splitlines()
    for name in ("DDO064", "DDO154"):
        row = trapped["results"][name]
        rings = sky.rings_in(raw, [{"line": line} for line in row["lines"]])
        if row["leading_compatible"]:
            params = list(map(F, row["parameters"]))
            shell = dict(
                zip(sorted({r[0] for r in rings}), row["shell_indices"], strict=True)
            )
            for radius, lo, hi, bar, _ in rings:
                estimate = (
                    sum(x * y for x, y in zip(bar, params[:3], strict=True))
                    + params[3] * profile[shell[radius]] / radius
                )
                assert lo <= estimate <= hi
        else:
            assert row["finite_envelope_excludes"]
            assert all(F(c["margin"]) > 0 for c in row["leading_rejections"])
            assert len(row["leading_rejections"]) == sum(
                1 for _ in skytest.embeddings(sorted({r[0] for r in rings}))
            )
