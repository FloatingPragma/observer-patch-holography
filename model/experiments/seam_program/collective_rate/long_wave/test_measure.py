"""Replay the two physical protocols and independently check the integer commutator."""

import hashlib
import json
from fractions import Fraction

import numpy as np

from model.candidate_laws.interval_echo.candidate import native_law
from model.executable_model.carrier import Carrier

from model.experiments.seam_program.collective_rate.long_wave import measure


def test_local_expression_against_the_native_full_hilbert_space():
    board = Carrier(
        "irregular", 6, [(0, 1), (0, 5), (1, 5), (1, 2), (2, 3), (3, 4), (4, 5)]
    )
    z = [-3, 5, -9, 4, 1, 2]
    g = np.array(z) / np.linalg.norm(z)
    labels = [native_law.SRC.index(x) for x in ((0, 0), (0, 1), (1, 0), (1, 1))]
    masks = np.arange(1 << board.n)
    bits = (masks[:, None] >> np.arange(board.n)) & 1
    observable = bits @ g
    f = measure.BASE
    dist = f ** bits.sum(1) * (1 - f) ** (board.n - bits.sum(1))
    a, b = measure.clock_seam()
    # The whole coherent perturbation has zero partial trace over either
    # endpoint, including correlations with every port on the other side.
    seam = np.zeros((len(masks), len(masks)))
    for i in masks[(bits[:, a] == 1) & (bits[:, b] == 0)]:
        j = i ^ (1 << a) ^ (1 << b)
        seam[i, j] = seam[j, i] = -dist[i] / (f * (1 - f))
    assert np.all(seam.diagonal() == 0)
    for site in (a, b):
        outside = masks[bits[:, site] == 0]
        np.testing.assert_array_equal(
            seam[np.ix_(outside, outside)]
            + seam[np.ix_(outside + (1 << site), outside + (1 << site))],
            0,
        )
    driven = np.zeros_like(seam)
    for site in range(board.n):
        lo = masks[bits[:, site] == 0]
        hi = lo + (1 << site)
        reduced = seam[np.ix_(lo, lo)] + seam[np.ix_(hi, hi)]
        driven[np.ix_(hi, hi)] += g[site] * reduced
        driven[np.ix_(lo, lo)] -= g[site] * reduced
    for kappa in (0, int(measure.KAPPA)):
        block = native_law.seam_term(kappa)[np.ix_(labels, labels)]
        h = np.diag(kappa * (bits @ (6 - board.degree))).astype(float)
        for i, j in board.seams:
            word = bits[:, i] * 2 + bits[:, j]
            for pos in range(4):
                rows = (
                    masks
                    ^ ((bits[:, i] ^ (pos >> 1)) << i)
                    ^ ((bits[:, j] ^ (pos % 2)) << j)
                )
                h[rows, masks] += block[pos, word].real
        assert abs(np.sum(h * seam) - 1) < 2e-15
        first = h * (observable[None, :] - observable[:, None])
        second = h @ first - first @ h
        m2, polynomial, local = measure.coefficients(board, z, kappa)
        chi = f * (1 - f)
        np.testing.assert_allclose(
            np.sum(dist * (first**2).sum(1)) / chi, float(m2), rtol=2e-14
        )
        np.testing.assert_allclose(
            np.sum(dist * (second**2).sum(1)) / chi,
            sum(float(v) * f**i for i, v in enumerate(polynomial)),
            rtol=2e-14,
        )
        np.testing.assert_allclose(
            -np.sum(second * driven), float(local), rtol=2e-13, atol=1e-14
        )


def test_registered_long_modes_and_exact_changes():
    archive = json.loads((measure.HERE / "results.json").read_text())
    assert (
        hashlib.sha256((measure.HERE / "profiles.json").read_bytes()).hexdigest()
        == "cbefbfe8955f4be658e47c7f556a9b3a666545b659addb95b169e54e6a3c23ee"
    )
    boards = [measure.cell(), *(measure.federation(i) for i in range(3))]
    last = 4
    for board, frozen in zip(boards, archive["boards"], strict=True):
        live = measure.compare_board(board)
        p = live["profile"]
        assert Fraction(p["exact_residual_squared"]) < Fraction(1, 10**9)
        eigen = float(Fraction(p["exact_rayleigh"]))
        assert 0 < eigen < last and p["selector"]["residual"] < 1e-9
        last = eigen
        assert p["exact_rayleigh"] == frozen["profile"]["exact_rayleigh"]
        assert live["by_contrast"] == frozen["by_contrast"]
        active = live["by_contrast"]["8"]
        control = live["by_contrast"]["0"]
        appointment = active["finite_appointment"]
        assert Fraction(appointment["dwell"]) > 0
        for name, interval in appointment["absolute_signal_differences"].items():
            lo, hi = map(Fraction, interval)
            assert (lo < hi < 0) if name == "collar_driven" else (0 < lo < hi)
        tangents = {
            k: Fraction(v) for k, v in active["mixed_reference_tangents"].items()
        }
        assert tangents["ordinary_passive"] > 0 and tangents["ordinary_driven"] > 0
        assert (
            tangents["fixed_marginal_passive"] > 0
            and tangents["fixed_marginal_driven"] < 0
        )
        assert tangents["observed_over_predicted"] < -1
        assert (
            tangents["calibrated_prediction"]
            == tangents["ordinary_driven"]
            * tangents["fixed_marginal_passive"]
            / tangents["ordinary_passive"]
        )
        a = {k: Fraction(v) for k, v in active["comparisons"].items()}
        assert a["ordinary_passive_change"] == a["ordinary_driven_change"] > 0
        assert (
            a["fixed_marginal_driven_change"] < 0 < a["fixed_marginal_passive_change"]
        )
        assert (
            a["localized_seam_driven_second_time_derivative"]
            != a["localized_seam_passive_change"]
            == 0
        )
        assert (
            abs(control["readings"]["fixed_marginal_driven_change"])
            < abs(float(a["fixed_marginal_driven_change"])) / 50
        )
        if board.n >= 240:
            assert active["baseline_fourth_over_second_squared"] > 2e4


def test_cell_spectrum_independent_flights_and_finite_counts():
    data = measure.run()
    archive = json.loads((measure.HERE / "results.json").read_text())
    read = data["common_stationary_seam_read"]
    assert sum(read["source_counts"]) == read["trials"]
    lo, hi = read["independent_99pct_band"]
    assert lo < read["signed_modular_slope"] < hi < 0
    shots = 2**50
    for kappa, live in data["cell_time_records"].items():
        frozen = archive["cell_time_records"][kappa]
        assert live["count_records"] == frozen["count_records"]
        envelope = data["boards"][0]["by_contrast"][kappa]
        second = float(Fraction(envelope["normalized_second_moment"]))
        poly = [float(Fraction(x)) for x in envelope["normalized_fourth_polynomial"]]
        for f, (m2, m4) in zip(
            live["product_fillings"], live["direct_spectral_moments"], strict=True
        ):
            np.testing.assert_allclose(m2, second, rtol=5e-13)
            np.testing.assert_allclose(
                m4, np.polynomial.polynomial.polyval(f, poly), rtol=5e-13
            )
        for name in ("reference", "ordinary"):
            pair = live["sources"][name]
            np.testing.assert_allclose(
                pair["passive_normalized"], pair["driven"], rtol=0, atol=1e-14
            )
        counted = live["count_records"]
        for click, expected in zip(
            counted["dimer_positive_counts"], live["material_clock_clicks"], strict=True
        ):
            assert (
                abs(click / shots - expected) < counted["click_probability_halfwidth"]
            )
        for name, arm in counted["observations"].items():
            panel = live["sources"].get(name, live["sources"]["reference"])
            for meter, values in arm.items():
                normalized = np.array(
                    panel["passive_normalized" if meter == "passive" else "driven"]
                )
                if name.startswith("local") and meter == "driven":
                    sign = 1 if name == "local_plus" else -1
                    normalized += (
                        sign
                        * counted["local_signed_step"]
                        * np.array(live["localized_seam_driven_change"])
                    )
                assert (
                    np.max(np.abs(np.array(values["reconstructed"]) - normalized))
                    < values["simultaneous_99pct_halfwidth"]
                )
        np.testing.assert_allclose(
            live["localized_seam_driven_change"],
            frozen["localized_seam_driven_change"],
            rtol=1e-11,
            atol=1e-13,
        )
        if kappa == "8":
            rows = counted["observations"]
            for meter, sign in (("passive", 1), ("driven", -1)):
                held = rows["fixed_marginals"][meter]
                base = rows["reference"][meter]
                assert (
                    sign * (held["reconstructed"][1] - base["reconstructed"][1])
                    > held["simultaneous_99pct_halfwidth"]
                    + base["simultaneous_99pct_halfwidth"]
                )

            for prefix in ("ordinary", "collar"):
                for meter in ("passive", "driven"):
                    high = rows[prefix + "_upper"][meter]
                    low = rows[prefix + "_lower"][meter]
                    sign = -1 if (prefix, meter) == ("collar", "driven") else 1
                    assert sign * (
                        high["reconstructed"][1] - low["reconstructed"][1]
                    ) > (
                        high["simultaneous_99pct_halfwidth"]
                        + low["simultaneous_99pct_halfwidth"]
                    )
