"""Check the reduced signed record against unreduced native seam evolution."""

import gzip
import json

import numpy as np
import pytest
from scipy.linalg import eigh

from model.executable_model.carrier import cell
from ...row_31_orbitals.occupied_write import core
from ...row_31_orbitals.occupied_write.test_occupied import raw_hamiltonian
from . import flight, measure


def receipt():
    return json.loads(gzip.decompress(measure.OUTPUT.read_bytes()))


@pytest.mark.parametrize("strength", [0, 8])
def test_full_132_state_sign_reader_and_support(strength):
    board = cell()
    seam = (0, 1)
    native = raw_hamiltonian(board, seam, strength)
    dim = board.n
    source = np.zeros((11, dim))
    wrote = flight.preparation().reshape(11, 9)
    source[:, 0], source[:, 1] = wrote[:, 3], wrote[:, 1]
    keep = np.array([True] + [abs(d) <= 1 for c, d in core.LAW.LABELS]).repeat(dim)
    np.testing.assert_allclose(native[np.ix_(keep, ~keep)], 0, atol=2e-12)
    spectrum = 6 * strength + eigh(native[np.ix_(keep, keep)], eigvals_only=True)
    assert spectrum.min() >= 5 * strength - 1e-11 and spectrum.max() <= 7 * strength + 7
    e, u = eigh(native)
    times = np.array([0, 0.117, np.pi, 73 * np.pi, 400 * np.pi])
    state = ((np.exp(-1j * times[:, None] * e) * (u.T @ source.ravel())) @ u.T).reshape(
        len(times), 11, dim
    )
    direction = np.array([False] + [d == 1 for c, d in core.LAW.LABELS])
    direct = np.sum(abs(state[:, direction, :]) ** 2, axis=(1, 2))
    got, _ = flight.clock_trace(board, seam, strength, times)
    np.testing.assert_allclose(got, direct, atol=3e-10, rtol=0)
    # The written orbital reference in raw_hamiltonian removed the uniform-load energy.
    if strength == 8:
        form = measure.formation(board, seam, flight.preparation(), measure.Sampler())
        actual = 48 + source.ravel() @ native @ source.ravel()
        assert abs(actual - form["operator_mean"]) < 1e-11
        assert abs(actual - form["reading"]) < form["simultaneous_sampler_halfwidth"]


def test_all_law_probabilities_and_independent_binomial_totals():
    d = receipt()
    core.integer_blocks()
    assert len(d["sampling_weights_over_2_to_34"]) <= 160
    rng = np.random.Generator(np.random.PCG64(d["seed"]))
    expected = {}
    for key, weights in d["sampling_weights_over_2_to_34"].items():
        expected[key] = rng.multinomial(
            d["shots_per_setting"], np.array(weights) / 2**34
        ).tolist()
    for size, panel in d["boards"].items():
        native = panel["flight"]["native"]
        for label, run in panel["flight"].items():
            banks = run["banks"]
            for i, t in enumerate(d["ticks"]):
                p = banks["positive"][i]["calculated"][0]
                assert abs(banks["negative"][i]["calculated"][0] - (1 - p)) < 1e-12
                assert banks["dephased"][i]["calculated"] == [0.5, 0.5, 0]
                for source, bank in banks.items():
                    row = bank[i]
                    assert row["counts"] == expected[f"{size}_{label}_{source}_{t}"]
                    assert row["counts"][2] == 0 and sum(row["counts"]) == measure.SHOTS
                    assert (
                        abs(row["counts"][0] / measure.SHOTS - row["calculated"][0])
                        < measure.RADIUS
                    )
        assert panel["full_swing"] and panel["fitted"]["held_pass"]
        assert panel["fitted"]["held_max_error"] < 0.0013
        gear = native["ground_audit"]
        assert abs(panel["fitted"]["angular_frequency"] - gear["gap"]) < 4e-6
        assert min(gear["ground_weights"]) > 0.994
        assert gear["signed_gap_coefficient"] > 0.994
        assert max(gear["ground_residuals"]) < 1e-11
        for i, r in enumerate(panel["formation"]["components"]):
            assert r["counts"] == expected[f"formation_{size}_{i}"]
    for i in range(3):
        assert d["planted_signed_record"][i]["counts"] == expected[f"planted_{i}"]
        assert d["timer_calibration"][i]["counts"] == expected[f"timer_{i}"]


def test_degree_bounds_all_three_federations_and_signed_source():
    for level in range(3):
        board, seam, _ = core.chart(level)
        assert min(board.degree) == 5 and max(board.degree) == 6
        n = board.n
        for odd in (False, True):
            h, chi = flight.block(board, seam, 8, odd)
            np.testing.assert_allclose(
                chi @ (h @ chi) + 48, 40 + (sum(board.degree[list(seam)]) - 2) / 4
            )
            assert np.flatnonzero(chi).tolist() == [seam[1], n + seam[0]]
            # Cauchy on each edge bounds the quadratic energy of L/2 by max degree.
            assert h.shape == (2 * n, 2 * n)
