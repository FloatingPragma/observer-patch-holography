from fractions import Fraction
from itertools import product
import json

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import Carrier
from . import measure


def test_sharp_moments_and_curvature_against_complete_native_evolution():
    board = Carrier("four-rail check", 4, [(0, 1), (0, 2), (1, 2), (2, 3)])
    words = list(product(range(3), repeat=board.n))
    lookup = {w: i for i, w in enumerate(words)}
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    k = candidate.KAPPA
    for rate in (1, 2):
        term = rate * exchange + k * contrast
        np.testing.assert_allclose(term[9:, :9], 0, atol=1e-13)
        h = np.zeros((len(words), len(words)))
        for i, word in enumerate(words):
            n = np.asarray(word)
            h[i, i] = k * ((6 - board.degree) @ (n**2 - 1))
            for a, b in board.seams:
                j = law.SRC.index((word[a], word[b]))
                for output, (x, y) in enumerate(law.SRC):
                    new = list(word)
                    new[a], new[b] = x, y
                    h[lookup[tuple(new)], i] += term[output, j]
            row = measure.sharp_moments(board, n, rate)
            np.testing.assert_allclose(2 * h[i, i], row["twice_mean"], atol=5e-13)
            residual = h[:, i].copy()
            residual[i] = 0
            np.testing.assert_allclose(
                4 * (residual @ residual), row["four_times_variance"], atol=1e-13
            )
        energy, states = eigh(h)
        delay = 1e-5
        evolution_minus_identity = (states * np.expm1(-1j * energy * delay)) @ states.T
        np.fill_diagonal(evolution_minus_identity, 0)
        departure = np.sum(abs(evolution_minus_identity) ** 2, axis=0)
        exact_curvature = np.array(
            [
                -measure.sharp_moments(board, np.asarray(w), rate)[
                    "four_times_variance"
                ]
                / 2
                for w in words
            ]
        )
        np.testing.assert_allclose(
            -2 * departure / delay**2, exact_curvature, rtol=1e-5, atol=2e-6
        )


def test_window_moment_bound_and_filled_sham():
    # The two-point law on {-1, 4} with weights {4/5,1/5} saturates Cantelli.
    mean, var = Fraction(0), Fraction(4)
    assert measure.window_bound(mean, var, 4, 5) == Fraction(1, 5)
    assert measure.window_bound(mean, var, -2, -1) == Fraction(4, 5)
    assert measure.window_bound(mean, var, -1, 1) == 1
    assert measure.window_bound(mean, 0, 4, 5) == 0
    result = measure.run()
    assert result == json.loads((measure.HERE / "results.json").read_text())
    assert (
        result["filled_sham"]["four_times_variance"]
        == result["filled_sham"]["twice_mean"]
        == 0
    )
    assert result["sharp_vacancy"]["four_times_variance"] == 5
    for panel in result["trials"].values():
        assert (
            panel["double_exchange_control"]["four_times_variance"]
            == 4 * panel["four_times_variance"]
        )
        assert (
            Fraction(panel["exact_rounded_target_weight_upper"])
            < Fraction(panel["one_unit_half_width_weight_upper"])
            < Fraction(1, 50)
        )
