import json
from itertools import product

import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm

from . import measure as m
from . import receiver as r


def test_native_three_coordinate_filter_and_free_clock_control():
    for q in (1, 2):
        for k in (4, 8):
            diag, _ = r.terms(q, k)
            for s in (0, 1, 2):
                assert diag[q, s] == k * (q - s) ** 2 + (q != s) / 2
            for sign in (-1, 1):
                for word in product((0, 1), repeat=3):
                    for t in (r.FLIGHT, 0.5 * r.FLIGHT):
                        p = np.mean(
                            r.probabilities(q, k, sign, word, (0, 1, 2), t, False)
                        )
                        d = q * (4 * k + 0.5)
                        z = ((1 + np.exp(1j * d * r.CLOCK)) / 2) ** (3 - sum(word))
                        a = r.exact_coefficient(q, k, t) * z.real
                        b = np.sin(t * np.hypot(1, d)) * z.imag / (2 * np.hypot(1, d))
                        assert_allclose(p, 0.5 + sign * (a + b), atol=2e-14, rtol=0)
                        if k == 8 and q == 1 and not all(word):
                            assert_allclose(p, 0.5, rtol=0, atol=1e-14)
    assert len(r.pad_routes()) == 36


def test_readout_uses_source_green_difference_and_keeps_every_control():
    data = json.loads((m.HERE / "results.json").read_text())
    states = np.zeros((10,) * 3)
    for c in m.green.CORNERS:
        states[tuple(i - 1 for i in c)] = 1 / 8
    u = m.green.poisson(states, (1 / 8,) * 3)[1:-1, 1:-1, 1:-1]
    assert len(data["panels"]) == 14 and data["last_drain"] <= m.H
    for record, config in zip(data["panels"], m.configurations(), strict=True):
        name, face, q, k, axes, time, parked = config
        predicted, bound = m.forecast(face, q, k, axes, time, parked, u)
        shift = (record["right_reads"] - data["all_receiver_reads_per_panel"] / 2) / m.N
        assert record["setting"] == name and record["face"] == face
        assert record["shift_per_release"] == shift
        assert_allclose(record["predicted"], predicted, rtol=0, atol=1e-12)
        assert abs(shift - predicted) <= m.confidence(bound)
        if name == "matched":
            x = np.array(face[:3]) - 1
            y = x.copy()
            y[face[-1]] += 1
            assert_allclose(
                predicted, 65 / 4229 * (u[tuple(y)] - u[tuple(x)]), atol=1e-13
            )


def test_full_board_free_source_matches_counts():
    saved = json.loads((m.HERE / "results.json").read_text())["free_flight"]
    board, _, _ = r.rails.layout()
    h = np.diag(board.degree * 0.5) - board.adjacency().toarray() * 0.5
    p = abs(expm(-3j * h)[:, r.LEFT]) ** 2
    for bank in saved["banks"]:
        counts = np.zeros(board.n)
        counts[bank["where"]] = bank["counts"]
        assert sum(counts) == bank["starts"]
        assert np.max(abs(counts / bank["starts"] - p)) < np.sqrt(
            np.log(4 * board.n / 0.01) / (2 * bank["starts"])
        )
