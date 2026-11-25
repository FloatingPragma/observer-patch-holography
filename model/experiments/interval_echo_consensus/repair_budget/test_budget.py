from fractions import Fraction
from itertools import combinations
import json
from math import comb

import numpy as np

from .measure import (
    HERE,
    SHOTS,
    forecast,
    native_patch,
    operators,
    rates,
    terminal_budget,
)
from model.candidate_laws.interval_echo.candidate import native_law


def test_wait_law_against_enumeration_of_all_five_offer_words():
    words = {(2, 0, 2, 0): 1}
    pairs = tuple(combinations(range(4), 2))
    for t in range(6):
        distribution = [
            sum(n for w, n in words.items() if w.count(2) == 2 - k) / 6**t
            for k in range(3)
        ]
        np.testing.assert_allclose(
            forecast(4, 2, 2, t), distribution, rtol=0, atol=2e-15
        )
        after = {}
        for word, n in words.items():
            for a, b in pairs:
                w = list(word)
                if abs(w[a] - w[b]) == 2:
                    w[a] = w[b] = 1
                key = tuple(w)
                after[key] = after.get(key, 0) + n
        words = after


def test_scalar_effect_on_coherences_and_actual_99_state_free_law():
    _, _, words, idx, _, echo = native_patch()
    ops = operators(words, idx, echo, 1)
    effect = sum(g.T @ g for _, jumps in ops for g in jumps) / len(ops)
    expected = [w.count(0) * w.count(2) / len(ops) for w in words]
    np.testing.assert_allclose(effect, np.diag(expected), rtol=0, atol=2e-15)
    ex, co = native_law.local_term_parts()
    for n in (0, 2):
        counts = np.tile([int(a == n) + int(b == n) for a, b in native_law.SRC], 11)
        for h in (ex, co):
            np.testing.assert_allclose(
                (counts[:, None] - counts[None, :]) * h, 0, atol=2e-15
            )


def test_bound_exact_at_endpoints_and_dominates_small_complete_chains():
    for ports in (3, 4, 7, 12):
        for v in range(ports + 1):
            for d in range(ports - v + 1):
                ps = rates(ports, v, d)
                result = terminal_budget(ports, v, d)
                assert Fraction(result["mean"]) < 2 * comb(ports, 2)
                assert forecast(ports, v, d, result["cap"])[-1] >= 0.99 - 1e-11
                if not ps:
                    continue
                m = min(v, d)
                a = abs(v - d)
                factor = Fraction(1)
                for p in ps[:-1]:
                    factor *= p / (p - ps[-1])
                assert factor == Fraction(m * (a + 2), m + a + 1)
                for t in (1, 3, 8, 27):
                    assert (
                        sum(forecast(ports, v, d, t)[:-1])
                        <= float(factor * (1 - ps[-1]) ** t) + 1e-13
                    )


def test_all_forecast_decisions_from_archived_counts():
    doc = json.loads((HERE / "counts.json").read_bytes())
    for row in doc["banks"]:
        freq = np.array(row["writes"]) / SHOTS
        assert sum(row["writes"]) == SHOTS
        accepted = (
            np.max(np.abs(freq - forecast(4, 2, 2, row["attempts"])))
            <= doc["simultaneous_half_width"]
        )
        assert accepted == row["accepts_native_forecast"]
        if row["mode"] == "native":
            assert accepted
        if row["mode"] != "native" and row["attempts"] in (5, 9, 16):
            assert not accepted
