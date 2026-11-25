"""Native word rotations, conjugate gating, qutrit phase records and stock balance."""

import json
import numpy as np

from .apparatus import stock
from .measure import HERE, SHOTS, WORDS, gate


def test_exact_domain_for_coherent_half_pulses_and_unsuitable_stock():
    for sign in (-1, 1):
        g = 1 - sign
        c = WORDS.index((1 + sign, 1, g))
        n = WORDS.index((1 + sign, g, 1))
        ix = np.ix_([c, n], [c, n])
        for angle in (0.0, 0.437, np.pi / 2, np.pi):
            a = np.cos(angle / 2)
            b = np.sin(angle / 2)
            np.testing.assert_allclose(
                gate(angle, g)[ix], [[a, -b], [b, a]], atol=2e-13
            )
        half = gate(np.pi / 2, g)
        np.testing.assert_allclose(
            [(abs(half @ half) ** 2)[n, c], (abs(half) ** 2 @ abs(half) ** 2)[n, c]],
            [1, 0.5],
            atol=2e-13,
        )
        assert abs(gate(np.pi, 2 - g)[n, c]) < 2e-13
        assert stock(sign)["target"] == 167


def test_actual_converted_ternary_word_and_qutrit_records_for_each_sign():
    r = json.loads((HERE / "results.json").read_text())
    for experiment in r["experiments"]:
        for basis, i, counts in experiment["load_banks"]:
            assert sum(counts) == SHOTS
            live = [r["ternary_word_order"][j] for j, k in enumerate(counts[:27]) if k]
            if experiment["setting"] == "full":
                assert live == [experiment["neutral_word"]]
            elif experiment["setting"] == "off":
                assert live == [experiment["charged_word"]]
            else:
                assert len(live) == 2
        for _, _, _, c in experiment["counts"]:
            assert sum(c) == SHOTS
        if experiment["setting"] != "off":
            assert (
                experiment["neutral_qutrit_witness"]["entanglement_fidelity_lower_99"]
                > 2 / 3
            )
