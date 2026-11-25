from fractions import Fraction as Q
from itertools import product
import json
from math import factorial

import numpy as np

from . import measure as m


def test_sector_compression_and_equal_energy_sources():
    rail, pins, board = m.apparatus()
    law = m.mediator.candidate.native_law
    h = law.seam_term(8)
    diag = h.diagonal().reshape(3, 3)
    left, right = np.array(board.seams).T

    def energy(load):
        return np.sum(diag[load[left], load[right]]) + 8 * np.sum(
            (6 - board.degree) * load**2
        )

    filled = np.ones(board.n, dtype=int)
    costs = set()
    for pos, pair in pins.items():
        for pattern in ((0, 2), (2, 0)):
            ref = filled.copy()
            ref[list(pair)] = pattern
            costs.add(energy(ref) - energy(filled))
            if pos != 0:
                continue
            for n, (block, bits) in enumerate(m.sectors(8, None)):
                for j, word in enumerate(bits):
                    test = ref.copy()
                    test[list(rail)] -= word
                    field = m.mediator.contact(8)[0, pattern[0]] / 2
                    assert (
                        energy(test) - energy(ref) - m.MU * n
                        == block[j, j] + field * word[0]
                    )
    assert len(costs) == 1


def trace_log_series(n, z, a, b):
    """Exact unsplit matrix on the n active binary sites at fixed fugacity."""
    h = np.zeros((2**n, 2**n), dtype=object)
    for u in range(len(h)):
        bits = [(u >> j) & 1 for j in range(n)]
        h[u, u] = (
            Q(a) * bits[0]
            + Q(b) * bits[-1]
            + sum((j + 1) * x for j, x in enumerate(bits))
        )
        for j in range(n - 1):
            h[u, u] -= 17 * bits[j] * bits[j + 1]
            if bits[j] != bits[j + 1]:
                h[u, u ^ (3 << j)] = Q(-1, 2)
    weights = [z ** u.bit_count() / (1 + z) ** n for u in range(len(h))]
    power = np.eye(len(h), dtype=int).astype(object)
    partition = [Q(1)]
    log = [Q(0)]
    for k in range(1, n + 2):
        power = power @ (-h)
        partition.append(
            sum(w * x for w, x in zip(weights, power.diagonal(), strict=True))
            / factorial(k)
        )
        log.append(
            partition[k]
            - sum((j * log[j] * partition[k - j] for j in range(1, k)), Q(0)) / k
        )
    return log


def test_leading_mixed_response_at_three_fillings():
    for n, z in product((2, 3, 4), (Q(1, 3), Q(1), Q(2))):
        a, b = Q(-17), Q(31, 2)
        arms = [
            trace_log_series(n, z, u, v) for u, v in ((a, b), (0, 0), (a, 0), (0, b))
        ]
        mixed = [x + y - u - v for x, y, u, v in zip(*arms, strict=True)]
        assert mixed[:-1] == [0] * (n + 1)
        p = z / (1 + z)
        assert mixed[-1] == a * b * 17 ** (n - 1) * (p * (1 - p)) ** n


def test_population_odds_and_disconnected_control():
    saved = json.loads((m.HERE / "results.json").read_text())
    row = m.panel(8, 1)
    odds = np.array(row["empty"])
    connected = (
        np.log(odds[1:, 1:] * odds[0, 0] / (odds[1:, :1] * odds[:1, 1:])) / m.BETA
    )
    np.testing.assert_allclose(
        connected.ravel(), row["signed_interactions"], atol=8e-14, rtol=0
    )
    assert m.acquire(row, 915546) == saved["settings"][0]["acquisition"]
    assert max(abs(np.array(m.panel(8, 3, 1)["signed_interactions"]))) < 1e-12
    assert m.panel(8, 3, 1, 4)["quenched_mean"] > 0.0208
    bounds = json.loads((m.HERE / "bounds.json").read_text())
    assert bounds["second_difference"][0] > 0
