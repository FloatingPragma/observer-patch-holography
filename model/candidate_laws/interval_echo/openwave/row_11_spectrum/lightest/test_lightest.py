"""Microscopic sector restriction and integer certificates on the native boards."""

import json
from fractions import Fraction

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from . import measure


def test_full_native_columns_and_remote_positive_weights():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    block = exchange + 8 * contrast
    filled = law.SRC.index((1, 1))
    np.testing.assert_allclose(block[:, filled], 0, atol=2e-13)
    for incoming, outgoing in [((1, 0), (0, 1)), ((0, 1), (1, 0))]:
        sharp = law.SRC.index(incoming)
        partner = law.SRC.index(outgoing)
        expected = np.zeros(len(block))
        expected[[sharp, partner]] = [8.5, -0.5]
        np.testing.assert_allclose(block[:, sharp], expected, atol=2e-13)

    # Assemble one independent native vacancy Hamiltonian from its incidence
    # terms and its six-minus-degree private reference relative to a filled host.
    board = federation(0)
    h = np.diag(-8.0 * (6 - board.degree))
    for i, j in board.seams:
        h[i, i] += block[sharp, sharp]
        h[j, j] += block[sharp, sharp]
        h[i, j] += block[sharp, partner]
        h[j, i] += block[sharp, partner]
    probe = list(range(board.n))
    np.testing.assert_allclose(
        2 * h @ probe, measure.double_h_times(board, probe), atol=1e-10
    )
    energies, states = eigh(h, subset_by_index=(0, 1))
    assert energies[1] > energies[0] + 1e-7
    positive = states[:, 0] * np.sign(states[0, 0])
    assert np.min(positive) > 0
    row = json.loads((measure.HERE / "results.json").read_text())["boards"][0]
    assert Fraction(row["exact_componentwise_lower"]) < energies[0]
    assert energies[0] <= Fraction(row["exact_rayleigh_upper"]) + 5e-13


def test_certificates_all_boards_and_offset_invariance():
    stored = json.loads((measure.HERE / "results.json").read_text())
    assert measure.run() == stored
    for row in stored["boards"]:
        lower = Fraction(row["exact_componentwise_lower"])
        upper = Fraction(row["exact_rayleigh_upper"])
        assert 32 < lower <= upper < Fraction(69, 2) < 35
        lo, hi = row["compact_minus_ground_interval_ticks"]
        for offset in [Fraction(0), Fraction(73, 7)]:
            gap = 35 + offset - (upper + offset)
            assert Fraction(lo, measure.TICKS) <= gap <= Fraction(hi, measure.TICKS)
