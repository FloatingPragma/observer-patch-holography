"""Check the charge standards, the boundary spectrum and the finite-time enclosure."""

from fractions import Fraction
from types import SimpleNamespace

import numpy as np
from scipy.linalg import expm

from . import exterior


def test_opened_and_closed_boundary_against_dense_flight():
    # Two one-load dimers joined by a written connector. The rest is an exterior dimer.
    edges = [(0, 1), (1, 2), (2, 3), (0, 4), (3, 5), (4, 5)]
    board = SimpleNamespace(n=6, degree=[5, 6, 6, 5, 6, 5], seams=edges)
    for opened in (None, (0, 4), (3, 5)):
        terms = exterior.BoardTerms(board, [(0, 1), (2, 3)], [(1, 2)], opened)
        native = np.zeros((64, 64))
        for word in range(64):
            load = [1 - ((word >> p) & 1) for p in range(6)]
            native[word, word] = 8 * sum(load[i] for i in (0, 3, 5)) + 32 * (
                1 - load[1] * load[2]
            )
            for a, b in [(0, 1), (2, 3), (4, 5), *([opened] if opened else [])]:
                if load[a] != load[b]:
                    native[word, word] += 8.5
                    native[word ^ (1 << a) ^ (1 << b), word] -= 0.5
        for denom in exterior.DELAYS:
            flight = expm(-1j * native / denom)
            for initial in range(16):
                exact = float(np.sum(abs(flight[16:, initial]) ** 2))
                lower, upper = exterior.vacancy_bounds(
                    terms.powers(initial), terms.outside, denom
                )
                assert (
                    lower / exterior.RESOLUTION - 1e-15
                    <= exact
                    <= upper / exterior.RESOLUTION + 1e-15
                )
                if not opened:
                    assert exact == 0
                else:
                    diag, flips = terms.column(initial)
                    col = np.zeros(64)
                    col[initial] = diag
                    col[flips] = -1
                    ref = native[:, initial].copy()
                    ref[initial] -= native[0, 0]
                    np.testing.assert_array_equal(2 * ref, col)


def test_terminal_unit_standard_detects_polarization_and_aliases_ternary_charge():
    words = exterior.boundary_words()
    assert words["01"] == list(reversed(words["10"]))
    assert words["01"].count("blank") == 1
    assert words["11"].count("blank") == 2
    assert words["00"] == words["02"] == words["22"]
    exterior.check_contacts()
    assert (
        exterior.root_bounds(Fraction(2))[0] ** 2
        <= 2
        <= exterior.root_bounds(Fraction(2))[1] ** 2
    )
