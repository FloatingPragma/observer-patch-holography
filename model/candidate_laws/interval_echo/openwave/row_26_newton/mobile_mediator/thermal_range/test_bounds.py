"""Independent source, path and equilibrium checks on the range estimate."""

from fractions import Fraction
from itertools import product
import json

import numpy as np
from scipy.linalg import expm

from .. import matched_energy, measure
from . import bounds


def test_native_interpolation_bounds_and_critical_activity():
    board = measure.federation(0)
    result = bounds.run()
    assert json.loads((bounds.HERE / "results.json").read_text()) == result
    assert [(p["support_distance"], p["signs"][1]["connected_grand_potential"])
            for p in result["panels"]] == [
                (4, "negative"), (7, "negative"), (10, "positive"), (12, "positive")]
    assert all(p["signs"][0]["connected_grand_potential"] == "negative"
               for p in result["panels"])
    assert bounds.uniform_range(6) / bounds.uniform_range(5) == Fraction(6, 35)**2
    assert set(board.degree) == {5, 6}
    for distance, _, ends in measure.sources(board):
        if distance < 2:
            continue
        a = matched_energy.regular_pair(board, ends[0])
        b = matched_energy.regular_pair(board, ends[2])
        for q, reverse_a, reverse_b in product((-1, 1), (False, True), (False, True)):
            pins = (*(a[::-1] if reverse_a else a), *b)
            graph, h, v, w = measure.operators(board, pins, 8, q, reverse_b)
            assert len(graph) == bounds.MOBILE_PORTS
            assert graph.sum(axis=1).max() / 2 <= bounds.HOP_ROW_BOUND
            for field in (v, w):
                assert 0 < np.count_nonzero(field) <= bounds.SUPPORT_BOUND
                assert np.max(abs(field)) / 2 <= bounds.PIN_BOUND
            assert not np.any(v * w)
            assert h.diagonal().max() / 2 <= 133 / 2
            for i, j in measure.SETTINGS:
                assert min(h.diagonal() + i * v + j * w) / 2 >= bounds.DIAGONAL_FLOOR


def test_fair_path_average_and_grand_sign_against_direct_partition_functions():
    # Two allowed vertices each touch a single dipole constituent; two hops close a tour.
    h = np.array([[5.0, -0.5], [-0.5, 5.0]])
    for beta in (0.05, 0.2):
        pooled = np.zeros((2, 2))
        frozen = np.zeros_like(pooled)
        for first, second in product((-17, 15.5), repeat=2):
            for i, j in measure.SETTINGS:
                arm = h + np.diag([i * first, j * second])
                pooled[i, j] += np.trace(expm(-beta * arm)) / 4
                frozen[i, j] += np.exp(-beta * arm.diagonal()).sum() / 4
        z0 = pooled[0, 0]
        a = pooled[1, 0] - z0
        b = pooled[0, 1] - z0
        c = pooled[1, 1] - a - b - z0
        assert min(a, b, c) > 0
        assert abs(frozen[1, 1] + frozen[0, 0] - frozen[1, 0] - frozen[0, 1]) < 1e-13
        crossing = c / (a * b - z0 * c)
        for factor in (0.5, 2.0):
            z = factor * crossing
            direct = np.log1p(z * pooled)
            energy = -(direct[1, 1] + direct[0, 0] - direct[1, 0] - direct[0, 1]) / beta
            assert (energy < 0) == (factor < 1)
            assert max(0, -energy) <= z * c / beta
            assert max(0, energy) <= z**2 * a * b / (beta * (1 + z * z0)**2)


def test_vacancy_density_and_uniform_competition_exponent():
    for level in (0, 1, 2):
        board = measure.federation(level)
        m = 2**level
        assert board.n == 240 * m**2
        assert np.count_nonzero(board.degree == 5) == 9 * board.cells
        assert 2 * (9 * board.cells - 4) >= board.n - 4
        # A shortest cell path to either endpoint via home costs at most 12m.
        assert board.n - 4 >= (12 * m)**2
    gap = bounds.DIAGONAL_FLOOR - bounds.HOP_ROW_BOUND
    wait_energy = Fraction(69, 2)
    assert 2 * gap - Fraction(4, 5) * wait_energy == Fraction(7, 5)
