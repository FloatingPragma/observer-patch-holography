from itertools import combinations

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate as main_law
from model.evaluation import repair_band
from model.executable_model import cell, federation, laws
from . import candidate, replay


@pytest.mark.parametrize("level", ["cell", 0, 1, 2])
def test_one_load_is_laplacian_and_differs_only_by_the_main_law_energy_origin(level):
    board = cell() if level == "cell" else federation(level)
    hopping = candidate.generator(board)
    assert (hopping - board.laplacian() / 2).nnz == 0
    difference = main_law.generator(board) - hopping
    assert (difference - 6 * main_law.KAPPA * sparse.eye(board.n)).nnz == 0


def test_full_seam_restricts_to_exchange_and_cell_carries_the_gram_projector():
    exchange, _ = law.local_term_parts()
    for col, (a, b) in enumerate(law.SRC):
        image = np.zeros(law.DIM)
        image[col] += 0.5
        image[law.SRC.index((b, a))] -= 0.5
        np.testing.assert_allclose(exchange[:, col], image, rtol=0, atol=1e-15)
    band = repair_band.measure(candidate)
    assert band["one_load_scalar_offset"] == 0
    assert band["lifted_rank"] == 3
    assert laws.main_candidate() == main_law.NAME
    assert laws.load(candidate.NAME).ROLE == "control"


def test_different_two_load_records_from_independently_assembled_swaps():
    board = cell()
    pairs = list(combinations(range(board.n), 2))
    indices = {pair: i for i, pair in enumerate(pairs)}
    seam_indices = [indices[tuple(edge)] for edge in board.seams]
    exchange = np.zeros((len(pairs), len(pairs)))
    for a, b in board.seams:
        transposition = {a: b, b: a}
        for col, pair in enumerate(pairs):
            image = tuple(sorted(transposition.get(p, p) for p in pair))
            exchange[col, col] += 0.5
            exchange[indices[image], col] -= 0.5

    readouts = replay.experiment()
    for kappa in (0, 8):
        matrix = exchange.copy()
        matrix[seam_indices, seam_indices] -= 2 * kappa
        amplitudes = expm(-1j * matrix)[:, indices[(0, 1)]]
        probabilities = abs(amplitudes) ** 2
        separated = sum(p for pair, p in zip(pairs, probabilities) if pair not in board.seams)
        actual = readouts[f"kappa_{kappa}"]
        np.testing.assert_allclose(actual["by_seam"], probabilities[seam_indices], atol=1e-12, rtol=0)
        assert actual["no_pair"] == pytest.approx(separated, abs=1e-12)
    assert readouts["kappa_0"]["no_pair"] - readouts["kappa_8"]["no_pair"] > 0.45
    assert candidate.run(["--check"]) == 0
