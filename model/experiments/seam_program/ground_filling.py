"""Check the zero-temperature grand-canonical ground space on the blank binary federation."""

from itertools import combinations

import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.sparse.csgraph import maximum_flow

from model.candidate_laws.contrast_memory_spherical_v1 import candidate as reference
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation


def local_bound():
    """Return the worst partial-cell value of 4(C - 3N/4), counting internal seams."""
    car = cell()
    a, b = np.array(car.seams).T
    words = np.arange(1, 2**car.n - 1)[:, None]
    occupations = (words >> np.arange(car.n)) & 1
    cuts = np.sum(occupations[:, a] != occupations[:, b], axis=1)
    number = occupations.sum(axis=1)
    worst = 4 * len(car.seams)
    for gates in combinations(range(car.n), 3):
        references = number - occupations[:, gates].sum(axis=1)
        scaled_excess = 4 * (cuts + references) - 3 * number
        worst = min(worst, int(scaled_excess.min()))
    return int(cuts.min()), worst


def check_boards():
    # These source checks bind the diagonal cut energy and positive exchange to the native law.
    exchange, contrast = candidate.native_law.local_term_parts()
    ids = [candidate.native_law.SRC.index(pair) for pair in ((0, 0), (0, 1), (1, 0), (1, 1))]
    cols = np.eye(len(exchange))[:, ids]
    odd = np.array([0., 1., -1., 0.])
    np.testing.assert_allclose(exchange @ cols, cols @ np.outer(odd, odd) / 2, atol=1e-14)
    np.testing.assert_allclose(contrast @ cols, cols @ np.diag(odd**2), atol=1e-14)
    cell_cuts = []
    for level in range(3):
        board = federation(level)
        degree = board.degree.reshape(-1, 12)
        assert np.all(np.sort(degree, axis=1) == [5] * 9 + [6] * 3)
        internal = [set() for _ in range(board.cells)]
        for a, b in board.seams:
            if a // 12 == b // 12:
                internal[a // 12].add((a % 12, b % 12))
        pattern = set(cell().seams)
        assert all(block == pattern for block in internal)
        a, b = np.array(board.external).T // 12
        directed = coo_matrix((np.ones(2 * len(a), dtype=int), (np.r_[a, b], np.r_[b, a])),
                              shape=(board.cells, board.cells)).tocsr()
        # Any cut separates the root from some other cell. All capacities count actual bridges.
        cell_cuts.append(int(min(maximum_flow(directed, 0, j).flow_value for j in range(1, board.cells))))
        difference = candidate.generator(board) - board.laplacian()/2 - 6*candidate.KAPPA*eye(board.n)
        np.testing.assert_allclose(difference.data, 0., atol=1e-13)
        # With contrast removed and chemical potential zero, a third exact ground vector exists.
        np.testing.assert_allclose(reference.generator(board, 0) @ np.ones(board.n), 0., atol=1e-13)
    return cell_cuts


if __name__ == "__main__":
    cell_cuts = check_boards()
    cut, excess = local_bound()
    assert cut == 5 and excess == 17 and cell_cuts == [3, 3, 3]
    print(f"Minimum partial-cell cut {cut}, minimum partial-cell excess {excess}/4.")
    print(f"Intercell minimum cuts {cell_cuts}.")
    print("Gap between 24 and 25.5 on these boards at coexistence, in native energy units.")
