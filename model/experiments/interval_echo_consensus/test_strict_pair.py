from fractions import Fraction

import numpy as np

from model.executable_model import federation
from . import strict_pair as pair


def test_native_repair_and_all_binary_readers():
    pair.verify_columns()
    board = federation(0)
    assert pair.reader_census(board) == {"internal": {8: 308, 9: 284, 10: 8}, "intercell": {10: 30}}
    sources = pair.strict_scheduler(board)
    assert sources[156] == [156, 157, 2]
    assert {sign for _, _, sign in sources} == {-2, 2}


def test_return_estimate_including_vacuous_regime():
    assert pair.return_floor(10, 8) == Fraction(717, 727) ** 2
    assert pair.return_floor(9, 8) == Fraction(84, 85) ** 2
    assert pair.return_floor(10, 6) == Fraction(333, 343) ** 2
    assert pair.return_floor(10, 2) == 0
    assert pair.return_floor(10, 0) == 0
    assert pair.return_floor(10, 8, odd=1) > pair.return_floor(10, 8)


def test_weak_and_strong_flights_through_both_readers():
    board = federation(0)
    tag = (156, 157)
    home = [edge for edge in board.seams if edge[0] // 12 == edge[1] // 12 == 13]
    rows = pair.strict_flight(board, tag, home, strengths=(0, 8), times=(0, 1))
    reference = {0: (0.09966404598, 0.32320319077), 8: (0.99456025639, 0.00331920038)}
    marked = home.index(tag)
    for row in rows:
        direct, flags = (np.asarray(row[k]) for k in ("six_probe", "pair_flag"))
        assert direct.shape == flags.shape == (31,)
        assert flags[marked] == 0
        np.testing.assert_allclose(np.delete(direct[:-1], marked), np.delete(flags[:-1], marked))
        np.testing.assert_allclose(direct.sum(), 1, atol=1e-12, rtol=0)
        np.testing.assert_allclose(flags[-1], direct[-1] + direct[marked], atol=1e-12, rtol=0)
        expected = (1, 0) if row["t"] == 0 else reference[row["kappa"]]
        np.testing.assert_allclose((direct[marked], flags[:-1].sum()), expected, atol=2e-10, rtol=0)
