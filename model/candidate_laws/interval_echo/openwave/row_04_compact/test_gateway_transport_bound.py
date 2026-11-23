from fractions import Fraction

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from model.executable_model.carrier import Carrier
from . import gateway_transport_bound as bound


def test_walk_and_poisson_bound():
    walk, early, ceiling = bound.swap_home_ceiling()
    assert walk == Fraction(327415820588960477, 592297667290202112)
    assert 0 < early < Fraction(15, 1000)
    assert walk < ceiling < Fraction(559219, 10**6)

    # An off-protocol exposure compares the bound to the continuous-time law
    # directly, including initial ports outside this pilot's home.
    board = federation(0)
    in_home = np.repeat(np.arange(board.cells) == 13, 12).astype(float)
    observed = expm_multiply(-0.5 * bound.EXPOSURE * board.laplacian(), in_home)
    assert 0 < np.min(observed) <= np.max(observed) < float(ceiling)


def test_isolating_all_cells_preserves_every_home(monkeypatch):
    board = federation(0)
    isolated = Carrier("isolated control", board.n,
                       [(a, b) for a, b in board.seams if a // 12 == b // 12])
    monkeypatch.setattr(bound, "federation", lambda level: isolated)
    walk, _, ceiling = bound.swap_home_ceiling()
    assert walk == ceiling == 1


def test_nonstochastic_uniformization_rejected():
    board = federation(0)
    # A seventh seam at an occupied gateway makes its proposed holding weight negative.
    gateway = int(np.flatnonzero(board.degree == 6)[0])
    other = next(p for p in range(board.n) if p != gateway and p not in board.neighbours[gateway])
    changed = Carrier("extra seam", board.n, [*board.seams, (gateway, other)])
    with pytest.raises(ArithmeticError, match="stochastic"):
        bound._walk_ceiling(changed)
