"""Separate dense-spectrum and receiver-only checks of the internal oscillation."""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.executable_model import federation
from .readout import decode


HERE = Path(__file__).resolve().parent


def test_stationary_dimension_and_gap_against_full_port_matrix():
    board = federation(0)
    h = board.laplacian().toarray() / 2
    eig, vectors = eigh(h)
    at_rest = abs(eig - 3) < 1e-9
    panel = json.loads((HERE / "results.json").read_text())["panels"][0]
    assert at_rest.sum() == panel["stationary_dimension"] == 70
    # Rebuild the single-cell source from port-to-port distance, without
    # the production projector or any cell-to-cell continuation algorithm.
    from scipy.sparse.csgraph import shortest_path

    a = board.adjacency()[:12, :12]
    opposites = shortest_path(a, directed=False, unweighted=True) == 3
    gateway = board.degree[:12] == 6
    positive = gateway | (opposites @ gateway)
    shaped = np.where(positive, 1, -sum(positive) / sum(~positive))
    measured = np.linalg.norm(vectors[:12, at_rest].T @ shaped) ** 2 / sum(shaped**2)
    np.testing.assert_allclose(measured, panel["local_source_stationary"], atol=6e-11)


def test_counted_rest_phase_changed_law_and_conservation():
    result = json.loads((HERE / "results.json").read_text())
    for panel in result["panels"]:
        for rate in (1, 2):
            center, (lo, hi) = decode(panel, rate)
            assert lo < 3 * rate < hi and abs(center - 3 * rate) < 0.04
        assert panel["spatial"][0]["stationary_weight"] == 1
        for r in panel["interferometry"]:
            assert sum(r["counts"]) == result["shots"]
        for mode in panel["spatial"][1:]:
            assert 0 < mode["stationary_weight"] < 1 and mode["odd"] > 0
            assert all(cs[1] + cs[3] > 0 for _, cs in mode["visits"])
        for t, cs in panel["occupation_by_hops"]:
            assert sum(cs) == result["shots"]
            if t == 0:
                assert sum(cs[1:]) == 0
        assert sum(panel["occupation_by_hops"][-1][1][2:]) > 100_000


def test_complement_clock_uses_all_outcomes():
    for board in json.loads((HERE / "complement.json").read_text()):
        for row in board["rows"]:
            hist = np.array(row["clock"])
            assert sum(hist) == 2**23
            assert np.max(hist + np.roll(hist, -1)) / sum(hist) < 0.79
