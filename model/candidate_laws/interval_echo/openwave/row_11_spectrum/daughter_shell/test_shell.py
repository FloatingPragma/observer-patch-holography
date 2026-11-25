from fractions import Fraction as Q
import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from . import measure


def test_every_daughter_row_and_microscopic_potentials():
    board = federation(0)
    d = board.degree
    adjacency = board.adjacency().astype(int)
    adjacent = adjacency.toarray()
    valid = ~np.eye(board.n, dtype=bool)
    assert valid.sum() == 57360
    weight = np.where(
        (d[:, None] == 6) & (d[None, :] == 5) & adjacent.astype(bool), 7, 1
    )
    np.fill_diagonal(weight, 0)
    neighbour_weights = (
        adjacency @ weight + (adjacency @ weight.T).T + adjacent * weight.T
    )
    difference = d[:, None] - d[None, :]
    potential = 8 * (12 + 2 * difference + 2 * adjacent)
    diagonal = 2 * potential + d[:, None] + d[None, :] - adjacent
    twice_upper = diagonal * weight + neighbour_weights
    assert np.all(7 * twice_upper[valid] <= 2 * 942 * weight[valid])
    edges = np.asarray(board.seams)
    exchange, contrast = candidate.native_law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    np.testing.assert_allclose(h[9:, :9], 0, atol=2e-13)
    for dh in (5, 6):
        for dx in (5, 6):
            for link in (0, 1):
                cases = np.argwhere(
                    valid & (d[:, None] == dh) & (d[None, :] == dx) & (adjacent == link)
                )
                if not len(cases):
                    continue
                u, x = cases[0]
                word = np.ones(board.n, dtype=int)
                word[u], word[x] = 0, 2
                # Contract the native 99-state seam column for every physical bond.
                a, b = word[edges].T
                local_indices = 3 * a + b
                value = (
                    int(candidate.KAPPA) * ((6 - d) @ (word**2 - 1))
                    + h.diagonal()[local_indices].sum()
                )
                np.testing.assert_allclose(2 * value, diagonal[u, x], atol=1e-11)
                assert len(set(measure.pair_moves(board, (u, x)))) == d[u] + d[x] - link


def test_unrestricted_enclosures_against_larger_sparse_release_and_calibration():
    saved = json.loads((measure.HERE / "results.json").read_text())
    b = federation(0)
    source, target = tuple(saved["source"]), tuple(saved["target"])
    states = measure.triplet.ball(b, (source,), measure.RADIUS + 1)
    i, j = states.index(source), states.index(target)
    vector = np.zeros(len(states))
    vector[i] = 1
    for name, k in (("native", 8), ("exchange_only_control", 0)):
        h = measure.triplet.restriction(b, states, k)
        h.setdiag(h.diagonal() - h[i, i])
        for trial in saved["release"][name]["records"]:
            t = float(Q(trial["exposure"]))
            actual = (
                abs(
                    expm_multiply(
                        -1j * t * h, vector, traceA=-1j * t * h.diagonal().sum()
                    )[j]
                )
                ** 2
            )
            lo, hi = np.array(trial["joint_probability_ticks"]) / measure.TICKS
            assert lo <= actual <= hi
            assert sum(trial["simulated_click_and_other"]) == measure.SHOTS
    assert saved["calibration"]["joint_calibration_dimension"] == 1320
    for name, expected in (
        ("native_energy_read", 15.5),
        ("exchange_only_energy_read", -0.5),
    ):
        cal = saved["release"][name]
        assert cal["exact_initial_slope"] == expected
        for trial in cal["records"]:
            t = float(Q(trial["exposure"]))
            # Independent sampling errors for the opposite-sign launch settings.
            bound = 5 / (t * np.sqrt(2 * measure.SHOTS))
            assert (
                abs(trial["slope_from_counts"] - trial["expected_finite_time_slope"])
                < bound
            )
            assert all(
                sum(row) == measure.SHOTS
                for row in trial["plus_minus_counts_p_d_other"]
            )
        offsets = [
            abs(row["expected_finite_time_slope"] - expected) for row in cal["records"]
        ]
        assert offsets[1] < 0.26 * offsets[0]
    assert (
        saved["release"]["required_receiver_charge"]
        != saved["release"]["sealed_initial_charge"]
    )
    assert json.loads(json.dumps(measure.run())) == saved
