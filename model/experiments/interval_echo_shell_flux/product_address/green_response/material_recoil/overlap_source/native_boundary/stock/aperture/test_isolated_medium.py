import json

import numpy as np
from scipy import sparse

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    rest_limit,
    isolated_medium as experiment,
)


def test_certificate_encloses_exact_two_frequency_uniform_launch():
    h = sparse.csr_matrix([[0.0, 0.0], [0.0, 2.0]])
    # The blank reads the uniform wave. Late dephasing removes its coherence.
    obs = np.array([[1.0, 3.0], [3.0, -1.0]])
    bank = rest_limit.eigen_witness(h)
    bounds = rest_limit.limits(obs, bank, 2)
    assert len(bank) == 2
    assert bounds["limit_lower"] <= -3 <= bounds["limit_upper"] < -2.99
    for t in (0.05, 3.0, 80.0, bounds["sign_guaranteed_from_age"]):
        actual = 3 * (np.sinc(2 * t / np.pi) - 1)
        width = bounds["finite_window_error_numerator"] / t
        assert bounds["limit_lower"] - width <= actual <= bounds["limit_upper"] + width
    degenerate = rest_limit.eigen_witness(sparse.eye(3, format="csr"))
    assert degenerate == []


def test_pairs_match_single_site_and_nearest_neighbor_preparation_reads():
    board = experiment.lab.cold.federation(0)
    receipt = json.loads((experiment.HERE / "isolated_medium.json").read_bytes())
    differences = np.zeros((board.n, board.n))
    for a, b, phase in receipt["selected_pairs"]:
        differences[a, b] = differences[b, a] = phase
    assert np.array_equal(abs(differences).sum(axis=1), np.ones(board.n))
    assert not np.any(differences.diagonal())
    assert all(differences[a, b] == 0 for a, b in board.seams)
    row = receipt["decoded"]
    assert row["paired_plus"]["65536"][0][0] > 0 > row["paired_minus"]["65536"][0][1]
    assert row["rest"]["65536"][0][0] > 0 > row["rest"]["65536"][1][1]
    assert all(lo <= 0 <= hi for bank in row["random_port"].values() for lo, hi in bank)
