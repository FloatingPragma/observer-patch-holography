"""Independent direct propagation and a sign-blind receiver identity."""

import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import experiment as dipole
from . import measure


def test_complete_cube_and_fixed_analyzer():
    states = [np.outer(w, w.conj()) for w in measure.INPUTS]
    swap = np.eye(4)[[0, 2, 1, 3]]
    np.testing.assert_allclose(
        sum(np.kron(r, r) for r in states), np.eye(4) + swap, atol=1e-12
    )
    # A fixed nonzero state presented to any of the three measurement settings
    # gives the same total probability of landing in the plane.
    rng = np.random.default_rng(1721)
    cell = (rng.normal(size=(12, 2)) + 1j * rng.normal(size=(12, 2))) / 15
    f = dipole.frame(3)
    source = np.array([0.6, 0.8j])
    total = []
    for axis in measure.INPUTS.reshape(3, 2, 2):
        a = measure.detector(cell, np.eye(2), f, source, axis, False)
        b = measure.detector(cell, np.eye(2), f, source, axis[::-1], False)
        np.testing.assert_allclose([a[1], a[0], *a[2:]], b)
        total.append(sum(a[:2]))
    np.testing.assert_allclose(total, [total[0]] * 3)


def test_third_owner_unconditional_score_from_direct_propagation():
    board = dipole.carrier.federation(1)
    route, ports = dipole.path_and_gates(board)
    y = np.zeros((board.n, 2))
    y[:12] = dipole.frame(ports[0])[:, :2]
    at = expm_multiply(-32j * board.laplacian(), y)
    receiver = slice(12 * route[-1], 12 * route[-1] + 12)
    transfer = dipole.frame(ports[-1])[:, :2].T @ at[receiver]
    with np.load(dipole.HERE / "counts.npz") as pilot:
        quadratures = pilot["level1_full_unit_8.0_d3_quadratures"][7, :, :2, :2]
    polar = 2 * quadratures / dipole.SHOTS - 1
    u, _, vh = np.linalg.svd(polar[0] + 1j * polar[1])
    corrected = vh.conj().T @ u.conj().T @ transfer
    # Six-state score = 3*C-2*P(accept), including every failed start.
    expected_score = (
        abs(np.trace(corrected)) ** 2 - np.linalg.norm(corrected, "fro") ** 2
    ) / 2
    data = json.loads((measure.HERE / "records.json").read_text())
    lo, hi = data["panels"][2]["readings"][7]["arms"]["native"]["score_interval"]
    assert 0 < lo < expected_score < hi
