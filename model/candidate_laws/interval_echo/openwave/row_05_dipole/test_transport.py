"""Independent spectral and count checks for mobile signed-current readings."""

import numpy as np

from model.executable_model import federation
from .experiment import gates_of
from .test_dipole import read
from .transport import modes


def test_native_dense_spectrum_predicts_receiver_orbitals():
    data = read("transport")
    shots = data["shots"]
    board = federation(0)
    for run in (r for r in data["boards"] if r["owners"] == 20):
        k = run["kappa"]
        h = board.laplacian().toarray() / 2 + np.diag(k * (2 * board.degree - 6))
        energy, eigenvectors = np.linalg.eigh(h)
        _, source, _ = modes(gates_of(board, 0))
        coeffs = eigenvectors[:12].T @ np.column_stack([source, source.conj()])
        for panel in run["reads"]:
            o = panel["path"][-1]
            block = eigenvectors[12 * o : 12 * (o + 1)]
            states = block @ (np.exp(-1j * energy * panel["time"])[:, None] * coeffs)
            compact, plus, _ = modes(gates_of(board, o))
            assert np.max(abs(compact.T @ states)) < 1e-10
            for j, state in enumerate(states.T):
                positive = abs(np.vdot(plus, state)) ** 2 if plus is not None else 0
                negative = abs(plus @ state) ** 2 if plus is not None else 0
                counts = panel["sources"][j]["mode_counts"]
                np.testing.assert_allclose(
                    np.array(counts[:2]) / shots, [positive, negative], atol=0.0002
                )


def test_receiver_counts_dark_modes_and_phase_kick():
    data = read("transport")
    shots = data["shots"]
    for board in data["boards"]:
        for row in board["reads"]:
            bright, _, dark, conjugate, mixed = row["sources"]
            for source in row["sources"]:
                assert sum(source["mode_counts"]) == shots
                assert source["mode_counts"][2] == 0
                assert np.all(np.sum(source["current_counts"], axis=1) == shots)
            assert (
                dark["mode_counts"] == conjugate["mode_counts"] == [0, 0, 0, 0, shots]
            )
            np.testing.assert_allclose(
                np.array(mixed["mode_counts"][:4]) / shots,
                np.array(bright["mode_counts"][:4]) / (2 * shots),
                atol=0.0002,
            )
    run = next(b for b in data["boards"] if b["owners"] == 320 and b["kappa"] == 8)
    first = next(r for r in run["reads"] if r["time"] == 128 and len(r["path"]) == 2)
    pos, neg = first["sources"][:2]
    assert (
        pos["moment"] > data["moment_radius"] and neg["moment"] < -data["moment_radius"]
    )
    second = next(r for r in run["reads"] if r["time"] == 64 and len(r["path"]) == 3)
    firsttwo = np.array([s["mode_counts"][:2] for s in second["sources"][:2]])
    sign = firsttwo @ [1, -1] / shots
    assert sign[0] > data["mode_score_radius"] and sign[1] < -data["mode_score_radius"]
    bright, compact = run["kicks"]
    assert np.allclose(bright["matrix"], [-1 / 12, 0])
    assert abs(np.ptp(bright["mean_energies"]) - 1 / 6) < 1e-9
    assert min(bright["widths"]) > 8 and np.max(abs(np.array(compact["widths"]))) == 0
    for kick in bright["reads"]:
        xy = np.array(kick["ramsey_xy"]) @ [1, -1] / shots
        assert abs(-np.angle(complex(*xy)) / kick["time"] - kick["energy_read"]) < 1e-9
