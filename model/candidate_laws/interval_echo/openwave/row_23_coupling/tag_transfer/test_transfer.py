import json

import numpy as np
import pytest
from scipy.sparse import eye
from scipy.sparse.linalg import splu

from model.candidate_laws.interval_echo.openwave.row_23_coupling import experiment
from model.candidate_laws.interval_echo.openwave.row_23_coupling.tag_transfer import (
    measure,
)


@pytest.mark.parametrize("level", [0, 1])
def test_direct_board_response(level):
    board, seam, sources = measure.setup(level)
    profiles = sources[:4]
    blank = board.laplacian().tocsc() * 0.5
    endpoint = profiles[:2]
    with np.load(measure.HERE / "counts.npz") as bank:
        for stock, contrast in [("blank", 8), ("written", 8), ("written", 4)]:
            for amount in (1, 2):
                local = (
                    contrast * amount**2 * np.eye(2) + np.outer([1, -1], [1, -1]) / 2
                )
                h = blank.copy().tolil()
                if stock == "written":
                    h[np.ix_(seam, seam)] = h[np.ix_(seam, seam)].toarray() - local
                raw = bank[f"{board.cells}_{stock}_{contrast}_{amount}"]
                for z in [*measure.FREQUENCIES, -16 + 1j]:
                    solver = splu(z * eye(board.n, format="csc") - h.tocsc())
                    measured = measure.decode(raw, z)
                    exact = profiles @ solver.solve(profiles.T)
                    radii = (2 - np.eye(4)) * measure.precision(z)
                    assert np.all(abs(measured - exact) <= radii)
                    if stock == "written" and amount == 1:
                        # A second direct solve with double source tests the complete resolvent identity.
                        h[np.ix_(seam, seam)] -= 3 * contrast * np.eye(2)
                        double = splu(z * eye(board.n, format="csc") - h.tocsc())
                        target = profiles[2:] @ double.solve(profiles[2:].T)
                        h[np.ix_(seam, seam)] += 3 * contrast * np.eye(2)
                        pred, _ = measure.forecast(exact, -3 * contrast, 0, 0)
                        np.testing.assert_allclose(pred, target, atol=2e-13)
    np.testing.assert_array_equal(endpoint[:, seam], np.eye(2))


def test_cutoff_subtraction():
    board, seam, _ = measure.setup(0)
    energy, modes = np.linalg.eigh(board.laplacian().toarray() / 2)
    odd = np.outer([1, -1], [1, -1]) / 2
    single, double = [-8 * q * q * np.eye(2) - odd for q in (1, 2)]
    expected = np.linalg.inv(double) - np.linalg.inv(single)
    even_odd = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    np.testing.assert_allclose(
        even_odd.T @ expected @ even_odd, np.diag([3 / 32, 8 / 99]), atol=1e-14
    )
    for z in measure.FREQUENCIES:
        for cutoff in (0.75, 1.5, 2.25):
            keep = 2 * energy > cutoff**2
            vectors = modes[np.ix_(seam, np.flatnonzero(keep))]
            omitted = (vectors / (z - energy[keep])) @ vectors.T
            vertices = [experiment.dressed(v, omitted) for v in (single, double)]
            np.testing.assert_allclose(
                np.linalg.inv(vertices[1]) - np.linalg.inv(vertices[0]),
                expected,
                atol=6e-14,
            )


def test_held_decisions():
    with np.load(measure.HERE / "counts.npz") as bank:
        data = measure.summarize(bank)
    from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
        compare,
    )

    compare(
        data, json.loads((measure.HERE / "results.json").read_text()), "tag transfer"
    )
    fit = data["calibration"]
    assert abs(fit["value"] + 24) < fit["error"]
    assert fit["scalar_residual"] < 2 * fit["error"]
    for stock, k, excluded in [("blank", 8, 0), ("written", 8, 20), ("written", 4, 20)]:
        held = [r for r in data["held"] if (r["stock"], r["contrast"]) == (stock, k)]
        assert len(held) == 8
        assert all(np.all(np.array(r["residual"]) <= r["joint_radius"]) for r in held)
        assert sum(sum(r["copy_excluded"]) for r in held) == excluded
        if k == 4:
            assert sum(sum(r["full_contrast_excluded"]) for r in held) == 20
