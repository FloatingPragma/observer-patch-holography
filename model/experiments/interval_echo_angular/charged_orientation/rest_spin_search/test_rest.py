import json

import numpy as np
from scipy.linalg import expm

from . import measure as m, release as r


def test_geometry_and_exact_dark_modes_on_connected_graph():
    b = m.fixture()
    c = m.fixture(True)
    assert b["open_bridge"] == c["open_bridge"] == [1, 82]
    assert c["shutter_count"] == b["shutter_count"] + 4
    independent = np.zeros((12, 12))
    for i, j in m.EDGES:
        independent[i, i] += 8.5
        independent[j, j] += 8.5
        independent[i, j] = independent[j, i] = -0.5
    np.testing.assert_allclose(independent, m.native()[0], atol=1e-13)
    for n in (1, 2):
        v = np.column_stack([m.mode(n), m.mode(-n)])
        ev = 25.5 - np.cos(2 * np.pi * n / 5)
        np.testing.assert_allclose(independent @ v, ev * v, atol=1e-12)
        np.testing.assert_allclose(v.conj().T @ v, np.eye(2), atol=1e-12)
        np.testing.assert_allclose(
            m.rotation() @ v,
            v @ np.diag(np.exp([-2j * np.pi * n / 5, 2j * np.pi * n / 5])),
            atol=1e-12,
        )
    saved = json.loads((m.HERE / "results.json").read_text())
    for row in saved["stationary"]:
        lo, hi = row["cost_meter"]["simultaneous99"]
        assert lo < row["formation"] < hi
        assert sum(row["mirror_plus_minus"]) == m.SHOTS
        assert row["interferometer"][3]["return_hits"] == m.SHOTS
    assert saved["supplied_spinor_reference"]["5"] == 0


def test_recapture_misses_sampled_witness_and_independent_free_propagation():
    saved = json.loads((m.HERE / "release.json").read_text())
    h = r.hamiltonian()
    assert (
        r.exact_channel(h)
        == saved["exact_determinant"]
        == dict(prime=31, fifth_root=2, order=8, nonzero_residue=20)
    )
    eigenvectors = np.column_stack(r.bases(0)[0])
    analyzer = np.column_stack(r.bases(1)[0]).conj().T
    y = np.linalg.matrix_power(expm(-0.5j * h), 4096) @ eigenvectors
    transfer = analyzer @ y
    row = next(a for a in saved["records"] if a["time"] == 2048)
    np.testing.assert_allclose(
        row["squared_singular_values"],
        np.linalg.svd(transfer, compute_uv=False) ** 2,
        atol=3e-10,
    )
    for b in row["reads"]:
        assert np.all(np.asarray(b["recapture_counts"]).sum(axis=1) == r.SHOTS)
        assert sum(b["current_bond"]) == r.SHOTS
    witness = row["sampled_witness"]
    exact = -abs(np.linalg.det(transfer)) / 2
    assert abs(witness["min_pt"] - exact) < witness["simultaneous99_radius"]
    assert witness["min_pt"] + witness["simultaneous99_radius"] < 0
    turns = saved["after_flight_turns"]["2048"]
    assert (
        abs(
            turns[0]["source_return_probability"]
            - turns[2]["source_return_probability"]
        )
        < 1e-11
    )
    assert (
        turns[-1]["source_return_probability"] - turns[2]["source_return_probability"]
        > 0.12
    )
