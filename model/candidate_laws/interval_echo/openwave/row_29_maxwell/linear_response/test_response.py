import gzip
import json

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from . import measure, parked


def test_entangled_collar_receiver_under_local_message():
    h0, h1 = parked.integer_blocks()
    marked = np.arange(11) * 9 + law.SRC.index((1, 1))
    local = ((h0 + 8 * h1) / 4)[np.ix_(marked, marked)]
    rng = np.random.default_rng(481)
    amplitudes = rng.normal(size=(11, 11)) + 1j * rng.normal(size=(11, 11))
    amplitudes /= np.linalg.norm(amplitudes)
    flip, _ = np.linalg.qr(rng.normal(size=(11, 11)) + 1j * rng.normal(size=(11, 11)))
    u = expm(-1j * 0.37 * local)
    left = u @ amplitudes @ u.T
    right = u @ flip @ amplitudes @ u.T
    np.testing.assert_allclose(left.conj().T @ left, right.conj().T @ right, atol=1e-13)
    assert (
        parked.certificate()
        == json.loads((measure.HERE / "collar.json").read_text())["blocks"]
    )


def test_direct_triangle_and_edge_currents_for_arbitrary_perturbations():
    rng = np.random.default_rng(783)
    u = (rng.normal(size=12) + 1j * rng.normal(size=12)) / 15
    g = 0.08
    diff = (measure.ports(g + u) - measure.ports(g - u)) / 2
    odd = diff[:, 0] - diff[:, 1]
    a, b = measure.EDGE.T
    np.testing.assert_allclose(odd[:30] / 2, g * (u[b].imag - u[a].imag), atol=1e-16)
    np.testing.assert_allclose(odd[30:], 0, atol=1e-16)
    p = measure.ports(u)
    x, y, z = measure.FACE.T
    face = np.imag(u[x].conj() * u[y] + u[y].conj() * u[z] + u[z].conj() * u[x])
    np.testing.assert_allclose((p[30:, 0] - p[30:, 1]) * np.sqrt(3) / 2, face)


def test_complete_receipt_and_held_currents():
    arrays = {}
    for key in ("flight", "ground"):
        with gzip.open(measure.HERE / (key + ".npy.gz"), "rb") as stream:
            arrays[key] = np.load(stream, allow_pickle=False)
        assert (arrays[key].sum(axis=-1) == measure.SHOTS).all()
    report = measure.report(arrays)
    assert report == json.loads((measure.HERE / "results.json").read_text())
    band = report["raw_current_bound"]
    for row in report["panels"]:
        assert row["max_odd_face"] < band
        assert row["max_circular_kick_odd_edge"] < band
        if row["gain"] == 1:
            assert row["max_odd_edge"] > 2 * band and row["max_even_face"] > 2 * band
        if row["gain"] == 0:
            assert row["max_odd_edge"] < band and row["max_even_face"] < 2 * band
