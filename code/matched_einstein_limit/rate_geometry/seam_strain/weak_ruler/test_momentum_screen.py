import json
from itertools import combinations

import numpy as np
from scipy import sparse
import sympy as sp

from . import momentum_screen as screen
from .. import live_source


def test_junction_sign_and_shared_tangent_expansion():
    coords = {1: (0, 0, 0), 38: (1, 0, 0), 5: (0, 1, 0), 0: (0, 0, -1), 39: (0, 0, 1)}
    edges = list(combinations(screen.STATIONS, 2))
    lengths = np.asarray(
        [np.linalg.norm(np.subtract(coords[a], coords[b])) for a, b in edges]
    )
    trace, shear = screen.face_momentum(lengths, 0.25 * lengths)
    np.testing.assert_allclose([trace, *shear], [0.5, 0, 0], atol=1e-14)
    vel = np.zeros(len(edges))
    k = [set(e) for e in edges].index({38, 0})
    vel[k] = -1 / lengths[k]
    trace, shear = screen.face_momentum(lengths, vel)
    np.testing.assert_allclose([trace, *shear], [0, 0.5, 0], atol=1e-14)


def test_integer_gram_time_derivative_against_heisenberg_evolution():
    rng = np.random.default_rng(15)
    hi = rng.integers(-800, 800, size=(16, 16))
    lo = rng.integers(0, 1 << 36, size=hi.shape)
    h = rng.integers(-3, 3, size=(16, 16))
    h = h + h.T
    head = (
        sp.Matrix([[2, sp.I, 0, 0], [-sp.I, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]) / 6
    )
    actual = (
        live_source.rational_read(
            screen.derivative_grams(sparse.csr_matrix(h), hi, lo), head, (1, 2)
        )
        / 40
    )
    w = ((hi * float(1 << 36)) + lo) / float(1 << 82)
    v = w[:, :8] + 1j * w[:, 8:]
    dv = -1j * (h @ v) / 40

    def group(w):
        return w.reshape(2, 2, 4, 8).transpose(1, 0, 2, 3).reshape(2, 8, 8)

    sigma = np.kron(np.diag([1, 2]) / 3, np.array(head, complex))
    cross = sum(
        a @ sigma @ b.conj().T for a, b in zip(group(dv), group(v), strict=True)
    )
    np.testing.assert_allclose(
        np.array(actual, complex) * 1e23, (cross + cross.conj().T) * 1e23, atol=1e-14
    )


def test_physical_radar_face_and_independent_cut_contrasts():
    data = screen.run()
    saved = json.loads((screen.hinge.HERE / "momentum_screen.json").read_text())
    for record in (data, saved):
        assert list(record["face"]) == [1, 38, 5] and list(
            record["oriented_poles"]
        ) == [0, 39]
        assert record["material"]["ports"] == [230, 490, 240]
        assert record["material"]["boundary_bonds"] == [42, 66]
        for geom in record["geometry"].values():
            assert geom["normal_face_coefficient"] == [0, 0]
            assert len(geom["tangential_face_coefficients"]) == 2
        ref = record["geometry"]["reference"]
        assert ref["common_tangential_expansion"][1] < -0.03
        assert ref["tangential_face_coefficients"][0][0] > 0.0044
        assert ref["tangential_face_coefficients"][1][0] > 0.062
        for arm in ("ordinary", "decorrelated"):
            change = record["geometry"][arm]["tangential_minus_reference"]
            assert change[0][0] > 4e-6 and change[1][1] < -0.002
        cuts = record["material"]["oriented_cuts"]
        for cut in cuts.values():
            det = cut["energy_number_contrast_determinant"]
            assert det[1] < 0 or 0 < det[0]
            assert cut["any_energy_origin_max_gap_lower"] > 4e-5
        matched = cuts["half_shared_zone"]["matched_number_lottery"]
        assert (
            0
            < matched["ordinary_probability"][0]
            <= matched["ordinary_probability"][1]
            < 0.001
        )
        assert matched["number_minus_reference"] == [0, 0]
        assert matched["energy_minus_reference"][1] < -4e-5

    def overlap(a, b):
        assert max(a[0], b[0]) <= min(a[1], b[1])

    for key, panel in data["geometry"].items():
        archived = saved["geometry"][key]
        overlap(
            panel["common_tangential_expansion"],
            archived["common_tangential_expansion"],
        )
        for a, b in zip(
            panel["tangential_face_coefficients"],
            archived["tangential_face_coefficients"],
            strict=True,
        ):
            overlap(a, b)
    for name, cut in data["material"]["oriented_cuts"].items():
        archived = saved["material"]["oriented_cuts"][name]
        for measure in ("energy", "number"):
            for arm, reading in cut[measure].items():
                overlap(reading, archived[measure][arm])
