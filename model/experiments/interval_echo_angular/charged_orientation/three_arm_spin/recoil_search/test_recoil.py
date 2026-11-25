"""The original doublet recoils, but its frozen formation pole fails both length holdouts."""

import json

import numpy as np
import pytest
from scipy.linalg import expm

from ..apparatus import blocks
from .impulse import certify
from .measure import HERE, STARTS, run


@pytest.fixture(scope="module")
def replay():
    return run()


def test_local_analyzers():
    # Native half and full pulses and the endpoint phase generate the two
    # distinct instruments used on the same guide preparation.
    for key, hop, plus_energy, minus_energy in (
        ("blank", 0.5, 8.0, 9.0),
        ((1, 1), 4.0, -8.0, 0.0),
    ):
        matrix = blocks()[key] / 8
        i, j = next(
            (i, j)
            for i, j in zip(*np.triu_indices(len(matrix), 1), strict=True)
            if matrix[i, j] == -hop
        )
        pair = [i, j]
        restricted = matrix[np.ix_(pair, pair)]
        assert (
            np.count_nonzero(
                matrix[np.ix_(list(set(range(len(matrix))) - set(pair)), pair)]
            )
            == 0
        )
        H = restricted - np.eye(2) * matrix[4, 4]
        np.testing.assert_allclose(
            H @ np.array([1.0, 1.0]), [plus_energy] * 2, atol=1e-12
        )
        np.testing.assert_allclose(
            H @ np.array([1.0, -1.0]), [minus_energy, -minus_energy], atol=1e-12
        )
        half = expm(-1j * np.pi * H / (4 * hop))
        full = half @ half
        np.testing.assert_allclose(abs(full), [[0, 1], [1, 0]], atol=1e-12)
        axes = half @ np.diag([1, -1j]) @ np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        np.testing.assert_allclose(abs(axes), np.eye(2), atol=1e-12)
        probe = np.array([0.6, 0.8j])
        probs = abs(half @ probe) ** 2
        assert (
            abs(probs[1] - probs[0] - 2 * np.imag(probe[0].conj() * probe[1])) < 1e-13
        )
    for angle in (0.0, -0.25, np.pi, 0.37):
        assert certify(angle)["two_column_residual"] < 1e-12


def check_records(r):
    null = r["kicks"][2]
    for row in r["lengths"]:
        assert (
            sum(row[k] for k in ("accepted", "source_refused", "clock_refused"))
            == STARTS
        )
        lo, hi = row["probability_interval"]
        assert lo < row["source_probability"] < hi
    for k in r["kicks"]:
        assert (
            sum(k[j] for j in ("repeat_accept", "repeat_reject", "unprepared"))
            == STARTS
        )
        n = STARTS - k["unprepared"]
        assert abs(k["repeat_accept"] / n - k["conditional_repeat_filter"]) < np.sqrt(
            12 / n
        )
        for key, bins in [
            ("energy", "energy_bins"),
            ("measured_current", "current_bins"),
        ]:
            bank = k[key]
            assert abs(bank["estimate"] - bank["exact"]) < bank["radius"]
            assert all(
                sum(b[bins]) + b["source_misses"] == STARTS for b in bank["banks"]
            )
        assert all(
            sum(p["port_ordinal_counts"]) + p["misses"] == STARTS for p in k["position"]
        )
        p0, p1 = k["position"][0], k["position"][3]
        delta = p1["estimate"] - p0["estimate"]
        current = k["measured_current"]
        if k["k"]:
            assert np.sign(k["k"]) * delta > p0["radius"] + p1["radius"]
            assert np.sign(k["k"]) * current["estimate"] > current["radius"]
            assert (
                k["energy"]["estimate"] - null["energy"]["estimate"]
                > k["energy"]["radius"] + null["energy"]["radius"]
            )
        else:
            assert abs(delta) < p0["radius"] + p1["radius"]
            assert abs(current["estimate"]) < current["radius"]


def test_stored_and_fresh_readers(replay):
    saved = json.loads((HERE / "readings.json").read_text())
    for data in (saved, replay):
        check_records(data)
        original, first, second = data["lengths"]
        assert [r["arm_steps"] for r in data["lengths"]] == [19, 20, 21]
        assert original["source_probability"] > 6e-5
        assert (
            first["source_probability"] < 4e-13 and second["source_probability"] < 6e-12
        )
        assert all(
            row["distance"] - row["exclusion_allowance"] > 1e-3
            for row in [first, second]
        )
        assert data["projected_current"] < 1e-12 and data["pinned_control"] < 1e-12
        assert data["ideal_inverse_longitudinal_curvature"] == pytest.approx(
            5.69136, abs=1e-5
        )
    for stored, fresh in zip(saved["kicks"], replay["kicks"], strict=True):
        for key in ["velocity", "ideal_line_survival", "conditional_repeat_filter"]:
            np.testing.assert_allclose(stored[key], fresh[key], atol=3e-9)
    for stored, fresh in zip(saved["lengths"], replay["lengths"], strict=True):
        np.testing.assert_allclose(
            stored["source_probability"], fresh["source_probability"], rtol=1e-6
        )
        assert stored["added_ports"] == fresh["added_ports"]
