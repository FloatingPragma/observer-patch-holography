import json

import numpy as np
import pytest

from . import instrument, measure


@pytest.mark.parametrize("load", [1, 2])
def test_native_two_readers_agree_on_every_density(load):
    h, marks, sorter = instrument.instruments(load)
    above_ground = h - np.linalg.eigvalsh(h)[0] * np.eye(4)
    np.testing.assert_allclose(marks.sum(axis=0), above_ground, atol=2e-13)
    np.testing.assert_allclose(
        sorter.conj().T @ np.diag(instrument.BITS) @ sorter, above_ground, atol=2e-13
    )
    for flips in (1, 2):
        np.testing.assert_allclose(
            instrument.beat_probability(load, flips, 0.6),
            (1 + np.cos(0.6 * flips)) / 2,
            atol=1e-13,
        )


def test_archived_counts_and_held_source_calibration():
    d = json.loads((measure.HERE / "results.json").read_text())
    n, w = d["starts"], d["source_probability_width"]
    for s in d["spectra"]:
        assert s["interval"][0] < s["flips"] < s["interval"][1]
    cal = d["calibrated_clicks_per_excitation"]
    lo, hi = (cal - w) / (1 + 2 * w), (cal + w) / (1 - 2 * w)
    for s in d["sources"]:
        core = s["full"]
        assert sum(core["edge_trials"]) == sum(s["energy_ports"]) == n
        assert all(
            0 <= k <= m
            for k, m in zip(core["edge_clicks"], core["edge_trials"], strict=True)
        )
        observed = sum(core["edge_clicks"]) / n
        energy = np.dot(s["energy_ports"], instrument.BITS) / n
        assert lo * (energy - 2 * w) - w <= observed <= hi * (energy + 2 * w) + w
        if "missing_first_edge" in s:
            ctrl = s["missing_first_edge"]
            assert sum(ctrl["edge_trials"]) == n and ctrl["edge_clicks"][0] == 0
            target = 1 / 8 if s["stock"] == "first dimer" else 3 / 8
            assert abs(sum(ctrl["edge_clicks"]) / n - target) < w
    for s in d["downstream"]:
        total = sum(sum(a["starts"]) for a in s["arms"])
        right = sum(sum(a["right"]) for a in s["arms"])
        assert (right - total / 2) / n == s["per_start"]
        assert abs(s["per_start"] - s["predicted"]) < s["width"]
        assert s["last_drain"] <= 2048
        if s["stock"] == "ground":
            assert s["releases"] == total == 0
        else:
            assert s["per_start"] + s["width"] < 0
