"""Check the fractional phase and the recorded Ramsey fringes in the full port basis."""

from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).parent


def full_port_model(load, flux, momentum):
    words = list(combinations(range(10), load))
    index = {s: i for i, s in enumerate(words)}
    offset = (0, 4, 7) if load == 3 else (0, 1, 4, 7)
    source = np.zeros(len(words), complex)
    for x in range(10):
        word = tuple(sorted((k + x) % 10 for k in offset))
        source[index[word]] = np.exp(2j * np.pi * momentum * x / 10) / np.sqrt(10)
    h = np.zeros((len(words), len(words)), complex)
    for j, word in enumerate(words):
        for site in word:
            for step in (-2, -1, 1, 2):
                target = (site + step) % 10
                if target in word:
                    # Each occupied contact is visited in both directions.
                    h[j, j] -= 8
                else:
                    dest = tuple(sorted((set(word) - {site}) | {target}))
                    h[index[dest], j] -= np.exp(1j * step * flux / 10) / 2
                    h[j, j] += 0.5
    return h, source


@pytest.mark.parametrize("load", [3, 4])
def test_full_space_and_independent_count_decoder(load):
    result = next(
        r for r in json.loads((HERE / "dressing.json").read_text()) if r["load"] == load
    )
    for r in result["fresh"]:
        low, high = r["predicted_bounds"]
        counts = np.array(r["read"]["real_imag_plus_counts"])
        real, imag = 2 * counts / 262144 - 1
        slope = -np.polyfit(
            np.arange(len(real)), np.unwrap(np.angle(real + 1j * imag)), 1
        )[0]
        assert low <= r["spectral"] <= high
        assert high - low < 1e-6
        assert low - 0.00035 <= slope <= high + 0.00035
    r = result["fresh"][0]
    samples = [1, 8, 23]
    # No momentum block reduction or Schur decomposition enters this comparator.
    arms = []
    for phi in (r["flux"], -r["flux"]):
        h, source = full_port_model(load, phi, r["m"])
        energy, vectors = np.linalg.eigh(h)
        arms.append(
            (np.exp(-1j * np.outer(samples, energy)) * (vectors.conj().T @ source))
            @ vectors.T
        )
    overlap = np.sum(arms[1].conj() * arms[0], axis=1)
    counts = np.array(r["read"]["real_imag_plus_counts"])[:, samples]
    np.testing.assert_allclose(
        [overlap.real, overlap.imag], 2 * counts / 262144 - 1, atol=0.012
    )


def test_fractional_translation_and_competing_charges():
    studies = json.loads((HERE / "results.json").read_text())["studies"]
    dressed = json.loads((HERE / "dressing.json").read_text())
    for s, d, q in zip(studies[1:3], dressed, (-1 / 3, 2 / 3), strict=True):
        assert d["load"] == s["load"]
        np.testing.assert_allclose(
            [p["direct"] for p in d["rows"]],
            [p["rigid_direct"] for p in d["rows"]],
            atol=1e-12,
        )
        for proposal in s["errors"]:
            err = proposal["maximum_measured_error"]
            if proposal["charge"] == q:
                assert 0.01 < err < 0.0161
            elif proposal["charge"] in (-1, 0, 1):
                assert err > 0.33
    assert studies[0]["errors"][-1]["maximum_measured_error"] < 0.0001


def test_flux_cycle_gauge_permutation_and_independent_gap_sample():
    for panel in json.loads((HERE / "cycle.json").read_text()):
        load = panel["load"]
        for ramp in panel["ramps"]:
            expected = (ramp["initial_momentum"] - load) % 10
            for arm in ramp["inputs"].values():
                counts = np.array(arm["momentum_wall_other_counts"])
                assert counts[expected].sum() == counts.sum() == 262144
            if panel["contrast"] != 0:
                assert ramp["inputs"]["wall"]["wall_probability"] > 0.9775
                assert ramp["inputs"]["compact"]["wall_probability"] < 1e-6
        if panel["contrast"] == 0:
            continue
        for phi in (0.03, 1.234, 6.18):
            h, _ = full_port_model(load, phi, 1)
            energy = np.linalg.eigvalsh(h)
            gap = energy[-10] - energy[-11]
            assert gap >= panel["certificate"]["uniform_band_gap"] > 13.85
