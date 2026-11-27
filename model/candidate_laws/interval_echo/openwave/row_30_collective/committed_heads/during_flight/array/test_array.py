import json

import numpy as np
import pytest

from .. import measure as parent
from . import material, measure, readback


def test_covariance_vs_native_flight():
    """Two sources, two sites, every coin word, both click outcomes, both laws."""
    for k in (8, 0):
        waits = np.array([0.0, 0.3, 2.0])
        truth = material.joint(k, 5, 2, (3, 4), (1, 0), waits)
        for i, site in enumerate((3, 4)):
            for j, (_, bank) in enumerate(
                parent.histories(k, 5, 2, site, (1, 0), waits + 2)
            ):
                np.testing.assert_allclose(
                    truth[j, i], bank["spaced"].transpose(2, 0, 1), rtol=0, atol=3e-12
                )


def test_readback_census_and_record_failures():
    data = json.loads((measure.HERE / "results.json").read_text())
    with np.load(measure.HERE / "counts.npz") as f:
        raw = dict(f)
    rebuilt = readback.decode(**raw)
    readback.verify(rebuilt, data)
    assert data["total_wall_starts"] == 120259084288
    decision = data["decisions"]
    assert decision["both_frequency_bins_resolved"]
    assert decision["profiles_require_more_than_one_spatial_vector"]
    assert decision["single_real_wavenumber_necessary_test"] == ["excluded", "excluded"]
    assert decision["changed_law_zero_response"] == "compatible"
    assert decision["nominal_cut_only_explanation"] == "excluded"
    assert decision["nominal_intact_only_explanation"] == "unresolved"
    assert decision["cut_only_on_meter_box"] == "excluded"
    raw["counts"][0, 0, 3, 0, 0, 0] += 1
    with pytest.raises(ValueError, match="incomplete"):
        readback.decode(**raw)


def test_positive_wave_and_frequency_null_controls():
    x = np.arange(7)
    for k in (0.0, 0.43, 1.2, np.pi):
        wave = 0.2 * np.exp(1j * k * x) + (0.43 - 0.31j) * np.exp(-1j * k * x)
        interval = readback.wavelength(wave, 1e-6)
        assert interval[0] <= k <= interval[1]
    rng = np.random.default_rng(311)
    phase = rng.uniform(-4, 4, 49)
    m = 257
    w = (1 - np.cos(2 * np.pi * np.arange(m) / m)) / m
    direct = abs(np.exp(1j * phase[:, None] * np.arange(m)) @ w)
    np.testing.assert_allclose(readback.leakage(phase, m), direct, rtol=0, atol=2e-14)
    clean = measure.integer_draws(
        rng, np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]), 321
    )
    np.testing.assert_array_equal(clean, [[321, 0, 0], [0, 0, 321]])


def test_uniform_cut_enclosure_against_corner_eigenvectors():
    report = json.loads((measure.HERE / "results.json").read_text())
    box = report["cut_contact_sensitivity"]
    target = report["sector"]["intact"]["best_frequency"]
    clock = report["clocks"]
    for j in box["contact_J"]:
        for h in box["exchange_h"]:
            links = np.tile([2 * j, -2 * h], 10)[:-1]
            gen = np.diag(links, 1) - np.diag(links, -1)
            freq, modes = np.linalg.eigh(1j * gen)
            entry = abs(modes) ** 2
            assert entry.max() <= box["mode_entry_squared_upper"]
            for sign in (-1, 1):
                dt = clock["spacing"] + sign * clock["spacing_error"]
                response = readback.leakage(
                    (freq[:, None] - freq) * dt - target / 8, 4096
                )
                exact = np.sqrt(
                    20 * np.einsum("ap,aq,pq->a", entry[::2], entry[1::2], response**2)
                )
                assert exact.max() < box["cut_only_ceiling"]
