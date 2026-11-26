import json

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.row_22_confinement.mobile_string import (
    material,
)
from . import measure, wall


def test_native_word_propagation_at_zone_edge():
    records = measure.load()
    for ring in records["materials"]:
        n = len(ring["ring"])
        bits = material.bits(n)
        walls = np.sum(bits != np.roll(bits, 1, axis=1), axis=1)
        for strength in (4, 8):
            hamiltonian = material.generate(
                strength * (2 * walls + bits.sum(axis=1)).astype(float), n, 0.5, True
            )
            for sign in (1, -1):
                carriers = [sign * (n // 2 - 1) % n, n // 2]
                sites = 1 << np.arange(n)
                wave = np.zeros(2**n, dtype=complex)
                wave[sites] = sum(measure.character(n, k) for k in carriers) / np.sqrt(
                    2
                )
                evolved = expm_multiply(-2j * hamiltonian, wave)
                law = abs(evolved[sites]) ** 2
                row = next(
                    r
                    for r in ring["flights"][str(strength)]["rows"]
                    if r["carriers"] == carriers
                )
                counts = np.array(row["banks"][measure.WAITS.index(2)])
                np.testing.assert_allclose(
                    counts / measure.SHOTS,
                    np.r_[law, 1 - law.sum()],
                    rtol=0,
                    atol=0.0015,
                )


def test_all_sources_count_completeness_and_velocity_results():
    data = measure.load()
    result = measure.read(data)
    report = json.loads((measure.HERE / "results.json").read_text())
    assert set(result) == set(report["polar"])
    assert sum(len(row["speed"]) for row in result.values()) == 44
    for key, row in result.items():
        np.testing.assert_allclose(
            row["speed"], report["polar"][key]["speed"], rtol=0, atol=1e-12
        )
        assert row["near_rest_increase_at_least"] > 0.0219
        assert row["near_edge_decrease_at_least"] > 0.022
    data["materials"][0]["flights"]["8"]["rows"].pop()
    with pytest.raises(ValueError, match="Missing outward"):
        measure.read(data)
    data["materials"][0]["sources"]["1"].reverse()
    with pytest.raises(ValueError, match="Source preparation"):
        measure.read(data)


def test_all_archived_wall_flights():
    results = wall.compare()
    report = json.loads((measure.HERE / "results.json").read_text())["wall"]
    assert set(report) == set(results)
    for key, result in results.items():
        assert len(result["rows"]) == 7
        for row, expected in zip(result["rows"], report[key]["rows"], strict=True):
            for field, val in row.items():
                np.testing.assert_allclose(val, expected[field], rtol=0, atol=1e-12)
    assert sum(row["low_continuum"] for row in results.values()) == 9
    assert sum(row["high_continuum"] for row in results.values()) == 0
    assert sum(row["low_lattice"] for row in results.values()) == 12
    assert sum(row["high_lattice"] for row in results.values()) == 15
