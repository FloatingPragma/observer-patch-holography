import json
from fractions import Fraction
from math import sqrt

import numpy as np
import pytest

from .measure import (
    BENCH,
    HERE,
    N_CLOCK,
    RANGE,
    clock,
    fields,
    hamiltonian,
    instrument,
    station,
)
from .readback import read


def test_physical_capture_and_independent_clock():
    assert station(12) == (152, 150, 154, 212, 207)
    assert station(17) == (212, 207, 206, 106, 96)
    balanced = [(51 - sqrt(297)) / 4, 8.5, (51 + sqrt(297)) / 4]
    np.testing.assert_allclose(instrument(False)[0], balanced, atol=5e-14)
    # The half-integer entries are exact. det(H - x) vanishes at the three
    # distinct roots 0, 8 and 9; LAPACK returns them within its backward error.
    (a, b, c), (d, e, f), (g, h, i) = (
        [Fraction(x) for x in row] for row in hamiltonian(True).tolist()
    )
    for x in (0, 8, 9):
        assert (a - x) * ((e - x) * (i - x) - f * h) - b * (d * (i - x) - f * g) + c * (
            d * h - (e - x) * g
        ) == 0
    np.testing.assert_allclose(
        np.linalg.eigvalsh(hamiltonian(True)), [0, 8, 9], rtol=0, atol=1e-14
    )
    delays = np.arange(N_CLOCK) * 2 * np.pi / RANGE
    amplitudes = np.exp(-1j * np.outer(balanced, delays))
    explicit = abs(np.fft.ifft(amplitudes, axis=1)) ** 2
    np.testing.assert_allclose(
        clock(np.array(balanced)), explicit, atol=5e-13, rtol=1e-12
    )


def test_charge_sector_and_held_controls():
    data = json.loads((HERE / "counts.json").read_text())
    results = read(data)
    arm_results = [r["distinct"] for r in results["pairs"] if r["key"][2] == "arm"]
    winding_results = [
        r["distinct"] for r in results["pairs"] if r["key"][2] == "character"
    ]
    assert len(arm_results) == len(winding_results) == 12
    assert all(arm_results) and not any(winding_results)
    assert results["comparisons"] == {
        "position": 17,
        "one_contact": 24,
        "contrast_four": 24,
    }
    assert results["totals"]["native"]["accepted"] == 510459
    assert results["totals"]["pickup_parked"]["accepted"] == 0
    assert all(c == [0, 0, 0, BENCH, 0, 0, 0] for c in results["calibration_zero"])
    data["banks"].pop()
    with pytest.raises(ValueError, match="count bank"):
        read(data)


def test_coincidence_is_contained_in_sharp_pickup_at_each_mask():
    for by_mask in fields(8).values():
        for daughter, pickup in by_mask:
            selected = (abs(daughter) ** 2).sum(axis=1)
            assert np.all((selected > 0) & (selected < pickup))
