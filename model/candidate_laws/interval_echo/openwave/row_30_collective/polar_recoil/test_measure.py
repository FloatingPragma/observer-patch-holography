import json

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from . import measure as testbed
from model.candidate_laws.interval_echo.openwave.row_22_confinement.mobile_string import material


def test_direct_word_evolution_at_two_wavelengths():
    raw = json.loads((testbed.HERE / "counts.json").read_text())
    n = 10
    bits = material.bits(n)
    walls = (bits != np.roll(bits, 1, axis=1)).sum(axis=1)
    h = material.generate(16.0 * walls + 8.0 * bits.sum(axis=1), n, .5, True)
    for mode in (2, -4):
        phase = np.exp(2j * np.pi * mode * np.arange(n) / n)
        initial = np.zeros(2**n, dtype=complex)
        sites = 1 << np.arange(n)
        initial[sites] = (1 + phase) / np.sqrt(2 * n)
        evolved = expm_multiply(-2j*h, initial)
        p = abs(evolved[sites])**2
        row = next(r for r in raw["rows"] if (r["n"], r["source"]) == (n, mode))
        bank = np.asarray(row["banks"][testbed.WAITS.index(2)]) / testbed.SHOTS
        np.testing.assert_allclose(bank, np.r_[p, 1-p.sum()], rtol=0, atol=.00065)


def test_unique_clock_winding_and_held_predictions():
    with pytest.raises(ValueError, match="windings"):
        testbed.unwrap((-9, 9), [1], [(1, .001)])
    raw = json.loads((testbed.HERE / "counts.json").read_text())
    result = testbed.read(raw)
    assert result["rest"][0] < 39.91317 < result["rest"][1]
    assert [(r["tests"],r["passed"]) for r in result["forecasts"]] == [(15,1),(15,3),(15,3)]
