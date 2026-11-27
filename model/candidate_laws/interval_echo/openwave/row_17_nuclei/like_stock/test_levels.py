import json
from itertools import combinations

import numpy as np
import pytest

from model.executable_model import federation
from . import levels


def test_finite_clock_against_coherent_delay_paths():
    t = np.arange(levels.DELAYS) * 2 * np.pi / levels.PERIOD
    e = np.array([0.0, 69.9453125, 77.092803, 511.99999])
    time_read = np.mean([abs(np.fft.ifft(np.exp(-1j * x * t))) ** 2 for x in e], axis=0)
    np.testing.assert_allclose(levels.population(e), time_read, atol=1e-12, rtol=1e-7)


def independent_spectrum(k):
    # Integer mesh assembled from addresses of the actual federation; no QPE,
    # local projector calibration or sign-fragment class enters this check.
    board = federation(0)
    edges = [(a, b) for a, b in board.seams if max(a, b) < 12]
    referenced = set(np.flatnonzero(board.degree[:12] == 5))
    words = [frozenset(s) for s in combinations(range(12), 3)]
    address = {w: i for i, w in enumerate(words)}
    twice = np.zeros((220, 220), dtype=np.int64)
    for i, w in enumerate(words):
        twice[i, i] = -2 * k * len(w & referenced)
        for a, b in edges:
            if (a in w) != (b in w):
                twice[i, i] += 2 * k + 1
                twice[address[w ^ {a, b}], i] -= 1
    assert np.array_equal(twice, twice.T)
    return np.linalg.eigvalsh(twice) * 0.5


def test_registered_spectrum_and_exclusions():
    record = levels.load()
    survey = json.loads((levels.HERE / "counts.json").read_text())
    results = levels.decode(record, survey)
    assert [len(r["windows"]) for r in results["rows"]] == [8, 6, 0]
    for row in results["rows"]:
        spectrum = independent_spectrum(row["k"])
        lo, hi = row["energy_support"]
        assert lo - 1e-10 <= spectrum[0] <= spectrum[-1] <= hi + 1e-10
        previous = -1.0
        for w in row["windows"]:
            lower, upper = w["window"]
            assert lower > previous and upper < row["breakup_lower"]
            assert np.any((spectrum >= lower) & (spectrum <= upper))
            assert w["bright"] >= results["min_counts"]
            previous = upper
    record["rows"][0]["counts"][0] += 1
    with pytest.raises(ValueError):
        levels.decode(record, survey)
