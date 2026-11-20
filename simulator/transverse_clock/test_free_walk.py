"""Test the free-law source map and the inclusive registered detector readouts."""
import numpy as np
import pytest

from simulator import law
from tools.strict_json import load
from . import equivalence, frontier
from .receipts import HERE, verify_freeze


def test_each_free_slot_column_under_fourier_transform():
    carrier = (2, 3, 4, 5)
    angles = np.stack(np.meshgrid(*(2 * np.pi * np.fft.fftfreq(n) for n in carrier),
                                  indexing="ij"), axis=-1)
    # Translation covariance fixes every other source address from these eight impulses.
    for slot in range(law.SLOTS):
        state = np.zeros((*carrier, law.SLOTS), complex)
        state[0, 0, 0, 0, slot] = 1
        evolved = np.fft.fftn(law.step(state, law.LAWS["oph"]), axes=(0, 1, 2, 3))
        expected = np.array([law.walk_matrix(k)[:, slot] for k in angles.reshape(-1, 4)])
        np.testing.assert_allclose(evolved, expected.reshape(state.shape), atol=2e-14, rtol=0)


@pytest.mark.parametrize("factor", [-1, np.nan])
def test_fiber_check_rejects_changed_pinned_action(monkeypatch, factor):
    original = law.walk_matrix
    monkeypatch.setattr(law, "walk_matrix", lambda k, phase: factor * original(k, phase))
    with pytest.raises(ValueError, match="Fourier action"):
        equivalence.verify_fourier({"sizes": [16], "panels": {"x": [3, 4]}, "halfwidth": 2},
                                   {"sizes": []})


def test_passive_complement_against_all_literal_slots():
    t, speeds = 9, [(2, 3), (1, 1)]
    tails = frontier.exact_tails([t], speeds)
    row = equivalence.passive_outcomes([t], speeds, tails)[str(t)]
    state = np.zeros((2 * t + 1, 1, 1, 1, law.SLOTS), complex)
    state[t, 0, 0, 0] = law.EQUAL
    for _ in range(t):
        state = law.step(state, law.LAWS["oph"])
    amps = state[:, 0, 0, 0]
    equal = np.abs(amps @ law.EQUAL) ** 2
    other = (np.abs(amps) ** 2).sum(axis=1) - equal
    for a, b in speeds:
        cut = t + (a * t + b - 1) // b
        counts = list(map(int, row["thresholds"][f"{a}/{b}"].values()))
        denominator = int(row["denominator"])
        assert sum(counts) == denominator
        np.testing.assert_allclose(np.array(counts) / denominator,
                                   [equal[cut:].sum(), equal[:cut].sum(),
                                    other[cut:].sum(), other[:cut].sum()], atol=2e-15)
    tails[str(t)]["1/1"]["numerator"] = "0"
    with pytest.raises(ValueError, match="passive tail differs"):
        equivalence.passive_outcomes([t], speeds, tails)


def test_frozen_decisions_and_detector_partitions():
    verify_freeze()
    counts = load(HERE / "records/result.json")
    for size in (128, 256, 512, 1024):
        assert all(counts["clock"][f"{size}_{label}"]["decision"] == "within tolerance"
                   for label in ("out", "reverse", "fast"))
        assert counts["clock"][f"{size}_short"]["decision"] == "outside tolerance"
        assert counts["clock"][f"{size}_erased"]["decision"] == "unresolved"
    protocol = load(HERE / "local_clock/protocol.json")
    results = load(HERE / "local_clock/records/result.json")
    with np.load(HERE / "local_clock/records/opened.npz", allow_pickle=False) as detectors:
        for size in protocol["sizes"]:
            shots = protocol["shots_per_four_volume"] * size ** 4
            for label in protocol["center"]:
                key = f"{size}_{label}"
                expected = "outside tolerance" if label == "short" or key == "96_reversed" else "within tolerance"
                assert results[key]["decision"] == expected
                assert results[key]["maximum_count_difference"] == 0
                for name, columns in ((key, 9), (key + "_unwritten", 2)):
                    bins = detectors[name]
                    assert bins.shape[1] == columns
                    assert np.all(bins >= 0) and np.all(bins.sum(axis=1) == shots)
    angled = results["384_angled"]["observed"]
    assert angled["counts"] == 815120 and 5.71e-10 < angled["fraction"] < 5.73e-10
