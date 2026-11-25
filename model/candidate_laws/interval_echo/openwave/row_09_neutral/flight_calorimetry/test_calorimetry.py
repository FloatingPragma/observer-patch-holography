"""Independent native generator, actual delay paths, root brackets and all count banks."""

import json
import numpy as np
from scipy.linalg import eigh

from .apparatus import APPOINTMENT, DT, RING, SENSORS, clock, flight, loads, sector
from .certificate import exact_support
from .measure import HERE, STARTS, LOCAL_STARTS


def test_energy_and_decay_certificate_from_native_release():
    states, h = flight.chart((160, 167, 0), False)
    words, native, _, _ = sector("pair", 8)
    order = [states.index((*w, 0)) for w in words]
    np.testing.assert_array_equal(native, h[np.ix_(order, order)])
    bounds = exact_support()
    assert bounds["roots_per_bracket"] == [1, 1]
    assert bounds["charged_daughter_plus_pair_contrast_histogram"] == {
        32: 140,
        48: 658,
        64: 280,
        80: 14,
    }
    assert bounds["all_pair_contrast_histogram"] == {32: 154, 48: 28}


def test_coherent_filter_against_full_space_at_small_clock_size():
    # Direct cache-free 182-state summation along the separately translated source,
    # with no orbit reduction or assumed scalar energy probabilities.
    n = 19
    accepted, table = clock("pair", 8, True, -1, n)
    words, h, _, _ = sector("pair", 8)
    indices = [words.index((a, b)) for a, b in zip(RING, RING[1:] + RING[:1])]
    assert (160, 167) in [words[j] for j in indices]
    seed = np.zeros(len(words))
    seed[indices] = 1 / len(RING)
    ev, U = eigh(h)
    states = U @ (
        np.exp(-1j * np.outer(ev, APPOINTMENT + DT * np.arange(n)))
        * (U.T @ seed)[:, None]
    )
    out = np.fft.ifft(states, axis=1)
    assert abs(accepted - 1 / len(RING)) < 2e-14
    for name, edge in SENSORS.items():
        d = np.array([np.subtract(*loads(w, edge)) for w in words])
        for value in (-2, -1, 0, 1, 2):
            np.testing.assert_allclose(
                table[name][:, value + 2],
                np.sum(abs(out[d == value]) ** 2, axis=0),
                atol=4e-13,
            )


def test_joint_counts_and_native_soft_line_controls():
    result = json.loads((HERE / "results.json").read_text())
    with np.load(HERE / "counts.npz", allow_pickle=False) as banks:
        assert len(banks.files) == len(result["timing"]) == 90
        for row in result["timing"]:
            name = row["key"]
            bank = banks[name]
            assert sum(bank) == row["starts"] == STARTS
            assert list(bank[-3:]) == row["refusals"]
            clocks = bank[:-3].reshape(-1, 5)
            assert clocks.sum() == row["accepted"]
            assert clocks[:, [0, 4]].sum() == row["contrast_two"]
            if name.startswith("native_") and "_1_4096_" in name:
                assert row["soft"] == 0 and row["high_share_99"][0] > 0.9998
                assert {b for b, n in row["peak_bins"]} == {49.0625, 50.0625}
            if name.startswith("exchange_") and "_4096_" in name:
                assert row["soft_share_99"][0] > 0.60
            if name.startswith("depleted_"):
                assert row["contrast_two"] == 0 and row["peak_bins"][0][0] == 16
    for entry in result["direct"]:
        assert entry["post_tag_contrast_99"][0] > 47.95
        assert all(sum(c["counts"]) == LOCAL_STARTS for c in entry["cohorts"])
