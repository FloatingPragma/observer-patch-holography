from math import acos, pi
import json

import numpy as np
import pytest
from mpmath import iv
from scipy.optimize import linprog

from . import geometry as g, readout, smooth


def test_actual_station_enclosures_and_source_replay():
    saved = json.loads(g.REPORT.read_text())
    iv.dps = 35
    original = g.raw_bounds()
    assert g.certify(original) == {
        k: saved[k] for k in ("centered_gram_lower", "original_four_volume", "fillings")
    }
    assert readout.clock_floor(original) == saved["clock_fit"]
    assert smooth.certificate() == json.loads((g.HERE / "smooth.json").read_text())
    sample = saved["sampled"]
    counts = np.asarray(sample["terminal_counts"])
    assert counts.shape == (30, 11)
    assert np.all(counts.sum(axis=1) == sample["starts_per_pair"]) and np.all(
        counts >= 0
    )
    for k, (i, j) in enumerate(readout.ORDER):
        mean = counts[k] @ np.arange(2, 13) / (4 * sample["starts_per_pair"])
        assert abs(mean - float(original[i, j].mid)) < sample["clock_radius"] + 1.81e-6
    for data in [saved, sample]:
        first, last = data["fillings"][1], data["fillings"][5]
        assert first["solid_angles"][0][1] < last["solid_angles"][0][0]
        assert first["total_volume"][0] > last["total_volume"][1]


def test_clock_cycle_floor_against_direct_minimax():
    original = g.raw_bounds()
    rows, bound = [], []
    for i, j in readout.ORDER:
        row = np.zeros(7)
        row[i], row[j], row[-1] = 1, -1, -2
        rows.append(row)
        bound.append(np.log(float(original[i, j].mid) / float(original[j, i].mid)))
    fit = linprog(
        [0] * 6 + [1],
        A_ub=rows,
        b_ub=bound,
        bounds=[(0, 0)] + [(None, None)] * 5 + [(0, None)],
        method="highs",
    )
    assert fit.success
    lo, hi = readout.clock_floor(original)["log_error"]
    assert lo < fit.fun < hi
    # Changing the first station's unit of time adds an exact gradient to the log ratios.
    scaled = original.copy()
    scaled[0, :] *= iv.mpf("1.001")
    assert readout.clock_floor(scaled)["log_error"] == [lo, hi]


def test_regular_simplex_and_impossible_length():
    means = np.array(
        [[iv.mpf(9) if i != j else iv.mpf(0) for j in g.ALL] for i in g.ALL],
        dtype=object,
    )
    rows, volume = g.fillings(means)
    reference_angle = 6 * (3 * acos(1 / 3) - pi)
    assert len(rows) == 6
    for row in rows:
        assert row["solid_angles"][0][0] <= reference_angle <= row["solid_angles"][0][1]
    assert volume[0] <= 27 * 2**0.5 / 12 <= volume[1]
    means[1, 0] = means[0, 1] = iv.mpf(100)
    with pytest.raises(ValueError, match="simplex"):
        g.certify(means)


def test_fresh_clock_banks_and_parked_law():
    assert readout.collect() == json.loads(g.REPORT.read_text())
