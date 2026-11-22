"""Exercise label blindness, complete clock outcomes and the changed-law control."""

import json

import numpy as np
import pytest

from .blind_chart import BOOK, HERE, infer


def test_decoder_invariant_under_anonymous_relabeling():
    data = json.loads((HERE / "blind_chart.json").read_text())
    counts = np.array(data["calibration_counts"])
    coords, compressed = infer(counts)
    perm = np.array([3, 0, 2, 7, 1, 4, 6, 5])
    relabeled, matrix = infer(counts[np.ix_(perm, perm)])

    def distances(p):
        return np.count_nonzero(p[:, None] != p[None, :], axis=2)

    np.testing.assert_array_equal(
        distances(relabeled), distances(coords)[np.ix_(perm, perm)]
    )
    np.testing.assert_array_equal(matrix, compressed[np.ix_(perm, perm)])
    assert coords.shape == (8, 3)
    with pytest.raises(ValueError):
        infer(np.ones_like(counts))


def test_fresh_clock_outcomes_and_shift_only_control():
    data = json.loads((HERE / "blind_chart.json").read_text())
    protocol = json.loads(BOOK.read_text())
    assert min(data["removed_direction_error_lower_bounds"]) > 0.11
    for row in data["held"]:
        counts = np.array(row["detector_counts"])
        assert counts.shape == (8, 10)
        assert np.all(counts.sum(axis=1) == protocol["shots"])
        assert (row["maximum_count_error"] <= protocol["tolerance"]) == (
            row["arm"] == "nominal"
        )
