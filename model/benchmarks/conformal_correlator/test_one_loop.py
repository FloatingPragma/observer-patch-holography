from itertools import permutations
import json
from pathlib import Path

import numpy as np
import pytest

from model.benchmarks.conformal_correlator.one_loop import (
    REGULAR,
    SKEW,
    box,
    correlator,
    cross_ratio_box,
    run,
)


def test_local_insertion_and_symmetry_controls():
    for points in (REGULAR, SKEW):
        reference = cross_ratio_box(points)
        for perm in permutations(range(4)):
            assert box(points[list(perm)], 24) == pytest.approx(reference, rel=1e-12)
        assert correlator(points, 0)["connected_one_loop"] == 0
        translated = points + np.array([2, 3, 1, -1])
        np.testing.assert_allclose(
            correlator(translated, 0.01)["two_point"],
            correlator(points, 0.01)["two_point"],
            rtol=1e-13,
        )
    with pytest.raises(ValueError):
        box(np.zeros((4, 4)), 12)


def test_full_polarized_component_and_archive():
    saved = json.loads(Path(__file__).with_name("one_loop.json").read_text())
    for fresh, old in zip(run()["fixtures"], saved["fixtures"], strict=True):
        assert fresh["points"] == old["points"]
        for key in ("free", "interacting"):
            for observable, value in fresh[key].items():
                np.testing.assert_allclose(
                    value, old[key][observable], rtol=1e-12, atol=0
                )
        for key in (
            "relative_dilation_error",
            "relative_inversion_error",
            "relative_permutation_error",
            "relative_dilog_error",
            "relative_quadrature_change",
        ):
            assert max(fresh[key], old[key]) < 1e-9
        # Independent ratio from the two protected factors and the isolated R-symmetry monomial.
        points = np.array(fresh["points"])
        d12, d34 = [np.sum((points[i] - points[j]) ** 2) for i, j in ((0, 1), (2, 3))]
        loop = -cross_ratio_box(points) / 4
        prefactor = 2 * 3 / (4 * np.pi**2) ** 4
        expected = 0.01 * prefactor * loop / (d12 * d34)
        assert fresh["interacting"]["connected_one_loop"] == pytest.approx(
            expected, rel=1e-12, abs=0
        )
