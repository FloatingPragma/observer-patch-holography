"""Reconstruct the optimal event and two-copy bound from all conditional word bins."""

from fractions import Fraction
from itertools import product
import json

import numpy as np

from ... import shared_variance, volume
from . import area_read, measure


def test_all_bins_and_joint_reader_bound():
    _, _, _, src = volume.radar_source_and_hand(measure.federation(0))
    laws, _ = shared_variance.source_words(src)
    reverse, fresh, shared = np.rint(laws * 2592).astype(int)
    delta = shared - fresh
    np.testing.assert_array_equal(delta, fresh - reverse)
    assert int(reverse[delta > 0].sum()) == 19 * 36
    assert sum(max(int(x), 0) for x in delta.flat) == 36
    saved = json.loads((area_read.HERE / "area_read.json").read_text())
    assert saved == area_read.run()
    assert Fraction(saved["conditional_total_variation_slope"]) == Fraction(1, 72)
    for row in saved["finite_panels"]:
        p = Fraction(row["eligible_probability"])
        # Three complete count bins: no informative arrival, informative with
        # positive D, and the remaining informative outcomes.
        ref = (1 - p, p * Fraction(19, 72), p * Fraction(53, 72))
        slope = (0, p / 72, -p / 72)
        signed = [
            (i, j, slope[i] * ref[j] + ref[i] * slope[j])
            for i, j in product(range(3), repeat=2)
        ]
        assert sum(abs(x) for _, _, x in signed) / 2 <= 2 * p / 72
        event = sum(x for i, j, x in signed if 1 in (i, j))
        assert event == 2 * p / 72 * (1 - ref[1])
        assert Fraction(row["optimal_event_slope"]) == p / 72
