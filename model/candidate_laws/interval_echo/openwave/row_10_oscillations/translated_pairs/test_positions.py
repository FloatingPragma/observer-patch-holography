import json
from copy import deepcopy

import numpy as np
import pytest

from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release.flight import (
    chart,
)
from . import measure


def test_native_ring_matches_existing_release_generator():
    ring = measure.layout()[0]
    basis, h = measure.generator(ring, 8)
    physical, reference = chart((ring[0], ring[1], 0), False)
    indices = [basis.index((ring.index(a), ring.index(b))) for a, b, _ in physical]
    np.testing.assert_array_equal(h[np.ix_(indices, indices)], 8 * reference)


@pytest.mark.parametrize(
    "fixed,name,expected",
    [
        (8, "counts.json", [True, False, True]),
        (7, "second_counts.json", [True, True, True]),
    ],
)
def test_physical_ports_and_all_outcomes(fixed, name, expected):
    record = json.loads((measure.HERE / name).read_text())
    assert record["fixed_bond"] == fixed
    for result in measure.decode(record)["strengths"].values():
        assert result["covered"] == {"fixed": 36, "translated": 36}
        assert result["sham_events"] == 0
        assert [
            pair["distinguishable"] for pair in result["source_comparisons"]
        ] == expected


def test_translation_prediction_can_fail_with_valid_counts():
    record = json.loads((measure.HERE / "second_counts.json").read_text())
    changed = deepcopy(record)
    bank = next(
        row for row in changed["banks"] if row[:5] == [8, 64, "translated", 1, 8]
    )
    bank[-1] = [measure.EXPOSURE, 0, 0]
    assert measure.decode(changed)["strengths"]["8"]["covered"]["translated"] < 36
    changed["banks"].pop()
    with pytest.raises(ValueError, match="missing"):
        measure.decode(changed)
