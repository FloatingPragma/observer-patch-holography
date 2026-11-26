"""Independent reconstruction, blinded-instrument controls and record integrity."""

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
import pytest

from ..measure import carrier, chart
from .read import DIRECTIONS, ENROLLMENT, read_bank, report


def records():
    return json.loads(Path(__file__).with_name("records.json").read_text())


def test_reconstructed_intervals_enclose_the_native_directional_symbols():
    archive = records()
    outcome = report(archive)
    pos, axes = chart()
    for bank in archive["banks"]:
        i, j = carrier.CELL_SEAMS[bank["seam"]]
        axis = bank["axis"]
        phase = 0 if axis == "zero" else 4 * (pos[i] - pos[j]) @ axes[axis]
        truth = (1 - np.cos(phase)) / 12
        status, (lo, hi), passed = read_bank(bank["batches"])
        assert passed and lo - 1e-14 <= truth <= hi + 1e-14
        assert status == ("unresolved" if bank["blind"] else "resolved")
    assert set(outcome["energies"]) == set(DIRECTIONS)
    assert outcome["face_minus_port"][0] > 0
    assert outcome["blind"] == ("unresolved", (0, 1), True)
    assert outcome["half_checks_passing"] == 151


def test_changing_only_source_labels_can_counterfeit_the_entire_count_law():
    p = Fraction

    def response(exit_, exposed_offset):
        dark, gain = p(1, 8), p(1, 2)
        return [
            dark,
            dark + gain,
            dark + gain * exit_ + exposed_offset,
            dark + gain / 2,
        ]

    assert response(p(1, 4), 0) == response(p(1, 8), p(1, 16))


def test_complement_counts_and_failed_standards_are_retained():
    raw = records()
    with pytest.raises(ValueError, match="Incomplete"):
        report({**raw, "banks": raw["banks"][:-1]})
    with pytest.raises(ValueError, match="Duplicate"):
        report({**raw, "banks": raw["banks"] + raw["banks"][:1]})
    broken = deepcopy(raw["banks"][0]["batches"])
    broken[0][1][1] += 1
    with pytest.raises(ValueError, match="enrollment"):
        read_bank(broken)
    for batch in raw["banks"][0]["batches"]:
        n = sum(batch[3])
        batch[3] = [n, 0]
    assert report(raw)["energies"]["zero"] is None
    impossible = [
        [[0, ENROLLMENT // 32], [ENROLLMENT // 32, 0], [ENROLLMENT // 16, 0], [0, 0]]
        for _ in range(8)
    ]
    assert read_bank(impossible)[0] == "inconsistent"
