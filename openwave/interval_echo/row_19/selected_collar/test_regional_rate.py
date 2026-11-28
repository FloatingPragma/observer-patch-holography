"""Check the independent whole-cell instrument and its count-only decisions."""

from copy import deepcopy
import json

import numpy as np
import pytest

from openwave.interval_echo.row_19.selected_collar import (
    regional_precision,
    regional_rate as region,
)


def test_shared_bond_cancels_and_other_bonds_are_measured():
    formation = region.camera.BEAM.front.optics.pairing.formation
    board = formation.federation(1)
    route = formation.region(board)
    left, right = (p // 12 for p in route[1:3])
    for parity in (False, True):
        effect = region.regional_difference(board, route, parity)
        for a, b in board.external:
            left_touch = left in (a // 12, b // 12)
            right_touch = right in (a // 12, b // 12)
            # A physical bridge has hopping -1/2. Each endpoint assigns one half
            # to its cell, and the meter reads half of the oriented difference.
            expected = (int(left_touch) - int(right_touch)) * np.eye(4) / 8
            idx = region.camera.BEAM.indices
            np.testing.assert_array_equal(
                effect[idx([a]), :][:, idx([b])].toarray(), expected
            )


def test_unresolved_first_bank_and_energy_sham():
    primary = json.loads(region.RECEIPT.read_text())
    for experiment in primary["banks"].values():
        assert region.decode(experiment["counts"]) == experiment["decoded"]
        assert all(
            not c["one_percent_established"]
            for c in experiment["decoded"]["comparisons"].values()
        )
    bank = deepcopy(primary["banks"]["8"]["counts"])
    bank["stocks"]["ordinary_plus"] = deepcopy(bank["stocks"]["ordinary_minus"])
    assert region.decode(bank)["decision"] == "unidentified ordinary calibration"
    bank["stocks"]["ordinary_plus"]["energy_plus"].pop()
    with pytest.raises(ValueError, match="Incomplete regional energy"):
        region.decode(bank)


def test_precision_and_population_on_unchanged_sources():
    primary = json.loads(region.RECEIPT.read_text())
    data = json.loads(regional_precision.RECEIPT.read_text())
    assert regional_precision.decode(primary, data) == data["decoded"]
    for name, comparison in data["decoded"]["comparisons"].items():
        for key, value in data["population"]["8"]["comparisons"][name].items():
            lo, hi = map(float, comparison[key].strip("[]").split(","))
            assert lo <= value <= hi
