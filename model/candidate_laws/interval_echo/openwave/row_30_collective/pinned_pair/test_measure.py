import copy
import json

import numpy as np
import pytest

from . import measure
from ...row_04 import open_pair


def test_phase_prediction_uses_original_all_port_source():
    sector, exchange, contrast, _, source = open_pair.full_parts()
    dark = np.zeros(sector.P)
    dark[source[[0, 2]]] = [1, -1]
    np.testing.assert_array_equal(exchange @ dark, dark)
    np.testing.assert_array_equal(contrast * dark, 3 * dark)
    vector = np.array([-13, 14]) / np.sqrt(365)
    for phase in measure.PHASES:
        kicked = measure.tip_kick(phase) @ vector
        assert abs(abs(kicked.sum())**2 / 2 - measure.bright_weight(phase)) < 2e-13


def test_original_receiver_and_holds():
    data = json.loads((measure.HERE / "counts.json").read_text())
    report = measure.read(data)
    assert report["banks"] == 120
    assert report["largest_error_in_allowances"] < 1
    assert max(report["native_pi"]) < min(report["exchange_only_pi"])
    assert min(report["native_pi"]) > max(report["native_zero"])
    # Move accepted quiet outcomes to a source word. Totals still balance.
    bad = copy.deepcopy(data)
    bank = bad["panels"][0]["banks"][3]
    bank[0] += bank[-2]
    bank[-2] = 0
    with pytest.raises(ArithmeticError, match="rejected"):
        measure.read(bad)
