"""Check an independent geometric example and both instruments' experimental decisions."""

from copy import deepcopy
import json

import numpy as np
import pytest

from openwave.interval_echo.row_19.screen_rulers import measure, radar_precision


def test_triangle_with_known_area_and_orientation():
    # Every return at time 3 gives side length 1.5. An identity detector's
    # three centroids coincide with the three target vertices.
    returns = [float(k == 6) for k in range(2, 13)]
    panel = {
        "clock": [float(np.sin(t / 2) ** 2) for t in measure.TIMES],
        "range": [returns] * 6,
        "image": np.eye(4)[:3].tolist(),
    }
    result = measure.panel_read(panel, False)
    expected = 1.5**2 * np.sqrt(3) / 4
    for k in ("radar", "image"):
        assert float(result[k].mid) == pytest.approx(expected, abs=1e-13)
    panel["image"].reverse()
    assert float(measure.panel_read(panel, False)["image"].mid) == pytest.approx(
        -expected, abs=1e-13
    )


def test_complete_held_image_and_clock_swap():
    data = json.loads(measure.RECEIPT.read_text())
    banks = data["banks"]
    for b in banks.values():
        assert measure.decode(b["counts"], data["half_steps"], True) == b["counted"]
    assert (
        banks["8"]["counted"]["tests"]["image"]["comparison"]["decision"] == "different"
    )
    assert (
        banks["0"]["counted"]["tests"]["image"]["comparison"]["decision"]
        == "unresolved"
    )
    active = deepcopy(banks["8"]["counts"])
    for stock in measure.STOCKS:
        vacancy, excess = active["vacancy"][stock], active["excess"][stock]
        c1, c2 = vacancy["clock"], excess["clock"]
        vacancy["clock"], excess["clock"] = c2, c1
    crossed = measure.decode(active, data["half_steps"], True)
    assert crossed["tests"]["image"]["comparison"]["decision"] == "different"
    # The same complete trial records in both arms should not distinguish probes.
    active["excess"] = deepcopy(active["vacancy"])
    sham = measure.decode(active, data["half_steps"], True)
    assert all(
        v["comparison"]["decision"] == "unresolved" for v in sham["tests"].values()
    )
    active["excess"]["held_plus"]["range"][2].pop()
    with pytest.raises(ValueError, match="Incomplete range"):
        measure.decode(active, data["half_steps"], True)


def test_precision_radar_and_held_source_marginals():
    data = json.loads(measure.RECEIPT.read_text())
    precise = json.loads(radar_precision.RECEIPT.read_text())
    assert radar_precision.decode(data, precise) == precise["decoded"]
    stocks, _, recipes = measure.preparations()
    rp = stocks["reference"]
    for name in ("held_plus", "held_minus"):
        np.testing.assert_allclose(
            measure.source.marginals(stocks[name]),
            measure.source.marginals(rp),
            atol=1e-14,
            rtol=0,
        )
        assert recipes[name]["min_kept"] > 0 and recipes[name]["min_dephased"] > 0
