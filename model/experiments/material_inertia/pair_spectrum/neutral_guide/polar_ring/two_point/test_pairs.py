import json
import numpy as np
from . import measure, pairs


def test_energy_bins_cover_sharp_square_and_false_tags():
    for offset in (0.0, 1.0, -1.0):
        np.testing.assert_allclose(
            pairs.window(np.array([80.0 + offset]), 80.0), 1.0, atol=1e-12
        )
    assert pairs.window(np.array([78.0, 82.0, 240.0]), 80.0).max() < 0.007


def test_complete_pair_banks_and_frozen_forecasts():
    data = pairs.collect()
    assert json.loads(json.dumps(data)) == measure.read(measure.HERE / "pairs.json")
    for p in data["panels"]:
        for row in p["spatial"]:
            assert all(
                sum(c) == pairs.SHOTS
                for c in [row["reference_left"], row["reference_right"], *row["pairs"]]
            )
        assert all(sum(c) == pairs.SHOTS for c in p["energy"])
        assert pairs.conditional(p["energy"][0])[1] > 0.9
        assert pairs.conditional(p["energy"][-1])[0] > 0.9
        tv = pairs.compare(p)
        assert tv[2, 0] > 0.05 and tv[2, -1] < 0.05
        assert tv[5, -1] > 0.05
