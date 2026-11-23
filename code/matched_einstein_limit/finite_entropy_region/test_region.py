import json
from pathlib import Path

import numpy as np

from . import measure


def test_repair_weights_and_equal_volume_are_independently_readable():
    archive = json.loads(Path(__file__).with_name("results.json").read_text())
    states = [(n, b, w) for n in measure.TILES for b in range(2) for w in range(2**n)]
    fair = measure.repair_kernel(states, 1.)
    np.testing.assert_array_equal(fair, fair.T)
    for sector in archive["geometry"]:
        assert abs(sector["volume_reading"]/archive["target_volume"]-1) < .0002
        assert abs(sector["area_reading"]/np.pi-sector["tiles"]) < .001
    for row in archive["observed"]:
        assert sum(row["terminal_counts"]) == archive["trials_per_launch"]
        assert sum(row["joint_counts_by_sector_material_bit_edge_word"]) == archive["trials_per_launch"]
        if row["repair"] == "reversible":
            assert row["maximum_absolute_derivative"] < .03
        else:
            assert row["maximum_absolute_derivative"] > .6
