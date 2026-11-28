import json

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    source_reader as reader,
)


def test_count_only_nested_comparison_and_inversion():
    data = json.loads((reader.HERE / "source_reader.json").read_bytes())
    for lab in ("A", "B"):
        for k in (0, 8):
            case = data["cases"][f"{k}_{lab}"]["decoded"]
            truth = 32 * (k + 0.5) * np.log(8 / 7)
            for tag in ("folded", "folded_reversed"):
                record = case[tag]
                first, second = record["balls"]
                assert first == second and first[1] < truth < first[2]
                assert record["outer_minus_twice_inner"][1] < 0
                np.testing.assert_allclose(record["two_radius_inversion"], first)
                assert (
                    record["diamond_weighted_balls"][1][1]
                    > 2 * record["diamond_weighted_balls"][0][2]
                )
            assert case["blank"]["balls"] == [[0, 0, 0]] * 2
            assert case["outward"]["balls"][0][2] < case["folded"]["balls"][0][1]
