from fractions import Fraction as F
import math

import pytest

from model.candidate_laws.interval_echo_partition_repair.regional import PartitionRepair
from model.candidate_laws.interval_echo_partition_repair.regional_check import EDGES, biased_control, finite_repair, refining_faces, region_receipt


def test_local_commands_and_actual_observer_restriction():
    law = PartitionRepair()
    assert all(law.advance_bond(law.advance_bond(i,0),3) == i for i in range(5))
    assert finite_repair()["two_step_probability_floor"] == "1/16"
    row = law.read_region((0, 2, 4, 1), (0, 1, 0, 1), EDGES, [0, 1])
    assert row.word == (0, 0, 1, 1, 0, 1) and row.open_boundary == 2 and row.volume == 1
    assert law.read_region((0, 2, 4, 1), (0, 1, 0, 1), EDGES, []).word == ()
    with pytest.raises(ValueError):
        law.advance_bond(-1, 0)
    with pytest.raises(ValueError):
        law.read_region((5,), (0,1), ((0,1),), (0,))


def test_all_regions_and_correlated_fixed_volume_transfers():
    result = region_receipt()
    assert result["global_words"] == 10000 and result["strict_restrictions_verified"] == 65
    assert result["proper_nonempty_subsets"] == 14
    for row in result["readbacks"]:
        for arm in row["arms"].values():
            assert F(arm["entropy_slope_in_bits"])+F(arm["aperture_slope"]) == F(arm["volume_slope"])
            assert arm["finite_deficit_nats"] + 1e-14 >= F(arm["pinsker_lower_nats"])
    arms = next(row for row in result["readbacks"] if row["cells"] == [0, 2])["arms"]
    assert F(arms["held_faces_and_pairing"]["aperture_slope"]) == F(8,5)
    corr = arms["correlation_only"]
    assert corr["aperture_slope"] == corr["volume_slope"] == corr["entropy_slope_in_bits"] == "0"
    assert corr["finite_deficit_nats"] > 0.02763
    assert result["calibrated_gain_bits_per_area_unit"] == "1"


def test_changed_repair_and_refining_finite_contrasts():
    control = biased_control()
    assert control["endpoint_frequency"] == ["1/6", "1/2", "1/3"]
    assert control["biased_pair_fixed_area_entropy_slope_nats"] < -0.4054
    panels = refining_faces()
    for panel in panels:
        assert panel["deficit_over_ell4"] >= float(F(panel["lower_over_ell4"]))
    assert math.isclose(panels[-1]["gain"], 1, abs_tol=5e-7)
    assert F(panels[-1]["lower_over_ell4"]) == 32
