from fractions import Fraction as Q
import json

import numpy as np

from . import exterior_marks as analysis
from . import poisson_transfer


def test_exterior_floor_handles_map_edge_signed_increases_and_unlimited_inner_stock():
    # Set G times mapped baryonic mass to unity in the test units.
    g = analysis.transport.budget.local.profiles.galaxies.G * 1e9
    rows = dict(
        radius_kpc=[1.0, 2.0, 4.0],
        observed_kms=[1.0, 2.0, 2.0],
        error_kms=[0.01] * 3,
        baryons_kms=[0.0] * 3,
    )
    short = analysis.panel(rows, 2.0, 1 / g, Q(30), 5)
    assert short["unlimited_positive_anywhere"] and not short["positive_within_map"]
    assert short["strongest_exterior_pair"]["minimum_net_ratio"] > 7
    long = analysis.panel(rows, 4.0, 1 / g, Q(30), 5)
    assert long["positive_within_map"] and len(long["deposits"]) == 3
    assert analysis.tail(40) < analysis.tail(10) < 1
    assert analysis.reach(Q(4), Q(2), Q(1), 5)["inventory_insufficient"]


def test_weak_poisson_galaxy_normalization_is_held_out():
    inputs = {
        n: dict(
            observed_kms=[40.0, 30.0], baryons_kms=[20.0, 15.0], error_kms=[1.0, 1.0]
        )
        for n in ("a", "b", "c")
    }
    cols = {n: np.array([400.0, 225.0]) for n in inputs}
    original = poisson_transfer.scale_test(inputs, cols)
    inputs["c"]["observed_kms"] = [150.0, 120.0]
    changed = poisson_transfer.scale_test(inputs, cols)
    assert (
        original["galaxies"]["c"]["forecast_kms"]
        == changed["galaxies"]["c"]["forecast_kms"]
    )
    assert changed["galaxies"]["c"]["within_two"] == 0


def test_archived_outer_pair_is_beyond_entire_conditional_inventory():
    panel = json.loads((analysis.HERE / "exterior_marks.json").read_bytes())
    target = panel["galaxies"]["NGC2841"]
    pair = target["strongest_exterior_pair"]
    assert target["edge_kpc"] < pair["radius_kpc"] < pair["outer_kpc"]
    assert pair["minimum_net_ratio"] > float(Q(panel["cosmic_ceiling"]))


def test_outer_gas_individual_bands_and_joint_parameters():
    from . import dragonfly_profile as profile

    source = json.loads((profile.HERE / "dragonfly_gas.json").read_bytes())
    result = profile.run(source)
    assert len(result["held"]) == 7 and all(v["overlap"] for v in result["held"])
    assert not result["all_rings_share_one_fit"]
    assert result["shell_gain_comparison"]["required_ratio"] > 1.70
    assert result["stellar_normalization_check"]["ratio"] > 60


def test_mapped_pair_rate_and_withheld_galaxy_calibration():
    from . import pair_galaxies, pair_occupancy

    root = pair_occupancy.HERE
    radio = json.loads((root / "pair_occupancy_image.json").read_bytes())
    for row in pair_occupancy.evaluate(radio)["panels"].values():
        assert row["minimum_outer_to_inner_conversion"] > 5
        assert not row["joint_has_fit"]
    sample = json.loads((root / "pair_galaxy_maps.json").read_bytes())
    data = pair_galaxies.run(sample)["tests"]
    assert data["atomic"]["rings"] == data["mapped_baryons"]["rings"] == 369
    assert data["atomic"]["within_two"] == 93
    assert data["mapped_baryons"]["possible_full_curves"] == 0
