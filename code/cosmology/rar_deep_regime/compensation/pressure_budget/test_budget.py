import numpy as np

from cosmology.rar_deep_regime.compensation.pressure_budget import measure


def test_boundary_and_isolated_profiles():
    for row in measure.profile():
        assert row["traction_mass"] <= row["exterior_mass"]
    for row in measure.outer_wall():
        np.testing.assert_allclose(row["lens_over_mass"], 1 + row["boundary_over_mass"])
    assert measure.outer_wall()[-1]["lens_over_mass"] == 2.5


def test_observed_all_bins_and_simultaneous_alternative():
    cases = measure.run()["cases"]
    usual = cases["two_systematic_double_baryons"]
    np.testing.assert_allclose(
        usual["additional_mass_over_baryons"],
        [6.3442, 6.1046, 8.0427, 5.3987],
        atol=0.00006,
    )
    assert usual["closed_stock_layer_equivalents"]["native_external_coherent"] == [
        15,
        15,
        19,
        13,
    ]
    assert (
        min(
            cases["bonferroni_four_normal_95_double_baryons"][
                "additional_mass_over_baryons"
            ]
        )
        > 5
    )


def test_all_commissioned_supplies_in_galaxy_stacks():
    trails = measure.run()["commissioned_trails"]
    assert [t["conserved_ratio"] for t in trails] == ["425/202", "1037/202", "2125/202"]
    assert [
        sum(
            t["cases"]["two_systematic_double_baryons"]["can_supply_optical_lower_mass"]
        )
        for t in trails
    ] == [0, 0, 4]
    assert (
        sum(
            trails[1]["cases"]["bonferroni_four_normal_95_double_baryons"][
                "can_supply_optical_lower_mass"
            ]
        )
        == 1
    )
