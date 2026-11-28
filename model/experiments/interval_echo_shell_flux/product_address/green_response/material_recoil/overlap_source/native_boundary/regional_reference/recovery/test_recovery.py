"""Check contraction under both channels, implemented preparation, and held read signs."""

import json

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.regional_reference.recovery import (
    measure as m,
)


def test_joint_orbit_and_native_conditional_route():
    states, _, d, u = m.preparation()
    assert len(m.bypass()["parent_bypass"]) == 18
    output = m.analysis(np.log(8 / 7), states, d, u)["states"]
    assert (
        output["E"]["thinning"][3]["unflagged_loss"]
        < output["P"]["thinning"][3]["unflagged_loss"]
        < output["FG0"]["thinning"][3]["unflagged_loss"]
    )
    for name in ("P", "E", "FG0", "FG4"):
        s = output[name]
        np.testing.assert_allclose(
            s["thinning"][-1]["flagged_loss"],
            s["thinning"][-1]["unflagged_loss"],
            atol=1e-13,
        )
    np.testing.assert_allclose(
        output["E"]["connected_modular"], output["FG4"]["connected_modular"], atol=1e-13
    )


def test_certified_opposite_material_signs_and_independent_clocks():
    actual = m.run()
    m.region.check(json.loads((m.HERE / "results.json").read_text()), actual)
    read = actual["probes"]
    assert read["certified_G_time_four"]["outward_billionths"] == [170790276, 170790329]
    assert read["flight"]["FG0"]["mean"][2] < -0.134
    assert read["flight"]["FG4"]["mean"][2] > 0.1707
    assert read["independent_clock_max_error"] < 1e-9
    assert read["clock_contrasts"]["FG0"]["intervals"][0][1] > 0
    assert read["clock_contrasts"]["FG4"]["intervals"][1][1] < 0
    for row in read["count_banks"].values():
        np.testing.assert_array_equal(np.array(row["same"]) + row["opposite"], 2**30)
