import json

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.regional_reference import (
    measure,
)


def test_all_symmetric_outcomes_and_resolvable_one_sided_cut():
    results = measure.run()
    for panel in results["temperature_panels"]:
        for key, cut in panel["cuts"].items():
            if key.startswith("both_"):
                np.testing.assert_allclose(
                    cut["probabilities"][1], cut["probabilities"][2], atol=5e-12
                )
    side = results["temperature_panels"][1]["cuts"]["first_1"]
    assert 0.37 < side["uniform_resolution_total_variation"]["0.01"] < 0.38


def test_independent_source_and_flight_readings():
    archive = json.loads((measure.HERE / "two_reads.json").read_text())
    assert archive["maximum_matrix_residual"] < 2e-12
    for row in archive["held"].values():
        error = np.abs(
            np.array(row["observed"]) - row["inferred_from_charge_and_coincidence"]
        )
        assert np.all(
            error < archive["modular_calibration"]["prediction_99_half_width"]
        )

    def at_four(row):
        return archive["held"][row]["inferred_from_source_counts"][2]

    assert at_four("target") > at_four("paired") > at_four("phase_feedback")
    m = archive["modular_calibration"]
    assert (
        m["dark_mean"] - m["bright_mean"] > 1.1 and m["unwrite_99_half_width"] < 0.001
    )
