import numpy as np

from .experiment import run


def test_calibrated_joint_source_on_native_refinement():
    result = run(include_finest=False)
    assert [patch["cells"] for patch in result["addressed_patches"]] == [20, 80, 320]
    assert [full["cells"] for full in result["all_links"]] == [20, 80]
    for case in (*result["addressed_patches"], *result["all_links"]):
        assert case["joint_record_total_variation"] > 1e-7
        rates = case["objects"]
        assert abs(rates["charge"]["contact_clock_rate_change"]) > 1e-8
        assert abs(rates["cage"]["contact_clock_rate_change"]) > 1e-8
        for record in rates.values():
            np.testing.assert_allclose(record["rest_clock"]["phase_clock_rate"], 1, atol=3e-9)
            np.testing.assert_allclose(record["initial_radial_velocity"], 0, atol=3e-9)
            assert record["initial_centroid_proper_time_rate"] == 1
            assert 0 < record["initial_mass_response_rate"] < 1
            assert abs(record["mean_radial_drift"]) > 1e-5
            assert abs(record["contact_phase_rate_change"]) > record["contact_phase_rate_change_error_bound"]
            assert 0 <= record["clock_velocity_disk_read"] <= 1 + 3e-9
