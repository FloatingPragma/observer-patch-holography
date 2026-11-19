import numpy as np

from .ring_clock import run


def test_clock_and_velocity_records_distinguish_finite_spacing_on_the_same_loop():
    result = run()
    for record in result["objects"].values():
        assert 0 < record["measured_circulation_velocity"] < 1
        assert 0 < record["rest_normalized_rate"] < 1
        assert record["finite_spacing_rate_squared_defect"] > 0
        flows = [row["current"] for row in record["half_link_records"]]
        np.testing.assert_allclose(flows, flows[0], atol=2e-10, rtol=0)
