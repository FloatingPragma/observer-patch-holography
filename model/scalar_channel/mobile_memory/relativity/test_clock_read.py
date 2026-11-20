import numpy as np
import pytest

from .clock_read import run, two_read_records


def test_two_read_clock_retains_first_outcomes_and_filter_failure():
    lower = np.array([0.5, 0.0], dtype=complex)
    upper = np.array([0.0, 0.0], dtype=complex)
    # Orthogonal tag arms with no overlap cannot support a phase read.
    with pytest.raises(ValueError, match="measurable phase"):
        two_read_records(lower, upper, observed_cell=0)
    upper[0] = 0.4 * np.exp(-0.71j)
    observed = two_read_records(lower, upper, observed_cell=0)
    np.testing.assert_allclose(observed["phase"], 0.71, atol=1e-12)
    np.testing.assert_allclose(observed["rejected_by_first_read"], [0.09, 0.09], atol=1e-12)
    assert not np.any(np.array(observed["joint_records"])[..., 1])


def test_common_incidence_clock_uses_its_measured_rest_rate():
    result = run()
    for record in result["objects"].values():
        np.testing.assert_allclose(record["squared_speed_from_records"], 1, atol=2e-10)
        assert 0 < record["rest_normalized_clock_rate"] < 1
        assert abs(record["equal_mass_control"]["observed_frequency"]) < 2e-10
        np.testing.assert_allclose(
            record["infinitesimal_rate_from_occupation"] ** 2
            + record["group_rate_in_declared_momentum"] ** 2, 1, atol=2e-10,
        )
