import numpy as np

from . import sector_return_lapse as lapse
from . import sector_return_redshift as sky


def test_mass_shell_bound_and_inconsistent_upper_endpoint():
    z = sky.shell_bound([1.0, 2.0, 4.0], [4.0, 0.0, 0.0], [4.0, 2.0, 1.0])
    assert z["failed_upper_constraints"] == 0
    assert np.isclose(z["minimum_log_rate_difference"] * sky.C**2, 3.0)
    assert (
        sky.shell_bound([1.0, 2.0], [4.0, 0.0], [4.0, 1.0])["failed_upper_constraints"]
        == 1
    )


def test_static_gain_respects_both_signs_and_blank_tile():
    assert lapse.gain_interval([1.0, 2.0], [2.0, 4.0]) == (1.0, 4.0)
    assert lapse.gain_interval([-2.0, -1.0], [-4.0, -2.0]) == (1.0, 4.0)
    assert lapse.gain_interval([0.0, 0.0], [1.0, 2.0]) == (1.0, 0.0)
    table = (np.array([1.0]), np.array([[0.5, 0.55]]), np.array([[-0.1, -0.05]]))
    ordinary, extra = lapse.clock_rates(np.array([0.0, 1.0]), table)
    np.testing.assert_array_equal(ordinary, 1.0)
    np.testing.assert_array_equal(extra, 0.0)
