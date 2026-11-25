import json
from pathlib import Path

import numpy as np

from . import invisible


def test_equal_twelve_ranges_have_a_signed_entropy_and_coin_response():
    data = invisible.run()
    assert data["range_rows"] == 12 and data["exact_range_kernel_dimension"] == 6
    assert data["true_kernel_correction_HS_upper"] < 4e-7
    assert data["channel_trace_error_upper"] < 1.91e-7
    assert data["finite_entropy_increase_nats_lower"] > 7e-13
    assert float(data["parity_blind_joint_YY_derivative"].strip("[]").split(",")[0]) > 1.53
    complete = data["fixed_marginal_tomography"]
    assert complete["parity_blind_products"] == ["XZ","YX","YY","YZ","ZY","ZZ"]
    assert complete["minimal_additional_scalar_means"] == 6
    assert complete["singular_value_lower"] > 9.7e-9
    assert complete["HS_reconstruction_error_per_l2_readout_error_upper"] < 1.04e8
    assert data["fixed_state_centers"]["denominator"] == 2**35
    held = json.loads(Path(invisible.__file__).with_suffix(".json").read_text())
    np.testing.assert_allclose(held["tangent_center_even_then_odd"],data["tangent_center_even_then_odd"],
                               rtol=0, atol=2e-7)
