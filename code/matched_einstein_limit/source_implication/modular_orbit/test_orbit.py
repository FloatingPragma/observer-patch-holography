import json

import numpy as np

from .controls import nominal_experiment
from .measure import verify


def test_true_station_bounds_cover_independent_replay():
    true = verify()
    read = nominal_experiment()
    for key, measure in [
        ("true_entropy_slope", "entropy_slope"),
        ("true_area_slope", "area_slope"),
        ("right_Z_slope", "right_Z_slope"),
    ]:
        low, high = json.loads(true[key])
        assert low < read[measure] < high
    assert json.loads(true["right_doubling_error_per_unit_response"])[0] > 0.25
    assert json.loads(true["right_forward_error_per_unit_response"])[0] > 1 / 3
    assert read["finite_step_commutator_norm"] < 2e-16
    assert read["minimum_preparation_weight"] > 0.003
    assert read["analyzer_reference_slope"] > 1e-5
    assert abs(read["normalized_analyzer_modular_derivative"] - 2 * np.pi) < 1e-8
    impulse = [e["measured_ratio"] for e in read["analyzer_endpoints"]]
    assert np.allclose(impulse[:3], [0.30683959, 1, 1.69312600], atol=2e-8, rtol=0)
    assert max(abs(impulse[2] - 2), abs(impulse[3] - 4)) > 1 / 3
    assert max(e["rail_operator_error"] for e in read["analyzer_endpoints"]) < 1e-13
