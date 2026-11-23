import json
from pathlib import Path

import numpy as np

from . import measure


def test_same_law_phase_readout_and_forgotten_outcome_control():
    archived = json.loads(Path(__file__).with_name('results.json').read_text())
    observed = measure.run()
    assert len(observed['records']) == len(archived['records']) == 48
    for row, stored in zip(observed['records'], archived['records'], strict=True):
        assert row['perimeter']['counts'] == stored['perimeter']['counts']
        assert all(sum(cohort) == measure.SHOTS for cohort in row['interior']['counts'])
        np.testing.assert_allclose(row['interior']['entropy_slope']['population'],
                                   stored['interior']['entropy_slope']['population'], atol=1e-11)
        assert row['forgotten_head_control']['held_law_tv'] < 1e-13
        if row['exposure'] == 0:
            assert row['perimeter']['held_law_tv'] < 1e-14
        elif row['exposure'] == 1:
            assert row['perimeter']['significance'] > 35
    assert observed['empirical_optimal_relative_residual'] > .99
    assert observed['optimal_relative_residual'] > .99
    assert observed['calibration_first_seam_coefficient'] < 0
