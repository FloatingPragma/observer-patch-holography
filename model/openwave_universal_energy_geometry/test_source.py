from fractions import Fraction as F
import numpy as np
from .source import exact_attraction, ground_interval, rest_block, autonomous_records

def test_every_native_pair_attracts_with_true_energy_reference():
    r = exact_attraction()
    assert r['all_147_positive_energy_pair_interactions_attractive']
    assert F(r['first_native_pair_inward_central_force_interval_exact'][0]) > 0
    for i, j in ((0, 0), (1, 7), (7, 3)):
        v = ground_interval(i, j, F(1, 3))
        e = np.linalg.eigvalsh(rest_block(i, j, F(1, 3)))[0]
        assert float(v.lo)-1e-13 <= e <= float(v.hi)+1e-13

def test_autonomous_energy_and_complete_readout_work():
    rows = autonomous_records()
    assert any(abs(x['geometry_and_matter_energy_rates'][0]) > 1e-4 for x in rows)
    assert all(abs(sum(x['complete_geometry_coordinate_probabilities'])-1) < 1e-12 for x in rows)
    assert any(abs(x['complete_geometry_measurement_work']) > 1e-5 for x in rows)
