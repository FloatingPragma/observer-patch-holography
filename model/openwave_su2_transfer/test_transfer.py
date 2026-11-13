import numpy as np
from .transfer import execute

def test_positive_transfer_and_complete_histories():
    r = execute()
    for row in r['rows']:
        assert min(row['complete_transfer_eigenvalues']) > 0
        assert np.linalg.eigvalsh(row['five_history_reflection_Gram'])[0] > -1e-13
        assert abs(sum(row['complete_13_outcome_plaquette_probabilities'])-1) < 1e-12

def test_common_clock_and_negative_reflection_control():
    r = execute()
    errors = [x['unit_time_semigroup_operator_error'] for x in r['rows']]
    assert errors[0] > 3*errors[1] > 9*errors[2]
    assert r['negative_transfer_control']['bond_reflection_quadratic_form_exact'] == '-1/2'
