import numpy as np

from .pair import measure


def test_modular_energy_comparison_retains_quiet_and_null_support():
    report, readings = measure()
    assert sum(report['cutoffs'][1]['ranks']) == 17
    for fit in report['cutoffs']:
        assert fit['max_outside_support_frobenius'] < 1e-11
        assert fit['port9_prediction'] < -0.01 and fit['energy_contrasts'][1] > 0.09
        assert fit['matches_held'] is False
    budget = report['reference_survey']['shots_per_setting']
    for name, counts in readings.items():
        if name == 'diagonal':
            assert counts.sum() == budget
        else:
            np.testing.assert_equal(counts.sum(axis=1), budget)
